// sablast: Spectral Burrows-Wheeler transform accelerated local alignment search
//
// Copyright 2024 Tommi Mäklin [tommi@maklin.fi].

// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.

// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
use log::info;
use needletail::Sequence;
use sbwt::SbwtIndexVariant;

pub mod derandomize;
pub mod format;
pub mod index;
pub mod translate;

/// Builds an SBWT index from some fasta or fastq files.
///
/// Reads all sequence data in `seq_files` and builds an SBWT index
/// with the parameters and resources specified in `build_opts` (see
/// [index::BuildOpts] for details).
///
/// All files and sequence data in `seq_files` are merged into the
/// same index. It is not possible extract the individual sequences
/// from the index after it has been built; use [map] with the TODO
/// options instead if you need to know which reference sequences the
/// alignments are for.
///
/// TODO Describe map syntax in lib.rs documentation.
///
/// Returns a tuple containing the built
/// [sbwt::SbwtIndexVariant](https://docs.rs/sbwt/latest/sbwt/enum.SbwtIndexVariant.html)
/// and
/// [sbwt::LcsArray](https://docs.rs/sbwt/latest/sbwt/struct.LcsArray.html).
///
/// Panics if a file in `seq_files` is not readable or a valid FASTX
/// file.
///
/// # Input format detection
/// The sequence data is read using
/// [needletail::parser::parse_fastx_file](https://docs.rs/needletail/latest/needletail/parser/fn.parse_fastx_file.html).
///
/// Input file format (fasta or fastq) is detected automatically and
/// the files may be compressed in a
/// [DEFLATE-based](https://en.wikipedia.org/wiki/Deflate) format (.gz
/// files).
///
/// [Bzip2](https://sourceware.org/bzip2/) and
/// [liblzma](https://tukaani.org/xz/) compression (.bz2 and .xz
/// files) can be enabled using the needletail features field in
/// sablast Cargo.toml if compiling from source.
///
/// # Examples
/// ```rust
/// use sablast::build;
/// use sablast::index::BuildOpts;
///
/// let inputs = vec!["tests/data/clbS.fna.gz".to_string(), "tests/data/NZ_CP058217.1_clbS.fna.gz".to_string()];
///
/// let (sbwt_index, lcs_array) = build(&inputs, BuildOpts::default());
/// ```
///
pub fn build(
    seq_files: &[String],
    build_opts: index::BuildOpts,
) -> (sbwt::SbwtIndexVariant, sbwt::LcsArray) {
    let mut seq_data: Vec<Vec<u8>> = Vec::new();
    seq_files.iter().for_each(|file| {
	let mut reader = needletail::parse_fastx_file(file).unwrap_or_else(|_| panic!("Expected valid fastX file at {}", file));
	loop {
	    let rec = reader.next();
	    match rec {
		Some(Ok(seqrec)) => {
		    seq_data.push(seqrec.normalize(true).as_ref().to_vec());
		},
		_ => break
	    }
	}
    });

    index::build_sbwt_from_vecs(&seq_data, &Some(build_opts))
}

pub fn map(
    query_file: &String,
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> (Vec<char>, Vec<char>) {
    let (k, threshold) = match sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, 0.0000001_f64))
	},
    };
    // TODO handle multiple files and `input_list`

    let mut reader = needletail::parse_fastx_file(query_file).expect("valid path/file");
    let Some(rec) = reader.next() else { panic!("Invalid query {}", query_file); };
    let seqrec = rec.expect("invalid_record");

    let seq_fwd = seqrec.normalize(true);
    let ms_fwd = index::query_sbwt(seq_fwd.sequence(), sbwt, lcs);

    let seq_rev = seq_fwd.reverse_complement();
    let ms_rev = index::query_sbwt(seq_rev.sequence(), sbwt, lcs);

    info!("Translating result...");
    let runs = (derandomize::derandomize_ms_vec(&ms_fwd, k, threshold),
		derandomize::derandomize_ms_vec(&ms_rev, k, threshold));

    (translate::translate_ms_vec(&runs.0, k, threshold),
     translate::translate_ms_vec(&runs.1, k, threshold))
}

pub fn find(
    query_file: &String,
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<(usize, usize, char, usize, usize)> {
    let aln = map(query_file, &sbwt, &lcs);

    let mut run_lengths: Vec<(usize, usize, char, usize, usize)> = format::run_lengths(&aln.0).iter().map(|x| (x.0, x.1, '+', x.2 + x.3, x.3)).collect();
    let mut run_lengths_rev: Vec<(usize, usize, char, usize, usize)> = format::run_lengths(&aln.1).iter().map(|x| (x.0, x.1, '+', x.2 + x.3, x.3)).collect();

    run_lengths.append(&mut run_lengths_rev);
    run_lengths.sort_by_key(|x| x.0);

    run_lengths
}
