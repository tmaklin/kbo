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
/// from the index after it has been built; use [matches] with the TODO
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

/// Matches a query fasta or fastq file against an SBWT index.
///
/// Queries the sequence data in `query_seq` against the SBWT index
/// `sbwt` and its LCS array `lcs` using [index::query_sbwt]. Then,
/// derandomizes the resulting _k_-bounded matching statistics vector
/// using [derandomize::derandomize_ms_vec] and translates the
/// matching statistics to a character representation of the alignment
/// using [translate::translate_ms_vec].
///
/// Returns a vector containing the character representation of the
/// alignment.
///
/// Panics if the query file is not readable or if it's not a valid
/// FASTX file.
///
/// # Output format
/// See the documentation for [translate].
///
/// # Example
/// ```rust
/// use sablast::build;
/// use sablast::matches;
/// use sablast::index::BuildOpts;
///
/// let reference = vec!["tests/data/clbS.fna.gz".to_string()];
/// let (sbwt, lcs) = build(&reference, BuildOpts::default());
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let ms_vectors = matches(&query, &sbwt, &lcs);
/// ```
///
pub fn matches(
    query_seq: &[u8],
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<char> {
    let (k, threshold) = match sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, 0.0000001_f64))
	},
    };

    let noisy_ms = index::query_sbwt(query_seq, sbwt, lcs);
    let derand_ms = derandomize::derandomize_ms_vec(&noisy_ms, k, threshold);

    translate::translate_ms_vec(&derand_ms, k, threshold)
}

/// Finds the _k_-mers from an SBWT index in a query fasta or fastq file.
///
/// Aligns the sequence data in `query_seq` against the SBWT index
/// `sbwt` and its LCS array `lcs` using [matches]. Then uses
/// [format::run_lengths] to extract the local alignments from the
/// matching statistics.
///
/// Returns a vector of tuples, where each element represents a local
/// alignment block and contains the following values:
/// 1. Start of local alignment block in query (1-based indexing).
/// 2. End of local alignment block in query.
/// 3. Number of matches in the block.
/// 4. Number of mismatches and 1-character insertions in the block.
///
/// # Examples
/// ```rust
/// use sablast::build;
/// use sablast::find;
/// use sablast::index::BuildOpts;
///
/// let reference = vec!["tests/data/clbS.fna.gz".to_string()];
/// let (sbwt, lcs) = build(&reference, BuildOpts::default());
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let local_alignments = find(&query, &sbwt, &lcs);
/// ```
///
pub fn find(
    query_seq: &[u8],
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<(usize, usize, usize, usize)> {
    let aln = matches(query_seq, sbwt, lcs);
    format::run_lengths(&aln)
}
