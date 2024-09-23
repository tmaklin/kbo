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
use std::ffi::OsString;
use std::io::Write;
use std::path::PathBuf;

use needletail::Sequence;
use sbwt::BitPackedKmerSorting;
use sbwt::SbwtIndexBuilder;
use sbwt::SbwtIndexVariant;

// Parameters for SBWT construction
#[derive(Clone)]
pub struct BuildOpts {
    pub k: usize,
    pub add_revcomp: bool,
    pub num_threads: usize,
    pub mem_gb: usize,
    pub temp_dir: Option<String>,
    pub prefix_precalc: usize,
}
// Defaults
impl Default for BuildOpts {
    fn default() -> BuildOpts {
        BuildOpts {
	    k: 31,
	    add_revcomp: false,
	    num_threads: 1,
	    mem_gb: 4,
	    prefix_precalc: 8,
	    temp_dir: None,
        }
    }
}

struct FastxStreamer {
    inner: Box<dyn needletail::parser::FastxReader>,
    record: Vec<u8>
}

impl sbwt::SeqStream for FastxStreamer {
    fn stream_next(&mut self) -> Option<&[u8]> {
	let rec = self.inner.next();
	match rec {
	    Some(Ok(seqrec)) => {
		// Remove newlines and non IUPAC characters
		let normalized = seqrec.normalize(true);
		self.record = normalized.as_ref().to_vec();
		Some(&self.record)
	    },
	    _ => None,
	}
    }
}

/// Builds an SBWT index and its LCS array from a fasta or fastq file.
///
/// Streams all valid DNA sequences from `infile` to the SBWT API
/// calls to build the SBWT index and LCS array. Use the [BuildOpts]
/// argument `build_options` to control the options and resources
/// passed to the index builder.
///
/// Returns a tuple containing the SBWT index and the LCS array.
///
/// Requires write access to some temporary directory. Path can be set
/// using temp_dir in BuildOpts; defaults to $TMPDIR on Unix if not set.
///
/// # Examples
/// TODO Add examples to build_sbwt documentation.
///
pub fn build_sbwt_from_file(
    infile: &str,
    build_options: &Option<BuildOpts>,
) -> (sbwt::SbwtIndex<sbwt::SubsetMatrix>, Option<sbwt::LcsArray>) {
    // Get temp dir path from build_options, otherwise use whatever std::env::temp_dir() returns
    let temp_dir = build_options.as_ref().unwrap().temp_dir.clone().unwrap_or(std::env::temp_dir().to_str().unwrap().to_string());

    let algorithm = BitPackedKmerSorting::new()
	.mem_gb(build_options.as_ref().unwrap().mem_gb)
	.dedup_batches(false)
	.temp_dir(PathBuf::from(OsString::from(temp_dir)).as_path());

    let reader = FastxStreamer{inner: needletail::parse_fastx_file(infile).expect("valid path/file"), record: Vec::new()};

    let (sbwt, lcs) = SbwtIndexBuilder::new()
	.k(build_options.as_ref().unwrap().k)
	.n_threads(build_options.as_ref().unwrap().num_threads)
	.add_rev_comp(build_options.as_ref().unwrap().add_revcomp)
	.algorithm(algorithm)
	.build_lcs(true)
	.precalc_length(build_options.as_ref().unwrap().prefix_precalc)
	.run(reader);

    (sbwt, lcs)
}

pub fn serialize_sbwt(
    outfile_prefix: &str,
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: &sbwt::LcsArray,
) {
    let sbwt_outfile = format!("{}.sbwt", outfile_prefix);
    let lcs_outfile = format!("{}.lcs", outfile_prefix);

    // Write sbwt
    let sbwt_conn = std::fs::File::create(&sbwt_outfile).unwrap_or_else(|_| panic!("Expected write access to {}", sbwt_outfile));
    let mut sbwt_out = std::io::BufWriter::new(sbwt_conn);
    sbwt_out.write_all(&(b"SubsetMatrix".len() as u64).to_le_bytes()).expect("Serialized SBWT header part 1.");
    sbwt_out.write_all(b"SubsetMatrix").expect("Serialized SBWT header part 2.");
    sbwt.serialize(&mut sbwt_out).expect("Serialized SBWT index.");

    // Write lcs array
    let lcs_conn = std::fs::File::create(&lcs_outfile).unwrap_or_else(|_| panic!("Expected write access to {}", lcs_outfile));
    let mut lcs_out = std::io::BufWriter::new(lcs_conn);
    lcs.serialize(&mut lcs_out).expect("Serialized LCS array.");
}

/// Loads a prebuilt SBWT index and its LCS array from disk.
///
/// Reads the SBWT index stored at `index_prefix` + ".sbwt" and the
/// LCS array at `index_prefix` + ".lcs".
///
/// Returns a tuple containing the SBWT index variant and the LCS
/// array.
///
/// Panics if the SBWT or the LCS file are not readable with
/// std::fs::File::open.
///
/// # Examples
/// TODO Add examples to load_sbwt documentation.
///
pub fn load_sbwt(
    index_prefix: &str,
) -> (sbwt::SbwtIndexVariant, sbwt::LcsArray) {
    let indexfile = format!("{}.sbwt", index_prefix);
    let lcsfile = format!("{}.lcs", index_prefix);

    // Load sbwt
    let sbwt_conn = std::fs::File::open(&indexfile).unwrap_or_else(|_| panic!("Expected SBWT at {}", indexfile));
    let mut index_reader = std::io::BufReader::new(sbwt_conn);
    let sbwt = sbwt::load_sbwt_index_variant(&mut index_reader).unwrap();

    // Load the lcs array
    let lcs_conn = std::fs::File::open(&lcsfile).unwrap_or_else(|_| panic!("Expected LCS array at {}", lcsfile));
    let mut lcs_reader = std::io::BufReader::new(lcs_conn);
    let lcs = sbwt::LcsArray::load(&mut lcs_reader).unwrap();

    (sbwt, lcs)
}

/// Queries an SBWT index for the _k_-bounded matching statistics.
///
/// Matches the _k_-mers in `query` against the SBWT index `index` and
/// its longest common suffix array `lcs`.
///
/// Returns a vector containing the _k_-bounded matching statistic at
/// the position of each element in the query.
///
/// # Examples
/// TODO Add examples to query_sbwt documentation
///
pub fn query_sbwt(
    query: &[u8],
    index: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<usize> {
    let ms = match index {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
	    let streaming_index = sbwt::StreamingIndex::new(sbwt, lcs);
	    streaming_index.matching_statistics(query)
	},
    };
    ms.iter().map(|x| x.0).collect()
}
