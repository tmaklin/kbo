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
use std::ops::Deref;
use std::path::PathBuf;

use sbwt::BitPackedKmerSorting;
use sbwt::SbwtIndexBuilder;
use sbwt::SbwtIndexVariant;

// Parameters for SBWT construction
#[derive(Clone)]
pub struct SBWTParams {
    pub k: usize,
    pub add_revcomp: bool,
    pub num_threads: usize,
    pub mem_gb: usize,
    pub prefix_precalc: usize,
    pub temp_dir: Option<PathBuf>,
    pub index_prefix: Option<String>,
}
// Defaults
impl Default for SBWTParams {
    fn default() -> SBWTParams {
        SBWTParams {
	    k: 31,
	    add_revcomp: false,
	    num_threads: 1,
	    mem_gb: 4,
	    prefix_precalc: 8,
	    temp_dir: None,
	    index_prefix: None,
        }
    }
}

// SBWT needs the reader like this
struct MySeqReader {
    inner: jseqio::reader::DynamicFastXReader,
}
impl sbwt::SeqStream for MySeqReader {
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|rec| rec.seq)
    }
}
pub fn build_sbwt(
    infile: &String,
    params_in: &Option<SBWTParams>,
) -> (sbwt::SbwtIndex<sbwt::SubsetMatrix>, Option<sbwt::LcsArray>) {
    let params = params_in.clone().unwrap_or(SBWTParams::default());

    let reader = MySeqReader{inner: jseqio::reader::DynamicFastXReader::from_file(infile).unwrap()};

    let temp_dir = params.temp_dir.unwrap_or(std::env::temp_dir());
    let algorithm = BitPackedKmerSorting::new()
	.mem_gb(params.mem_gb)
	.dedup_batches(false)
	.temp_dir(temp_dir.deref());

    let (sbwt, lcs) = SbwtIndexBuilder::new()
	.k(params.k)
	.n_threads(params.num_threads)
	.add_rev_comp(params.add_revcomp)
	.algorithm(algorithm)
	.build_lcs(true)
	.precalc_length(params.prefix_precalc)
	.run(reader);

    return (sbwt, lcs);
}

pub fn serialize_sbwt(
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: &Option<sbwt::LcsArray>,
    params_in: &Option<SBWTParams>,
) {
    let params = params_in.clone().unwrap_or(SBWTParams::default());

    let mut sbwt_outfile = params.index_prefix.clone().unwrap_or("sbwt".to_string());
    sbwt_outfile.push_str(".sbwt");
    let mut sbwt_out = std::io::BufWriter::new(std::fs::File::create(sbwt_outfile).unwrap());

    sbwt.n_kmers();
    sbwt::write_sbwt_index_variant(&SbwtIndexVariant::SubsetMatrix(sbwt), &mut sbwt_out).unwrap();

    if let Some(lcs) = lcs{
        let mut lcs_outfile = params.index_prefix.clone().unwrap_or("sbwt".to_string());
        lcs_outfile.push_str(".lcs");
        let mut lcs_out = std::io::BufWriter::new(std::fs::File::create(&lcs_outfile).unwrap());
        lcs.serialize(&mut lcs_out).unwrap();
    }
}
