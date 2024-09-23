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

use needletail::Sequence;
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

pub fn build_sbwt(
    infile: &String,
    params_in: &Option<SBWTParams>,
) -> (sbwt::SbwtIndex<sbwt::SubsetMatrix>, Option<sbwt::LcsArray>) {
    let params = params_in.clone().unwrap_or(SBWTParams::default());

    let temp_dir = params.temp_dir.unwrap_or(std::env::temp_dir());
    let algorithm = BitPackedKmerSorting::new()
	.mem_gb(params.mem_gb)
	.dedup_batches(false)
	.temp_dir(temp_dir.deref());

    let mut reader = needletail::parse_fastx_file(&infile.clone()).expect("valid path/file");

    let mut seqs = vec!();
    while let Some(rec) = reader.next()  {
	let seqrec = rec.expect("invalid_record");
	let seq = seqrec.normalize(true);
	seqs.push(seq.deref().to_owned());
    }

    let (sbwt, lcs) = SbwtIndexBuilder::new()
	.k(params.k)
	.n_threads(params.num_threads)
	.add_rev_comp(params.add_revcomp)
	.algorithm(algorithm)
	.build_lcs(true)
	.precalc_length(params.prefix_precalc)
	.run_from_vecs(&seqs);

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

pub fn load_sbwt(
    index_prefix: String,
) -> (sbwt::SbwtIndexVariant, sbwt::LcsArray) {
    let mut indexfile = index_prefix.clone();
    let mut lcsfile = indexfile.clone();
    indexfile.push_str(".sbwt");
    lcsfile.push_str(".lcs");

    // Read sbwt
    let mut index_reader = std::io::BufReader::new(std::fs::File::open(indexfile).unwrap());
    let sbwt = sbwt::load_sbwt_index_variant(&mut index_reader).unwrap();

    // Load the lcs array
    let lcs = match std::fs::File::open(&lcsfile) {
        Ok(f) => {
            let mut lcs_reader = std::io::BufReader::new(f);
            sbwt::LcsArray::load(&mut lcs_reader).unwrap()
        }
        Err(_) => {
            panic!("No LCS array found at {}", lcsfile);
        }
    };
    return (sbwt, lcs);
}

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
