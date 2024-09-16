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
use sbwt::SbwtIndexVariant;

// Parameters for SBWT construction
#[derive(Clone)]
pub struct QueryParams {
    pub index_prefix: Option<String>,
}
// Defaults
impl Default for QueryParams {
    fn default() -> QueryParams {
        QueryParams {
	    index_prefix: None,
        }
    }
}

pub fn load_sbwt(
    params_in: &Option<QueryParams>,
) -> (sbwt::SbwtIndexVariant, sbwt::LcsArray) {
    let params = params_in.clone().unwrap_or(QueryParams::default());

    let mut indexfile = params.index_prefix.clone().unwrap();
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
    query_file: &String,
    index: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<(usize, std::ops::Range<usize>)> {
    match index {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
	    let mut reader = jseqio::reader::DynamicFastXReader::from_file(&query_file).unwrap();
	    let streaming_index = sbwt::StreamingIndex::new(sbwt, lcs);

	    // TODO handle input with multiple sequences
	    // implement as querying 1 record at a time
	    let rec = reader.read_next().unwrap();
	    let ms = streaming_index.matching_statistics(rec.unwrap().seq);
	    return ms;
	},
    };
}
