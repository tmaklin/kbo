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
use sbwt::SbwtIndexVariant;

pub mod format;
pub mod index;
pub mod translate;

pub fn map(
    query_file: &String,
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> (Vec<char>, Vec<char>) {
    let (k, threshold) = match sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), translate::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4 as usize, 0.0000001 as f64))
	},
    };
    // TODO handle multiple files and `input_list`
    let ms = index::query_sbwt(&query_file, &sbwt, &lcs);

    info!("Translating result...");
    let ms_fw = ms.iter().map(|x| x.0).collect::<Vec<usize>>();
    let ms_rev = ms.iter().map(|x| x.1).collect::<Vec<usize>>();
    let runs = (translate::derandomize_ms_vec(&ms_fw, k, threshold),
		translate::derandomize_ms_vec(&ms_rev, k, threshold));
    let aln = (translate::translate_ms_vec(&ms_fw, &runs.0, k, threshold),
	       translate::translate_ms_vec(&ms_rev, &runs.1, k, threshold));

    return aln;
}
