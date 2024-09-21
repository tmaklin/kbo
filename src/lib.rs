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

pub mod index;
pub mod map;

pub fn map(
    query_file: &String,
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<(usize, usize, usize, usize, char)> {
    let (k, threshold) = match sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), map::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4 as usize, 0.0000001 as f64))
	},
    };
    let translate_params = map::TranslateParams { k: k, threshold: threshold };

    // TODO handle multiple files and `input_list`
    let ms = index::query_sbwt(&query_file, &sbwt, &lcs);

    info!("Translating result...");
    let ms_fw = ms.iter().map(|x| x.0).collect::<Vec<usize>>();
    let ms_rev = ms.iter().map(|x| x.1).collect::<Vec<usize>>();
    let runs = (map::derandomize_ms(&ms_fw, &Some(translate_params.clone())),
		map::derandomize_ms(&ms_rev, &Some(translate_params.clone())));
    let aln = (map::translate_runs(&ms_fw, &runs.0, &Some(translate_params.clone())),
	       map::translate_runs(&ms_rev, &runs.1, &Some(translate_params)));
    let mut run_lengths: Vec<(usize, usize, usize, usize, char)> = map::run_lengths(&aln.0).iter().map(|x| (x.0, x.1, x.2, x.3, '+')).collect();
    let mut run_lengths_rev: Vec<(usize, usize, usize, usize, char)> = map::run_lengths(&aln.1).iter().map(|x| (x.0, x.1, x.2, x.3, '-')).collect();
    run_lengths.append(&mut run_lengths_rev);

    return run_lengths;
}
