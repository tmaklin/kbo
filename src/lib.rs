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

pub mod build;
pub mod map;

pub fn map(
    query_file: &String,
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<(usize, usize, usize, usize)> {
    let translate_params = map::TranslateParams {
	k: match sbwt {
	    SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
		sbwt.k()
	    }
	},
	threshold: match sbwt {
	    SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
		// 0.0000001 should maybe be prop. to 1/(number of unique k-mers in query)
		map::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4 as usize, 0.0000001 as f64)
	    }
	},
    };

    // TODO handle multiple files and `input_list`
    let ms = map::query_sbwt(&query_file, &sbwt, &lcs);

    info!("Translating result...");
    let ms_vec = ms.iter().map(|x| x.0).collect::<Vec<usize>>();
    let runs = map::derandomize_ms(&ms_vec, &Some(translate_params.clone()));
    let aln = map::translate_runs(&ms_vec, &runs, &Some(translate_params));
    let run_lengths = map::run_lengths(&aln);

    return run_lengths;
}
