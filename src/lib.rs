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

pub fn map(
    query_file: &String,
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> (Vec<char>, Vec<char>) {
    let (k, threshold) = match sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4 as usize, 0.0000001 as f64))
	},
    };
    // TODO handle multiple files and `input_list`

    let mut reader = needletail::parse_fastx_file(query_file).expect("valid path/file");
    let Some(rec) = reader.next() else { panic!("Invalid query {}", query_file); };
    let seqrec = rec.expect("invalid_record");

    let seq_fwd = seqrec.normalize(true);
    let ms_fwd = index::query_sbwt(seq_fwd.sequence(), &sbwt, &lcs);

    let seq_rev = seq_fwd.reverse_complement();
    let ms_rev = index::query_sbwt(seq_rev.sequence(), &sbwt, &lcs);

    info!("Translating result...");
    let runs = (derandomize::derandomize_ms_vec(&ms_fwd, k, threshold),
		derandomize::derandomize_ms_vec(&ms_rev, k, threshold));
    let aln = (translate::translate_ms_vec(&runs.0, k, threshold),
	       translate::translate_ms_vec(&runs.1, k, threshold));

    return aln;
}
