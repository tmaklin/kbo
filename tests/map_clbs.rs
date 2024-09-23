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
#[test]
fn map_nissle_against_clbs() {
    let (sbwt, lcs) = sablast::index::build_sbwt(&"tests/data/clbS.fna.gz".to_string(), &None);

    let expected = vec![(455, 967, 512, 1, '+')];
    let aln = sablast::map(&"tests/data/NZ_CP058217.1_clbS.fna.gz".to_string(), &sbwt::SbwtIndexVariant::SubsetMatrix(sbwt), &lcs.unwrap());

    let mut got: Vec<(usize, usize, usize, usize, char)> = sablast::format::run_lengths(&aln.0).iter().map(|x| (x.0, x.1, x.2, x.3, '+')).collect();
    let mut run_lengths_rev: Vec<(usize, usize, usize, usize, char)> = sablast::format::run_lengths(&aln.1).iter().map(|x| (x.0, x.1, x.2, x.3, '-')).collect();
    got.append(&mut run_lengths_rev);

    assert_eq!(got, expected);
}
