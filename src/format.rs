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
//! Converting alignment representations into various output formats.
pub fn run_lengths(
    aln: &[char],
) -> Vec<(usize, usize, usize, usize)> {
    // Store run lengths as Vec<(start, end, matches, mismatches)>
    let mut encodings: Vec<(usize, usize, usize, usize)> = Vec::new();

    let mut i = 0;
    let mut match_start: bool = false;
    while i < aln.len() {
	match_start = (aln[i] != '-' && aln[i] != ' ') && !match_start;
	if match_start {
	    let start = i;
	    let mut matches: usize = 0;
	    while i < aln.len() && (aln[i] != '-' && aln[i] != ' ') {
		matches += (aln[i] == 'M' || aln[i] == 'R') as usize;
		i += 1;
	    }
	    encodings.push((start + 1, i, matches, i - start - matches));
	    match_start = false;
	} else {
	    i += 1;
	}
    }
    encodings
}

pub fn relative_to_ref(
    ref_seq: &[u8],
    alignment: &Vec<char>,
) -> Vec<u8> {
    ref_seq.iter().zip(alignment.iter()).map(|x| if *x.1 == 'M' || *x.1 == 'R' { *x.0 } else { b'-' }).collect()
}

pub fn run_lengths2(
    aln: &[char],
) -> Vec<(usize, usize, usize, usize)> {
    // Store run lengths as Vec<(start, end, matches, mismatches)>
    let mut encodings: Vec<(usize, usize, usize, usize)> = Vec::new();

    let mut i = 0;
    let mut match_start: bool = false;
    while i < aln.len() {
	match_start = (aln[i] != '-' && aln[i] != ' ') && !match_start;
	let mut discont: bool = false;
	if match_start {
	    let start = i;
	    let mut matches: usize = 0;
	    while i < aln.len() && (aln[i] != '-' && aln[i] != ' ' && !discont) {
		matches += (aln[i] == 'M' || aln[i] == 'R') as usize;
		i += 1;
		discont = aln[i - 1] == 'R';
	    }
	    encodings.push((start + 1, i, matches, i - start - matches));
	    match_start = false;
	} else {
	    i += 1;
	}
    }
    encodings
}

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {
    #[test]
    fn run_lengths() {
	let expected: Vec<(usize, usize, usize, usize)> = vec![(6,33,28,0),(82,207,126,0),(373,423,51,0),(488,512,25,0)];
	let input = vec!['-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M'];
	let got = super::run_lengths(&input);
	assert_eq!(got, expected);
    }
}
