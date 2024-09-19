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
pub struct TranslateParams {
    pub k: usize,
    pub threshold: usize,
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

fn log_rm_max_cdf(
    t: usize,
    alphabet_size: usize,
    n_kmers: usize,
) -> f64 {
    n_kmers as f64 * (- ((1.0/(alphabet_size as f64)) as f64).powi(t as i32 + 1)).ln_1p()
}

pub fn random_match_threshold(
    k: usize,
    n_kmers: usize,
    alphabet_size: usize,
    max_error_prob: f64,
) -> usize {
    for i in 1..k {
	if 1.0 - log_rm_max_cdf(i, alphabet_size, n_kmers).exp() < max_error_prob {
	    return i;
	}
    }
    return k;
}

fn ms_to_run(
    curr: usize,
    next: usize,
    next_run: i64,
    threshold: usize,
    k: usize,
) -> i64 {
    let run: i64 = if curr == k && next == k {
	k as i64
    } else if curr == k && next_run == 1 {
	k as i64
    } else if curr == k && next_run < 0 {
	k as i64
    } else if curr < threshold {
	next_run - 1
    } else if curr > threshold && next_run <= 0 {
	curr as i64
    } else if curr > threshold && next_run == 1 {
	curr as i64
    } else if curr > threshold && next_run < curr as i64 {
	curr as i64
    } else {
	next_run - 1
    };

    return run;
}

fn run_to_aln(
    runs: &Vec<i64>,
    curr_ms: usize,
    threshold: usize,
    k: usize,
    res: &mut Vec<char>,
    pos: &mut usize,
) {
    let prev: i64 = runs[*pos - 1];
    let curr: i64 = runs[*pos];
    let next: i64 = runs[*pos + 1];

    if curr == k as i64 && next == k as i64 {
	res[*pos] = 'M';
    } else if curr > threshold as i64 && (next > 0 && next < threshold as i64) {
	res[*pos] = 'R';
	res[*pos + 1] = 'R';
    } else if next == 1 && curr == curr_ms as i64 {
	res[*pos] = 'M';
    } else if curr > threshold as i64 {
	res[*pos] = 'M';
    } else if curr == next - 1 && curr > 0 {
	res[*pos] = 'M';
    } else if curr == 0 && next == 1 && prev > 0 {
	res[*pos] = 'X';
	res[*pos - 1] = 'M';
    } else if curr == 0 && next == 1 && prev == -1 {
	let mut next_gap: usize = pos.clone();
	let mut gap_len: usize = 0;
	while runs[next_gap - 1] < 0 && next_gap > 1 {
	    gap_len += 1;
	    next_gap -= 1;
	}
	// TODO determine based on k and threshold: k - 2*threshold ?
	while *pos < *pos + gap_len && *pos < runs.len() {
	    res[*pos] = if gap_len > 29 { '-' } else { 'I' };
	    *pos += 1;
	}
    } else {
	res[*pos] = ' ';
    };
}

pub fn derandomize_ms(
    ms: &Vec<usize>,
    params_in: &Option<TranslateParams>
) -> Vec<i64> {
    let params = params_in.clone().unwrap();
    let len = ms.len();

    let mut runs: Vec<i64> = vec![0; len];

    // Traverse the matching statistics in reverse
    runs[len - 1] = ms[len - 1] as i64;
    for i in 2..len {
	runs[len - i] = ms_to_run(ms[len - i], ms[len - i + 1], runs[len - i + 1], params.threshold, params.k);
    }

    return runs;
}


pub fn translate_runs(
    ms: &Vec<usize>,
    runs: &Vec<i64>,
    params_in: &Option<TranslateParams>,
) -> Vec<char> {
    let params = params_in.clone().unwrap();
    let len = runs.len();
    let mut aln = vec![' '; len];

    // Traverse the runs
    for mut i in 3..(len - 1) {
	run_to_aln(&runs, ms[i], params.threshold, params.k, &mut aln, &mut i);
    }

    return aln;
}

pub fn run_lengths(
    aln: &Vec<char>,
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
    return encodings;
}
