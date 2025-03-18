// kbo: Spectral Burrows-Wheeler transform accelerated local alignment search
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
//! Gap filling using matching statistics and SBWT interval lookups.
use std::ops::Range;

/// Find the nearest unique context that overlaps the gap.
pub fn overlap_gap(
    noisy_ms: &[(usize, Range<usize>)],
    ref_seq: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    k: usize,
    threshold: usize,
    start_index: usize,
    end_index: usize,
) -> (usize, usize, Vec<u8>) {
    let mut kmer: Vec<u8> = Vec::with_capacity(k);
    let kmer_idx_start = (start_index + k - threshold).min(ref_seq.len() - threshold);
    let mut kmer_idx = kmer_idx_start;
    while kmer_idx > start_index + threshold {
        let sbwt_interval = &noisy_ms[kmer_idx].1;
        if sbwt_interval.end - sbwt_interval.start == 1 {
            sbwt.push_kmer_to_vec(sbwt_interval.start, &mut kmer);

            // Check that the overlapping parts of the
            // candidate k-mer match the reference sequence
            let right_matches_want = kmer_idx_start - (end_index - 1) - (kmer_idx_start - kmer_idx);
            let right_matches_got = count_right_overlaps(&kmer, ref_seq, end_index + right_matches_want);

            let ref_start_pos = if start_index > threshold { start_index - threshold } else { 0 };
            let left_matches_got = count_left_overlaps(&kmer, ref_seq, ref_start_pos);

            // If we find a k-mer very early the right overlap can be longer
            // than k, hence the .min(k) here
            let right_overlap_matches = right_matches_got == right_matches_want.min(k);

            let left_overlap_matches = left_matches_got == threshold;

            if right_overlap_matches && left_overlap_matches {
                break;
            } else {
                kmer.clear();
            }
        }
        kmer_idx -= 1;
    }
    (kmer_idx, kmer_idx_start, kmer)
}

/// Find the nearest unique context leftwards of a starting point.
pub fn nearest_unique_context(
    ms: &[(usize, Range<usize>)],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    search_start: usize,
    search_end: usize,
) -> (usize, Vec<u8>) {
    let k = sbwt.k();
    let mut kmer: Vec<u8> = Vec::with_capacity(k);
    let mut kmer_idx = search_start;

    while kmer_idx >= search_end {
        let sbwt_interval = &ms[kmer_idx].1;
        if sbwt_interval.end - sbwt_interval.start == 1 {
            sbwt.push_kmer_to_vec(sbwt_interval.start, &mut kmer);
            break;
        }
        kmer_idx -= 1;
    }

    (kmer_idx, kmer)
}

/// Count overlaps between a sequence and the last elements of a k-mer.
pub fn count_right_overlaps(
    kmer: &[u8],
    ref_seq: &[u8],
    ref_match_end: usize,
) -> usize {
    let mut kmer_pos = kmer.len() - 1;
    let mut ref_pos = ref_match_end - 1;
    let mut matches: usize = 0;
    while kmer_pos > 0 {
        if ref_seq[ref_pos] == kmer[kmer_pos] {
            matches += 1;
        } else {
            break;
        }
        kmer_pos -= 1;
        ref_pos -= 1;
    }
    matches
}

/// Count overlaps between a sequence and the first elements of a k-mer.
pub fn count_left_overlaps(
    kmer: &[u8],
    ref_seq: &[u8],
    ref_match_start: usize,
) -> usize {
    let mut kmer_pos = 0;
    let mut ref_pos = ref_match_start;
    let mut matches: usize = 0;
    while kmer_pos < kmer.len() {
        if ref_seq[ref_pos] == kmer[kmer_pos] {
            matches += 1;
        } else {
            break;
        }
        kmer_pos += 1;
        ref_pos += 1;
    }
    matches
}

/// Left extends a k-mer until the SBWT interval becomes non-unique
pub fn left_extend_kmer(
    kmer_start: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    max_extension_len: usize,
) -> Vec<u8> {
    // Found a candidate k-mer
    // try to left extend it as far as possible
    let mut left_extension_len = 0;

    let mut extend_count: usize = 0;
    let mut kmer = kmer_start.to_vec().clone();
    while extend_count < max_extension_len {
        let mut new_kmers: Vec<(Vec<u8>, Range<usize>)> = Vec::new();
        for c in [b'A',b'C',b'G',b'T'] {
            let mut new_kmer: Vec<u8> = Vec::new();
            new_kmer.push(c);
            new_kmer.append(&mut kmer.clone());
            for i in 0..(left_extension_len + 1) {
                new_kmer.pop();
            }
            let res = sbwt.search(&new_kmer).clone();
            if res.as_ref().is_some() {
                new_kmers.push((new_kmer, res.unwrap()));
            }
        }
        if new_kmers.len() == 1 && new_kmers[0].1.end - new_kmers[0].1.start == 1 {
            left_extension_len += 1;
            let mut kmer_extended = vec![new_kmers[0].0[0]];
            kmer_extended.append(&mut kmer.clone());
            kmer = kmer_extended;
        } else {
            break;
        }

        extend_count += 1;
    }
    kmer
}

/// Find and extend a nearest unique context until it overlaps the gap.
pub fn left_extend_over_gap(
    noisy_ms: &[(usize, Range<usize>)],
    ref_seq: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    overlap_req: usize,
    gap_start_index: usize,
    gap_end_index: usize,
) -> Vec<u8> {
    let k = sbwt.k();
    let search_start = (gap_end_index + k).min(ref_seq.len() - 1);
    let search_end = gap_end_index + overlap_req - 1;

    let mut kmer: Vec<u8> = Vec::with_capacity(k);
    let mut kmer_idx = search_start;
    while kmer_idx >= search_end {
        (kmer_idx, kmer) = nearest_unique_context(noisy_ms, sbwt, kmer_idx, search_end);
        if !kmer.is_empty() {
            // Check that the right overlapping parts of the candidate k-mer
            // match the reference sequence nucleotides
            let right_matches_want = search_start - (gap_end_index - 1) - (search_start - kmer_idx);
            let right_matches_got = count_right_overlaps(&kmer, ref_seq, gap_end_index + right_matches_want);

            // If we find a k-mer very early the right overlap can be longer
            // than k, hence the .min(k) here to check the right overlap match
            if right_matches_got == right_matches_want.min(k) {
                // Try to extend `kmer` left until it contains `overlap_req`
                // bases before the gap bases in `ref_seq`.
                let left_extend_length = overlap_req + (gap_end_index - gap_start_index) + right_matches_got - k;
                kmer = left_extend_kmer(&kmer, sbwt, left_extend_length);

                // We're done if the first `overlap_req` bases in `kmer` match
                // the `overlap_req` bases in `ref_seq` preceding `gap_start_index`
                let ref_start_pos = if gap_start_index > overlap_req { gap_start_index - overlap_req } else { 0 };
                let left_matches_got = count_left_overlaps(&kmer, ref_seq, ref_start_pos);

                if left_matches_got == overlap_req {
                    break;
                }
            }
            kmer.clear();
        }
        kmer_idx -= 1;
    }

    kmer
}
