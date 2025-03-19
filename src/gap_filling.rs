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

/// Count overlaps between a sequence and the last elements of a k-mer.
fn count_right_overlaps(
    kmer: &[u8],
    ref_seq: &[u8],
    ref_match_end: usize,
) -> usize {
    assert!(!kmer.is_empty());
    assert!(!ref_seq.is_empty());
    assert!(ref_seq.len() >= ref_match_end);

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
fn count_left_overlaps(
    kmer: &[u8],
    ref_seq: &[u8],
    ref_match_start: usize,
) -> usize {
    assert!(!kmer.is_empty());
    assert!(!ref_seq.is_empty());
    assert!(ref_seq.len() > ref_match_start);

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

/// Find the nearest unique context leftwards of a starting point.
pub fn nearest_unique_context(
    ms: &[(usize, Range<usize>)],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    search_start: usize,
    search_end: usize,
) -> (usize, Vec<u8>) {
    let k = sbwt.k();
    assert!(k > 0);
    assert!(!ms.is_empty());
    assert!(search_start >= search_end);
    assert!(search_start < ms.len());

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

/// Left extends a k-mer until the SBWT interval becomes non-unique
pub fn left_extend_kmer(
    kmer_start: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    max_extension_len: usize,
) -> Vec<u8> {
    assert!(!kmer_start.is_empty());

    let mut left_extension_len = 0;
    let mut kmer = kmer_start.to_vec();
    while left_extension_len < max_extension_len {
        let new_kmers: Vec<(Vec<u8>, Range<usize>)> = sbwt.alphabet().iter().filter_map(|c| {
            let new_kmer: Vec<u8> = [&[*c], &kmer[0..(kmer.len() - (left_extension_len + 1))]].concat();
            let res = sbwt.search(&new_kmer);
            if res.as_ref().is_some() {
                Some((new_kmer, res.unwrap()))
            } else {
                None
            }
        }).collect();
        if new_kmers.len() == 1 && new_kmers[0].1.end - new_kmers[0].1.start == 1 {
            kmer = [&[new_kmers[0].0[0]], kmer.as_slice()].concat();
        } else {
            break;
        }
        left_extension_len += 1;
    }
    kmer
}

/// Find and extend a nearest unique context until it overlaps the gap.
pub fn left_extend_over_gap(
    noisy_ms: &[(usize, Range<usize>)],
    ref_seq: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    left_overlap_req: usize,
    right_overlap_req: usize,
    gap_index: Range<usize>,
    search_radius: usize,
) -> Vec<u8> {
    let k = sbwt.k();
    assert!(k > 0);
    assert!(noisy_ms.len() == ref_seq.len());
    assert!(left_overlap_req <= gap_index.start);
    assert!(right_overlap_req <= ref_seq.len() - gap_index.end);
    assert!(gap_index.end > gap_index.start);
    assert!(gap_index.end < noisy_ms.len());

    let search_start = (gap_index.end + search_radius).min(ref_seq.len() - 1);
    let search_end = gap_index.end + right_overlap_req;

    let mut kmer: Vec<u8> = Vec::with_capacity(k);
    let mut kmer_idx = search_start;
    while kmer_idx >= search_end {
        (kmer_idx, kmer) = nearest_unique_context(noisy_ms, sbwt, kmer_idx, search_end);
        if !kmer.is_empty() {
            // Check that the right overlapping parts of the candidate k-mer
            // match the reference sequence nucleotides
            let right_matches_want = search_start - (gap_index.end - 1) - (search_start - kmer_idx);
            let right_matches_got = count_right_overlaps(&kmer, ref_seq, gap_index.end + right_matches_want);

            // We're done if the first `left_overlap_req` bases in `kmer` match
            // the `left_overlap_req` bases in `ref_seq` preceding `gap_index.start`
            let ref_start_pos = if gap_index.start > left_overlap_req { gap_index.start - left_overlap_req } else { 0 };
            let left_matches_got = count_left_overlaps(&kmer, ref_seq, ref_start_pos);

            // There's no point in extending if the k-mer already overlaps the
            // gap to the left but contains no matches
            let should_extend = kmer.len() < left_overlap_req + (gap_index.end - gap_index.start) + right_matches_got;

            if right_matches_got >= right_matches_want.min(k) && left_matches_got >= left_overlap_req {
                let start = left_matches_got - left_overlap_req;
                let end = kmer.len() - (right_matches_got - right_overlap_req);
                kmer = kmer[start..end].to_vec();
                break;
            } else if should_extend && right_matches_got >= right_matches_want.min(k) && left_matches_got < left_overlap_req {
                // If we find a k-mer very early the right overlap can be longer
                // than k, hence the .min(k) above to check the right overlap match

                // Try to extend `kmer` left until it contains `left_overlap_req`
                // bases before the gap bases in `ref_seq`.
                let left_extend_length = left_overlap_req + (gap_index.end - gap_index.start) + right_matches_got - k;
                kmer = left_extend_kmer(&kmer, sbwt, left_extend_length);

                if count_left_overlaps(&kmer, ref_seq, ref_start_pos) >= left_overlap_req {
                    let start = count_left_overlaps(&kmer, ref_seq, ref_start_pos) - left_overlap_req;
                    let end = kmer.len() - (right_matches_got - right_overlap_req);
                    kmer = kmer[start..end].to_vec();
                    break;
                }
            }
            kmer.clear();
        }
        kmer_idx -= 1;
    }

    kmer
}
