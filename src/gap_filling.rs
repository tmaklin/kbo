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
///
/// Queries the `sbwt` for _k_-mers in `search_range` from `ms`, starting from
/// the end of `search_range`, and checks whether SBWT interval of a queried
/// _k_-mer is of length 1.
///
/// Returns the position and contents of the rightmost _k_-mer in `search_range`
/// that has a SBWT interval of length 1. If a _k_-mer that satisfies these
/// conditions cannot be found, returns the position the search terminated at
/// and an empty _k_-mer.
///
/// The returned _k_-mer is stored as an ASCII-encoded vector.
///
/// Note that this function does not check for dummy characters ('$') in the
/// return value.
///
/// # Examples
///
/// This example shows how nearest_unique_context is used by kbo to
/// determine the nucleotides of the query sequence in a section that does not
/// match the reference. Note the first character in the output sequence.
///
/// ```rust
/// use kbo::index::BuildOpts;
/// use kbo::build;
/// use sbwt::SbwtIndexVariant;
/// use kbo::index::query_sbwt;
/// use kbo::gap_filling::nearest_unique_context;
///
/// // Parameters       : k = 7
/// //
/// // Ref sequence     : T,T,G,A, T, C,T G,G,C,T,G,C,T,G,A,G,A,G,C,T,G
/// // Query sequence   : T,T,G,A, A, C,A,G,G,C,T,G,C,T,C,A,G,A,G,C,T,G
/// //
/// // Search range     :                 - - - - - - - -
/// // Output           :               A,G,G,C,T,G,C
///
/// let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'T',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
/// let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'C',b'T',b'G',b'G',b'C',b'T',b'G',b'C',b'T',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
///
/// let mut opts = BuildOpts::default();
/// opts.k = 7;
/// opts.build_select = true;
/// opts.add_revcomp = false;
/// let (sbwt, lcs) = build(&[query], opts);
///
/// let ms = query_sbwt(&reference, &sbwt, &lcs);
///
/// let context = match sbwt {
///     SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
///         nearest_unique_context(&ms, &sbwt, 8..14)
///     },
/// };
///
/// # let expected = (12, vec![b'A',b'G',b'G',b'C',b'T',b'G',b'C']);
/// # assert_eq!(context, expected);
/// ```
///
pub fn nearest_unique_context(
    ms: &[(usize, Range<usize>)],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    search_range: Range<usize>,
) -> (usize, Vec<u8>) {
    let k = sbwt.k();
    assert!(k > 0);
    assert!(!ms.is_empty());
    assert!(search_range.end >= search_range.start);
    assert!(search_range.end < ms.len());

    let mut kmer: Vec<u8> = Vec::with_capacity(k);
    let mut kmer_idx = search_range.end;

    while kmer_idx >= search_range.start {
        let sbwt_interval = &ms[kmer_idx].1;
        if sbwt_interval.end - sbwt_interval.start == 1 {
            sbwt.push_kmer_to_vec(sbwt_interval.start, &mut kmer);
            break;
        }
        kmer_idx -= 1;
    }

    (kmer_idx, kmer)
}

/// Left extends a k-mer until the SBWT interval becomes non-unique.
///
/// Starting from the sequence in `kmer_start`, attempts to extend it to the
/// left by querying `sbwt` for patterns of the form cP, where P contains the
/// sequence in `kmer_start` and c is a character from the alphabet used in
/// `sbwt`.
///
/// If `sbwt` contains a single pattern cP, and no others, and the SBWT
/// interval of cP is of length 1, `kmer_start` is extended to the left with c.
///
/// This process repeats until the condition is no longer true, or `kmer_start`
/// has been extended `max_extension_len` times.
///
/// Returns the ASCII-encoded extended kmer.
///
/// # Examples
///
/// ```rust
/// use kbo::index::BuildOpts;
/// use kbo::build;
/// use sbwt::SbwtIndexVariant;
/// use kbo::index::query_sbwt;
/// use kbo::gap_filling::left_extend_kmer;
///
/// // Parameters       : k = 7
/// //
/// // Sequence         : T,T,G,A,A,C,A,G,G,C,T,G,C,C,G,T,A,A,C,A,G,G
/// //
/// // Starting k-mer   :             A,G,G,C,T,G,C
/// // Search range     :   - - - - - - - - - - - -
/// // Output           :       A,A,C,A,G,G,C,T,G,C
///
/// let sequence: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C',b'C',b'G',b'T',b'A',b'A',b'C',b'A',b'G',b'G'];
///
/// let mut opts = BuildOpts::default();
/// opts.k = 7;
/// opts.build_select = true;
/// opts.add_revcomp = false;
/// let (sbwt, lcs) = build(&[sequence], opts);
///
/// let kmer = vec![b'A',b'G',b'G',b'C',b'T',b'G',b'C'];
/// let max_extension_len = 5;
///
/// let extended = match sbwt {
///     SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
///         left_extend_kmer(&kmer, &sbwt, max_extension_len)
///     },
/// };
///
/// # let expected = vec![b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C'];
/// # assert_eq!(extended, expected);
/// ```
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
        (kmer_idx, kmer) = nearest_unique_context(noisy_ms, sbwt, search_end..kmer_idx);
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

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {

    #[test]
    fn nearest_unique_context() {
        use crate::build;
        use crate::index::BuildOpts;
        use crate::index::query_sbwt;
        use crate::gap_filling::nearest_unique_context;
        use sbwt::SbwtIndexVariant;

        // Parameters       : k = 9, threshold = 4
        //
        // Ref sequence     : T,T,G,A,T,T,A,A,C,A,G,G,C,A,G,C,T,C,A,G,A,G,C,T,G
        // Query sequence   : T,T,G,A,T,G,T,A,C,A,G,A,C,A,G,C,T,G,A,G,A,G,C,T,G
        //
        // Search range     :                       - - - - - -
        // Output           :                 C,A,G,A,C,A,G,C,T

        let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'T',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'A',b'G',b'C',b'T',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'G',b'T',b'A',b'C',b'A',b'G',b'A',b'C',b'A',b'G',b'C',b'T',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ k: 9, build_select: true, ..Default::default() });

        let ms = query_sbwt(&reference, &sbwt, &lcs);

        let context = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                nearest_unique_context(&ms, &sbwt, 11..16)
            },
        };
        let expected = (16, vec![b'C',b'A',b'G',b'A',b'C',b'A',b'G',b'C',b'T']);
        assert_eq!(context, expected);
    }

    #[test]
    fn left_extend_kmer() {
        use crate::build;
        use crate::index::BuildOpts;
        use crate::gap_filling::left_extend_kmer;
        use sbwt::SbwtIndexVariant;

        // Parameters       : k = 6
        //
        // Sequence         : T,T,G,A,T,G,T,A,C,A,G,A,C,T,G,C,G,G,A,G,A,G,C,T,G
        //
        // Starting k-mer   :                     G,A,C,T,G,C
        // Search range     :     - - - - - - - - - - - - - -
        // Output           :     G,A,T,G,T,A,C,A,G,A,C,T,G,C

        let sequence: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'G',b'T',b'A',b'C',b'A',b'G',b'A',b'C',b'T',b'G',b'C',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let mut opts = BuildOpts::default();
        opts.k = 6;
        opts.build_select = true;
        opts.add_revcomp = false;
        let (sbwt, _) = build(&[sequence], opts);

        let query = vec![b'G',b'A',b'C',b'T',b'G',b'C'];
        let max_extension_len = 8;

        let extended = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                let kmer_interval = sbwt.search(&query);
                let kmer = sbwt.access_kmer(kmer_interval.unwrap().start);
                left_extend_kmer(&kmer, &sbwt, max_extension_len)
            },
        };
        let expected = vec![b'G',b'A',b'T',b'G',b'T',b'A',b'C',b'A',b'G',b'A',b'C',b'T',b'G',b'C'];
        assert_eq!(extended, expected);
    }

}
