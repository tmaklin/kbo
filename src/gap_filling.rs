// kbo: Spectral Burrows-Wheeler transform accelerated local alignment search
//
// Copyright 2025 Tommi Mäklin [tommi@maklin.fi].

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

use sbwt::SbwtIndexVariant;

/// Count overlaps between a sequence and the last elements of a _k_-mer.
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

/// Count overlaps between a sequence and the first elements of a _k_-mer.
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
/// use kbo::BuildOpts;
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

/// Left extends a _k_-mer until the SBWT interval becomes non-unique.
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
/// use kbo::BuildOpts;
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
///
/// Searches `sbwt` for a _k_-mer that overlaps the region in `ref_seq` at range
/// `gap_index` using the SBWT intervals in `ms`. The search starts from the
/// SBWT interval at `gap_index.end + search_radius` and terminates when a
/// _k_-mer that overlaps the contains `left_overlap_req` characters matching
/// `ref_seq` to left of `gap_index.start` and `right_overlap_req` matching
/// characters from `gap_index.end`.
///
/// If no such _k_-mer is found, the search terminates at `gap_index.end +
/// right_overlap_req`. The characters between `gap_index.start` and
/// `gap_index.end` are not required to match.
///
/// Returns a concatenated sequence consisting of the the `left_overlap_req`
/// matching characters, the characters corresponding to range `gap_index`, and
/// the `right_overlap_req` matching characters.
///
/// The size of the return value is always `left_overlap_req + gap_index.end -
/// gap_index.start + right_overlap_req`.
///
/// If no _k_-mer meeting the conditions is found, returns an empty sequence.
///
/// The returned _k_-mer is stored as a ASCII-encoded vector.
///
/// # Examples
/// ```rust
/// use kbo::BuildOpts;
/// use kbo::build;
/// use sbwt::SbwtIndexVariant;
/// use kbo::index::query_sbwt;
/// use kbo::gap_filling::left_extend_over_gap;
///
/// // Parameters       : k = 9, left_overlap_req = 4, right_overlap_req = 4
/// //
/// // Ref sequence     : T,T,G,A,T,T,A,A,C,A,G,G,C,T,G,C,G,C,A,G,A,G,C,T,G
/// // Query sequence   : T,T,G,A,T,G,T,A,C,A G,A,C,T,G,C,G,G,A,G,A,G,C,T,G
/// // Translation      : M,M,M,M,M,-,-,-,-,-,-,-,M,M,M,M,M,X,M,M,M,M,M,M,M
/// // Search range     :                               - - -
/// // Overlap sequence :   T,G,A,T,G,T,A,C,A,G,A,C,T,G,C
///
/// let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'T',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
/// let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'G',b'T',b'A',b'C',b'A',b'G',b'A',b'C',b'T',b'G',b'C',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
///
/// let mut opts = BuildOpts::default();
/// opts.k = 9;
/// opts.build_select = true;
/// opts.add_revcomp = false;
/// let (sbwt, lcs) = build(&[query], opts);
///
/// let ms = query_sbwt(&reference, &sbwt, &lcs);
///
/// let overlap_seq = match sbwt {
///     SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
///         left_extend_over_gap(&ms, &reference, &sbwt, 4, 4, 5..12, 6)
///     },
/// };
///
/// # let expected = vec![b'T',b'G',b'A',b'T',b'G',b'T',b'A',b'C',b'A',b'G',b'A',b'C',b'T',b'G',b'C'];
/// # assert_eq!(overlap_seq, expected);
/// ```
///
pub fn left_extend_over_gap(
    ms: &[(usize, Range<usize>)],
    ref_seq: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    left_overlap_req: usize,
    right_overlap_req: usize,
    gap_index: Range<usize>,
    search_radius: usize,
) -> Vec<u8> {
    let k = sbwt.k();
    assert!(k > 0);
    assert!(ms.len() == ref_seq.len());
    assert!(left_overlap_req <= gap_index.start);
    assert!(right_overlap_req <= ref_seq.len() - gap_index.end);
    assert!(gap_index.end > gap_index.start);
    assert!(gap_index.end < ms.len());

    let search_start = (gap_index.end + search_radius).min(ref_seq.len() - 1);
    let search_end = gap_index.end + right_overlap_req;

    let mut kmer: Vec<u8> = Vec::with_capacity(k);
    let mut kmer_idx = search_start;
    while kmer_idx >= search_end {
        (kmer_idx, kmer) = nearest_unique_context(ms, sbwt, search_end..kmer_idx);
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

/// Refines a translated alignment by filling in gaps.
///
/// Attempts to resolve gaps '-'s in `translation` by using the
/// colexicographic intervals stored in the second component of
/// `noisy_ms` tuples to extract the _k_-mers that overlap the gaps
/// from the SBWT index `query_sbwt`.
///
/// The SBWT index must have [select
/// support](https://docs.rs/sbwt/latest/sbwt/struct.SbwtIndexBuilder.html)
/// enabled.
///
/// If the gap can be overlapped with a single _k_-mer so that the _k_-mer
/// contains `threshold` bases matching the bases in `ref_seq` to both left and
/// right of the gap, the single _k_-mer will be considered valid and its
/// contents used to fill in the gap.
///
/// If the gap cannot be overlapped with a single _k_-mer, the function will
/// attempt to find a _k_-mer that contains `threshold` overlaps to right of the
/// gap in `ref_seq` and extend this _k_-mer leftwards until it also overlaps
/// `ref_seq` to left of the gap
///
/// In case extending the _k_-mer was required, the _k_-mer is used to fill the
/// gap if the _k_-mer:
///
/// - contains no indels in the gap bases, **and**
///
/// one of the following is true:
///
/// - mismatches `ref_seq` only in the first and last gap base, **or**
/// - contains significant matching sections (determined using `max_err_prob`).
///
/// Use [variant calling](crate::variant_calling) and [add_variants](crate::translate::add_variants) in addition to fill_gaps if
/// resolving indels is required.
///
/// Returns a refined translation where the gaps characters '-' have been
/// replaced whenever possible.
///
/// # Examples
/// ```rust
/// use kbo::build;
/// use kbo::BuildOpts;
/// use kbo::index::query_sbwt;
/// use kbo::derandomize::derandomize_ms_vec;
/// use kbo::derandomize::random_match_threshold;
/// use kbo::translate::translate_ms_vec;
/// use kbo::gap_filling::fill_gaps;
/// use sbwt::SbwtIndexVariant;
///
/// // Parameters       : k = 9, threshold = 4
/// //
/// // Ref sequence     : T,T,G,A,T,T,A,A,C,A,G,G,C,T,G,C,G,C,A,G,A,G,C,T,G
/// // Query sequence   : T,T,G,A,T,G,T,A,C,A,G,A,C,T,G,C,G,G,A,G,A,G,C,T,G
/// //
/// // Translation      : M,M,M,M,M,-,-,-,-,-,-,-,M,M,M,M,M,X,M,M,M,M,M,M,M
/// // With gap filling : M,M,M,M,M,G,T,M,M,M,M,A,M,M,M,M,M,G,M,M,M,M,M,M,M
///
/// let query =     b"TTGATGTACAGACTGCGGAGAGCTG".to_vec();
/// let reference = b"TTGATTAACAGGCTGCGCAGAGCTG".to_vec();
///
/// let threshold = 4;
/// let mut opts = BuildOpts::default();
/// opts.k = 9;
/// opts.build_select = true;
/// let (sbwt, lcs) = build(&[query], opts);
///
/// let k = match sbwt {
///     SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
///         sbwt.k()
///     },
/// };
///
/// let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
/// let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);
/// let translated = translate_ms_vec(&derand_ms, k, threshold);
///
/// let gaps_filled = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.001_f64);
///
/// # let expected = b"MMMMMGTMMMMAMMMMMGMMMMMMM".to_vec();
/// # assert_eq!(gaps_filled, expected.iter().map(|x| *x as char).collect::<Vec<char>>());
/// ```
///
pub fn fill_gaps(
    translation: &[char],
    noisy_ms: &[(usize, Range<usize>)],
    ref_seq: &[u8],
    query_sbwt: &SbwtIndexVariant,
    threshold: usize,
    max_err_prob: f64,
) -> Vec<char> {
    let n_elements = translation.len();
    assert!(!translation.is_empty());
    assert!(translation.len() == noisy_ms.len());
    let k = match query_sbwt {
        SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
            assert!(sbwt.k() > 0);
            sbwt.k()
        },
    };

    let mut refined = translation.to_vec().clone();
    match query_sbwt {
        SbwtIndexVariant::SubsetMatrix(ref sbwt) => {

            let mut i: usize = threshold + 1;
            while i < refined.len() - threshold {
                if refined[i - 1] == '-' || refined[i - 1] == 'X' {
                    // Figure out how long the gap is
                    let start_index = i - 1;
                    while i < n_elements && refined[i] == '-' {
                        i += 1;
                    }
                    let end_index = i.min(refined.len() - threshold);

                    let overlap_without_extend = end_index - start_index + 2*threshold <= k;
                    let search_radius = k - (threshold * overlap_without_extend as usize);
                    let kmer = left_extend_over_gap(noisy_ms, ref_seq, sbwt, threshold, threshold, start_index..end_index, search_radius);

                    // Check if we want to use this k-mer
                    let kmer_found = !kmer.is_empty() && !kmer.contains(&b'$');
                    let no_indels = kmer.len() == threshold + (end_index - start_index) + threshold;

                    let matching_bases: Vec<bool> = kmer[threshold.min(kmer.len())..(threshold + end_index - start_index).min(kmer.len())]
                        .iter()
                        .zip(ref_seq[start_index..end_index].iter())
                        .map(|(kmer_nt, ref_nt)| kmer_nt == ref_nt)
                        .collect();

                    let total_overlaps: usize = matching_bases.iter().filter(|x| **x).count();
                    let log_probs: f64 = matching_bases.windows(2).scan(0, |consecutive_overlaps, x| {
                        if x[0] && x[1] {
                            *consecutive_overlaps += 1;
                            Some(0.0)
                        } else {
                            let log_prob = if *consecutive_overlaps > 0 {
                                crate::derandomize::log_rm_max_cdf(&*consecutive_overlaps + 1, 4, 1)
                            } else {
                                0.0
                            };
                            *consecutive_overlaps = 0;
                            Some(log_prob)
                        }
                    }).sum();

                    let fill_overlaps = log_probs > (-max_err_prob).ln_1p();
                    let fill_flanked = !matching_bases.is_empty() && !matching_bases[0] && !matching_bases[matching_bases.len() - 1] && total_overlaps + 2 == end_index - start_index;

                    let pass_checks = kmer_found && no_indels && (overlap_without_extend || fill_overlaps || fill_flanked);

                    if pass_checks {
                        refined[start_index..end_index]
                            .iter_mut()
                            .zip(kmer[threshold..(threshold + end_index - start_index)].iter())
                            .zip(ref_seq[start_index..end_index].iter())
                            .for_each(|((refined_val, kmer_nt), ref_nt)| {
                                *refined_val = if kmer_nt == ref_nt { 'M' } else { *kmer_nt as char }
                            });
                    }
                }
                i += 1;
            }
        },
    };
    refined
}

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {

    #[test]
    fn nearest_unique_context() {
        use crate::build;
        use crate::BuildOpts;
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
        use crate::BuildOpts;
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

    #[test]
    fn left_extend_over_gap() {
        use crate::BuildOpts;
        use crate::build;
        use sbwt::SbwtIndexVariant;
        use crate::index::query_sbwt;
        use crate::gap_filling::left_extend_over_gap;

        // Parameters       : k = 5, left_overlap_req = 3, right_overlap_req = 3
        //
        // Ref sequence     : T,T,G,A,A,C,A,G,G,C,T,G,C,G,C,A,G,A,G,C,T,G
        // Query sequence   : T,T,G,A,T,C,T,G,G,C,T,G,C,G,G,A,G,A,G,C,T,G
        // Translation      : M,M,M,M,-,-,-,M,M,M,M,M,M,M,X,M,M,M,M,M,M,M
        // Search range     :                   - - - - -
        // Overlap sequence :   T,G,A,T,C,T,G,G,C

        let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'C',b'T',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];

        let mut opts = BuildOpts::default();
        opts.k = 5;
        opts.build_select = true;
        opts.add_revcomp = false;
        let (sbwt, lcs) = build(&[query], opts);

        let ms = query_sbwt(&reference, &sbwt, &lcs);

        let overlap_seq = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                left_extend_over_gap(&ms, &reference, &sbwt, 3, 3, 4..7, 4)
            },
        };

        let expected = vec![b'T',b'G',b'A',b'T',b'C',b'T',b'G',b'G',b'C'];
        assert_eq!(overlap_seq, expected);

    }

    #[test]
    fn fill_gaps() {
        use crate::build;
        use crate::BuildOpts;
        use crate::index::query_sbwt;
        use crate::derandomize::derandomize_ms_vec;
        use crate::translate::translate_ms_vec;
        use super::fill_gaps;
        use sbwt::SbwtIndexVariant;

        // Parameters       : k = 7, threshold = 3
        //
        // Ref sequence     : T,T,G,A, T,T,G,G,C,T,G,G,G,C,A,G,A,G,C,T,G
        // Query sequence   : T,T,G,A,     G,G,C,T,G,G,G,G,A,G,A,G,C,T,G
        //
        // Result MS vector : 1,2,3,4, 1,2,3,3,3,4,4,4,4,3,1,2,3,4,4,4,4
        // Derandomized MS  : 1,2,3,4,-1,0,1,2,3,4,4,4,4,0,1,2,3,4,4,4,4
        // Translation      : M,M,M,M, -,-,M,M,M,M,M,M,M,X,M,M,M,M,M,M,M
        // Refined          : M,M,M,M, -,-,M,M,M,M,M,M,M,G,M,M,M,M,M,M,M
        // Changed this pos :                            |

        let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'G',b'G',b'C',b'T',b'G',b'G',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'T',b'G',b'G',b'C',b'T',b'G',b'G',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ k: 7, build_select: true, ..Default::default() });

        let k = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                sbwt.k()
            },
        };
        let threshold = 3;

        let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);
        let translated = translate_ms_vec(&derand_ms, k, threshold);

        let refined = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.001_f64);

        let expected = vec!['M','M','M','M','-','-','M','M','M','M','M','M','M','G','M','M','M','M','M','M','M'];
        assert_eq!(refined, expected);
    }

    // Test for case with matches-substitution-match-substitution-matches
    #[test]
    fn fill_gaps_with_clustered_changes() {
        use crate::build;
        use crate::BuildOpts;
        use crate::index::query_sbwt;
        use crate::derandomize::derandomize_ms_vec;
        use crate::translate::translate_ms_vec;
        use super::fill_gaps;
        use sbwt::SbwtIndexVariant;

        // Parameters       : k = 9, threshold = 3
        //
        // Ref sequence     : T,T,G,A, A, C,A,G,G,C,T,G,C,G,C,A,G,A,G,C,T,G
        // Query sequence   : T,T,G,A, T, C,T G,G,C,T,G,C,G,G,A,G,A,G,C,T,G
        //
        // Result MS vector : 1,2,3,4, 1, 1,1,2,2,3,4,4,4,4,3,1,2,3,4,4,4,4
        // Derandomized MS  : 1,2,3,4,-2,-1,0,1,2,3,4,4,4,4,0,1,2,3,4,4,4,4
        // Translation      : M,M,M,M, -, -,-,M,M,M,M,M,M,M,X,M,M,M,M,M,M,M
        // Refined          : M,M,M,M, T, C,T,M,M,M,M,M,M,M,G,M,M,M,M,M,M,M

        let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'C',b'T',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ k: 9, build_select: true, ..Default::default() });

        let k = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                sbwt.k()
            },
        };
        let threshold = 3;

        let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

        let translated = translate_ms_vec(&derand_ms, k, threshold);

        let refined = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.001_f64);

        let expected = vec!['M','M','M','M','T','M','T','M','M','M','M','M','M','M','G','M','M','M','M','M','M','M'];
        assert_eq!(refined, expected);
    }

    // Test for case with failure to resolve the gap due to ambiguous results
    #[test]
    fn fill_gaps_with_clustered_changes2() {
        use crate::build;
        use crate::BuildOpts;
        use crate::index::query_sbwt;
        use crate::derandomize::derandomize_ms_vec;
        use crate::translate::translate_ms_vec;
        use super::fill_gaps;
        use sbwt::SbwtIndexVariant;

        // Parameters       : k = 9, threshold = 3
        //
        // Ref sequence     : T,T,G,G, A, C,A,G,G,C,T,G,G,G,C,A,G,A,G,C,T,G
        // Query sequence   : T,T,G,G, G, C,T G,G,C,T,G,G,G,G,A,G,A,G,C,T,G
        //
        // Result MS vector : 1,2,3,4, 1, 1,1,2,2,3,4,4,4,4,3,1,2,3,4,4,4,4
        // Derandomized MS  : 1,2,3,4,-2,-1,0,1,2,3,4,4,4,4,0,1,2,3,4,4,4,4
        // Translation      : M,M,M,M, -, -,-,M,M,M,M,M,M,M,X,M,M,M,M,M,M,M
        // Refined          : M,M,M,M, G, M,T,M,M,M,M,M,M,M,R,R,M,M,M,M,M,M

        let query: Vec<u8> = vec![b'T',b'T',b'G',b'G',b'G',b'C',b'T',b'G',b'G',b'C',b'T',b'G',b'G',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let reference: Vec<u8> = vec![b'T',b'T',b'G',b'G',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'G',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ k: 9, build_select: true, ..Default::default() });

        let k = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                sbwt.k()
            },
        };
        let threshold = 3;

        let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

        let translated = translate_ms_vec(&derand_ms, k, threshold);

        let refined = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.001_f64);

        let expected = vec!['M','M','M','M','G','M','T','M','M','M','M','M','M','M','R','R','M','M','M','M','M','M'];
        assert_eq!(refined, expected);
    }

    // Test for left-extending over a longer gap
    #[test]
    fn fill_gaps_left_extend_short() {
        use crate::build;
        use crate::BuildOpts;
        use crate::index::query_sbwt;
        use crate::derandomize::derandomize_ms_vec;
        use crate::translate::translate_ms_vec;
        use super::fill_gaps;
        use sbwt::SbwtIndexVariant;

        // Parameters       : k = 9, threshold = 3
        //
        // Ref sequence     : T,T,G,A, A, C,A,G,G,C,T,G,C,G,C,A,G,A,G,C,T,G
        // Query sequence   : T,T,G,A, T, C,A G,A,C,T,G,C,G,G,A,G,A,G,C,T,G
        //
        // Result MS vector : 1,2,3,4, 1, 1,1,2,2,3,4,4,4,4,3,1,2,3,4,4,4,4
        // Derandomized MS  : 1,2,3,4,-2,-1,0,1,2,3,4,4,4,4,0,1,2,3,4,4,4,4
        // Translation      : M,M,M,M, -, -,-,-,-,M,M,M,M,M,X,M,M,M,M,M,M,M
        // Refined          : M,M,M,M, T, C,A,G,A,M,M,M,M,M,G,M,M,M,M,M,M,M

        let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'C',b'A',b'G',b'A',b'C',b'T',b'G',b'C',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ k: 9, build_select: true, ..Default::default() });

        let k = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                sbwt.k()
            },
        };
        let threshold = 3;

        let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

        let translated = translate_ms_vec(&derand_ms, k, threshold);

        let refined = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.001_f64);

        let expected = vec!['M','M','M','M','T','M','M','M','A','M','M','M','M','M','G','M','M','M','M','M','M','M'];
        assert_eq!(refined, expected);
    }

    // Test for left-extending over a longer gap
    #[test]
    fn fill_gaps_left_extend_long() {
        use crate::build;
        use crate::BuildOpts;
        use crate::index::query_sbwt;
        use crate::derandomize::derandomize_ms_vec;
        use crate::translate::translate_ms_vec;
        use super::fill_gaps;
        use sbwt::SbwtIndexVariant;

        // Parameters       : k = 9, threshold = 4
        //
        // Ref sequence     : T,T,G,A,T,T,A,A,C,A,G,G,C,T,G,C,G,C,A,G,A,G,C,T,G
        // Query sequence   : T,T,G,A,T,G,T,A,C,A G,A,C,T,G,C,G,G,A,G,A,G,C,T,G
        //
        // Translation      : M,M,M,M,M,-,-,-,-,-,-,M,M,M,M,M,X,M,M,M,M,M,M,M
        // Refined          : M,M,M,M,M,T,A,C,A,G,A,M,M,M,M,M,G,M,M,M,M,M,M,M

        let query: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'G',b'T',b'A',b'C',b'A',b'G',b'A',b'C',b'T',b'G',b'C',b'G',b'G',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
        let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'T',b'A',b'A',b'C',b'A',b'G',b'G',b'C',b'T',b'G',b'C',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ k: 9, build_select: true, ..Default::default() });

        let k = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                sbwt.k()
            },
        };
        let threshold = 4;

        let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

        let translated = translate_ms_vec(&derand_ms, k, threshold);

        let refined = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.001_f64);

        let expected = vec!['M','M','M','M','M','G','T','M','M','M','M','A','M','M','M','M','M','G','M','M','M','M','M','M','M'];
        assert_eq!(refined, expected);
    }

    // More difficult test case matches-substitution-match-substitution-matches
    // Parameters       : k = 51, threshold = 23
    #[test]
    fn fill_gaps_with_clustered_changes_k51() {
        use crate::build;
        use crate::BuildOpts;
        use crate::index::query_sbwt;
        use crate::derandomize::derandomize_ms_vec;
        use crate::translate::translate_ms_vec;
        use super::fill_gaps;
        use sbwt::SbwtIndexVariant;

        let query: Vec<u8> = vec![b'T',b'C',b'A',b'A',b'G',b'A',b'T',b'G',b'C',b'T',b'T',b'G',b'G',b'T',b'A',b'T',b'G',b'G',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'A',b'G',b'A',b'C',b'A',b'T',b'C',b'A',b'T',b'C',b'G',b'G',b'A',b'T',b'A',b'T',b'T',b'A',b'C',b'A',b'T',b'G',b'T',b'T',b'A',b'A',b'G',b'T',b'G',b'T',b'A',b'T',b'T',b'A',b'A',b'G',b'T',b'C',b'T',b'T',b'G',b'A',b'A',b'G',b'A',b'T',b'G',b'A',b'A',b'T',b'T',b'T',b'A',b'A',b'A',b'C',b'T',b'G',b'G',b'A',b'A',b'G',b'A',b'A',b'A',b'T',b'T',b'C',b'A',b'A',b'G',b'A',b'G',b'A',b'A',b'T',b'A',b'A',b'T',b'G',b'A',b'T',b'A',b'G',b'T',b'T',b'T',b'C',b'T',b'T',b'A',b'T',b'T',b'A',b'G',b'A',b'T',b'T',b'T',b'A',b'A',b'A',b'T',b'G',b'A',b'A',b'G',b'A',b'A',b'G',b'A',b'A',b'G',b'G',b'T',b'C',b'T',b'A',b'A',b'T',b'C',b'G',b'C',b'A',b'C',b'G',b'T',b'G',b'T',b'T',b'A',b'A',b'C',b'T',b'T',b'T',b'A',b'G',b'T',b'A',b'C',b'G',b'A',b'T',b'T',b'G',b'T',b'G',b'C',b'A',b'G',b'G',b'A',b'A',b'A',b'C',b'A',b'G',b'G',b'A',b'T',b'T',b'T',b'G',b'T',b'A',b'A',b'C',b'T',b'G',b'G',b'T',b'T',b'A',b'T',b'A',b'T',b'C',b'G',b'C',b'T',b'G',b'T',b'G',b'T',b'T',b'A',b'C',b'A',b'T',b'G',b'A',b'C',b'G',b'T',b'A',b'A',b'C',b'T',b'G',b'A',b'A',b'C',b'A',b'A',b'C',b'A',b'A',b'C',b'A',b'A',b'G',b'T',b'T',b'G',b'A',b'A',b'C',b'G',b'T',b'G',b'A',b'G',b'C',b'G',b'T',b'C',b'G',b'T',b'G',b'A',b'A',b'T',b'T',b'T',b'G',b'T',b'T',b'G',b'C',b'C',b'A',b'A',b'T',b'G',b'T',b'A',b'T',b'C',b'A',b'C',b'A',b'T',b'G',b'A',b'G',b'T',b'T',b'A',b'C',b'G',b'T',b'G',b'C',b'T',b'C',b'C',b'T',b'T',b'T',b'A',b'A',b'C',b'T',b'T',b'C',b'T',b'A',b'T',b'G',b'A',b'A',b'T',b'A',b'G',b'T',b'T',b'A',b'C',b'A',b'T',b'T',b'G',b'A',b'A',b'G',b'C',b'A',b'C',b'T',b'T',b'G',b'A',b'A',b'G',b'A',b'A',b'G',b'G',b'T',b'G',b'C',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'G',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'A',b'C',b'T',b'T',b'G',b'C',b'G',b'C',b'C',b'A',b'C',b'A',b'A',b'T',b'T',b'T',b'T',b'T',b'A',b'T',b'C',b'T',b'G',b'T',b'T',b'A',b'C',b'C',b'C',b'G',b'T',b'G',b'A',b'A',b'G',b'A',b'A',b'A',b'C',b'A',b'G',b'A',b'A',b'C',b'G',b'A',b'A',b'T',b'G',b'A',b'T',b'T',b'C',b'G',b'A',b'C',b'T',b'G',b'G',b'T',b'C',b'A',b'A',b'T',b'G',b'A',b'C',b'T',b'T',b'G',b'C',b'T',b'A',b'C',b'A',b'G',b'T',b'T',b'A',b'T',b'C',b'T',b'A',b'A',b'A',b'A',b'T',b'G',b'G',b'A',b'T',b'A',b'A',b'T',b'G',b'A',b'G',b'T',b'C',b'T',b'G',b'A',b'T',b'C',b'A',b'A',b'A',b'T',b'C',b'A',b'A',b'C',b'A',b'A',b'A',b'G',b'A',b'A',b'A',b'T',b'T',b'A',b'T',b'C',b'G',b'A',b'C',b'T',b'T',b'T',b'A',b'A',b'C',b'A',b'T',b'G',b'T',b'T',b'C',b'A',b'T',b'T',b'A',b'A',b'T',b'A',b'A',b'A',b'A',b'T',b'T',b'A',b'T',b'T',b'A',b'A',b'T',b'C',b'G',b'A',b'C',b'A',b'T',b'G',b'A',b'A',b'A',b'T',b'G',b'T',b'C',b'T',b'G',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'T',b'A',b'C',b'A',b'A',b'C',b'A',b'T',b'T',b'T',b'A',b'T',b'T',b'C',b'G',b'A',b'G',b'A',b'T',b'A',b'T',b'T',b'C',b'C',b'G',b'A',b'A',b'A',b'A',b'A',b'G',b'A',b'C',b'G',b'A',b'T',b'T',b'T',b'T',b'C',b'A',b'C',b'A',b'G',b'A',b'A',b'T',b'T',b'T',b'G',b'A',b'T',b'C',b'C',b'T',b'G',b'A',b'T',b'A',b'A',b'A',b'A',b'T',b'G',b'A',b'C',b'G',b'C',b'A',b'A',b'G',b'T',b'A',b'T',b'T',b'T',b'G',b'A',b'T',b'A',b'A',b'T',b'G',b'T',b'C',b'A',b'T',b'T',b'A',b'C',b'A',b'A',b'A',b'T',b'G',b'C',b'G',b'A',b'T',b'G',b'A',b'A',b'A',b'T',b'A',b'T',b'T',b'C',b'T',b'A',b'G',b'A',b'G',b'G',b'C',b'G',b'A',b'T',b'A',b'A',b'A',b'C',b'G',b'T',b'G',b'T',b'C',b'G',b'A',b'G',b'T',b'T',b'C',b'C',b'A',b'C',b'G',b'T',b'G',b'A',b'A',b'A',b'C',b'A',b'A',b'A',b'A',b'T',b'C',b'C',b'A',b'C',b'T',b'T',b'T',b'A',b'T',b'A',b'A',b'T',b'C',b'G',b'A',b'A',b'T',b'G',b'A',b'C',b'G',b'A',b'T',b'T',b'C',b'G',b'T',b'A',b'T',b'T',b'A',b'A',b'A',b'G',b'A',b'T',b'A',b'A',b'T',b'G',b'G',b'C',b'A',b'T',b'T',b'G',b'G',b'T',b'A',b'T',b'T',b'C',b'C',b'T',b'A',b'T',b'C',b'A',b'A',b'T',b'A',b'A',b'A',b'G',b'T',b'C',b'G',b'A',b'T',b'A',b'A',b'G',b'A',b'T',b'A',b'T',b'T',b'C',b'G',b'A',b'C',b'C',b'G',b'A',b'T',b'T',b'C',b'T',b'A',b'T',b'C',b'G',b'T',b'G',b'T',b'A',b'G',b'A',b'T',b'A',b'A',b'G',b'G',b'C',b'A',b'C',b'G',b'T',b'A',b'C',b'G',b'C',b'G',b'T',b'A',b'A',b'A',b'A',b'T',b'G',b'G',b'G',b'T',b'G',b'G',b'T',b'A',b'C',b'T',b'G',b'G',b'A',b'T',b'T',b'A',b'G',b'G',b'A',b'C',b'T',b'A',b'G',b'C',b'C',b'A',b'T',b'T',b'T',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'G',b'A',b'T',b'T',b'G',b'T',b'G',b'G',b'A',b'A',b'G',b'C',b'G',b'C',b'A',b'C',b'A',b'A',b'T',b'G',b'G',b'T',b'C',b'G',b'T',b'A',b'T',b'T',b'T',b'G',b'G',b'G',b'C',b'A',b'A',b'A',b'C',b'A',b'G',b'T',b'G',b'T',b'A',b'G',b'A',b'A',b'G',b'G',b'T',b'C',b'A',b'A',b'G',b'G',b'T',b'A',b'C',b'A',b'T',b'C',b'T',b'A',b'T',b'C',b'T',b'T',b'T',b'A',b'T',b'C',b'A',b'C',b'A',b'C',b'T',b'T',b'C',b'C',b'A',b'T',b'G',b'T',b'G',b'A',b'A',b'G',b'T',b'C',b'A',b'T',b'T',b'G',b'A',b'A',b'G',b'A',b'C',b'G',b'G',b'T',b'G',b'A',b'T',b'T',b'G',b'G',b'G',b'A',b'T',b'G',b'A',b'A',b'T',b'A',b'A',b'T',b'A',b'A',b'G',b'G',b'A',b'G',b'C',b'A',b'T',b'A',b'T',b'T',b'A',b'A',b'A',b'T',b'C',b'T',b'G',b'T',b'C',b'A',b'T',b'T',b'T',b'T',b'A',b'G',b'C',b'A',b'C',b'T',b'A',b'C',b'T',b'C',b'G',b'T',b'C',b'T',b'T',b'G',b'A',b'T',b'G',b'A',b'G',b'T',b'G',b'T'];
        let reference: Vec<u8> = vec![b'T',b'C',b'A',b'A',b'G',b'A',b'T',b'G',b'C',b'T',b'T',b'G',b'G',b'T',b'A',b'T',b'G',b'G',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'A',b'G',b'A',b'C',b'A',b'T',b'C',b'A',b'T',b'C',b'G',b'G',b'A',b'T',b'A',b'T',b'T',b'A',b'C',b'A',b'T',b'G',b'T',b'T',b'A',b'A',b'G',b'T',b'G',b'T',b'A',b'T',b'T',b'A',b'A',b'G',b'T',b'C',b'T',b'T',b'G',b'A',b'A',b'G',b'A',b'T',b'G',b'A',b'A',b'T',b'T',b'T',b'A',b'A',b'A',b'C',b'T',b'G',b'G',b'A',b'A',b'G',b'A',b'A',b'A',b'T',b'T',b'C',b'A',b'A',b'G',b'A',b'G',b'A',b'A',b'T',b'A',b'A',b'T',b'G',b'A',b'T',b'A',b'G',b'T',b'T',b'T',b'C',b'T',b'T',b'A',b'T',b'T',b'A',b'G',b'A',b'T',b'T',b'T',b'A',b'A',b'A',b'T',b'G',b'A',b'A',b'G',b'A',b'A',b'G',b'A',b'A',b'G',b'G',b'T',b'C',b'T',b'A',b'A',b'T',b'C',b'G',b'C',b'A',b'C',b'G',b'T',b'G',b'T',b'T',b'A',b'A',b'C',b'T',b'T',b'T',b'A',b'G',b'T',b'A',b'C',b'G',b'A',b'T',b'T',b'G',b'T',b'G',b'C',b'A',b'G',b'G',b'A',b'A',b'A',b'C',b'A',b'G',b'G',b'A',b'T',b'T',b'T',b'G',b'T',b'A',b'A',b'C',b'T',b'G',b'G',b'T',b'T',b'A',b'T',b'A',b'T',b'C',b'G',b'C',b'T',b'G',b'T',b'G',b'T',b'T',b'A',b'C',b'A',b'T',b'G',b'A',b'C',b'G',b'T',b'A',b'A',b'C',b'T',b'G',b'A',b'A',b'C',b'A',b'A',b'C',b'A',b'A',b'C',b'A',b'A',b'G',b'T',b'T',b'G',b'A',b'A',b'C',b'G',b'T',b'G',b'A',b'G',b'C',b'G',b'T',b'C',b'G',b'T',b'G',b'A',b'A',b'T',b'T',b'T',b'G',b'T',b'T',b'G',b'C',b'C',b'A',b'A',b'T',b'G',b'T',b'A',b'T',b'C',b'A',b'C',b'A',b'T',b'G',b'A',b'G',b'T',b'T',b'A',b'C',b'G',b'T',b'A',b'C',b'A',b'C',b'C',b'T',b'T',b'T',b'A',b'A',b'C',b'T',b'T',b'C',b'T',b'A',b'T',b'G',b'A',b'A',b'T',b'A',b'G',b'T',b'T',b'A',b'C',b'A',b'T',b'T',b'G',b'A',b'A',b'G',b'C',b'A',b'C',b'T',b'T',b'G',b'A',b'A',b'G',b'A',b'A',b'G',b'G',b'T',b'G',b'C',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'G',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'A',b'C',b'T',b'T',b'G',b'C',b'G',b'C',b'C',b'A',b'C',b'A',b'A',b'T',b'T',b'T',b'T',b'T',b'A',b'T',b'C',b'T',b'G',b'T',b'T',b'A',b'C',b'C',b'C',b'G',b'T',b'G',b'A',b'A',b'G',b'A',b'A',b'A',b'C',b'A',b'G',b'A',b'A',b'C',b'G',b'A',b'A',b'T',b'G',b'A',b'T',b'T',b'C',b'G',b'A',b'C',b'T',b'G',b'G',b'T',b'C',b'A',b'A',b'T',b'G',b'A',b'C',b'T',b'T',b'G',b'C',b'T',b'A',b'C',b'A',b'G',b'T',b'T',b'A',b'T',b'C',b'T',b'A',b'A',b'A',b'A',b'T',b'G',b'G',b'A',b'T',b'A',b'A',b'T',b'G',b'A',b'G',b'T',b'C',b'T',b'G',b'A',b'T',b'C',b'A',b'A',b'A',b'T',b'C',b'A',b'A',b'C',b'A',b'A',b'A',b'G',b'A',b'A',b'A',b'T',b'T',b'A',b'T',b'C',b'G',b'A',b'C',b'T',b'T',b'T',b'A',b'A',b'C',b'A',b'T',b'G',b'T',b'T',b'C',b'A',b'T',b'T',b'A',b'A',b'T',b'A',b'A',b'A',b'A',b'T',b'T',b'A',b'T',b'T',b'A',b'A',b'T',b'C',b'G',b'A',b'C',b'A',b'T',b'G',b'A',b'A',b'A',b'T',b'G',b'T',b'C',b'T',b'G',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'T',b'A',b'C',b'A',b'A',b'C',b'A',b'T',b'T',b'T',b'A',b'T',b'T',b'C',b'G',b'A',b'G',b'A',b'T',b'A',b'T',b'T',b'C',b'C',b'G',b'A',b'A',b'A',b'A',b'A',b'G',b'A',b'C',b'G',b'A',b'T',b'T',b'T',b'T',b'C',b'A',b'C',b'A',b'G',b'A',b'A',b'T',b'T',b'T',b'G',b'A',b'T',b'C',b'C',b'T',b'G',b'A',b'T',b'A',b'A',b'A',b'A',b'T',b'G',b'A',b'C',b'G',b'C',b'A',b'A',b'G',b'T',b'A',b'T',b'T',b'T',b'G',b'A',b'T',b'A',b'A',b'T',b'G',b'T',b'C',b'A',b'T',b'T',b'A',b'C',b'A',b'A',b'A',b'T',b'G',b'C',b'G',b'A',b'T',b'G',b'A',b'A',b'A',b'T',b'A',b'T',b'T',b'C',b'T',b'A',b'G',b'A',b'G',b'G',b'C',b'G',b'A',b'T',b'A',b'A',b'A',b'C',b'G',b'T',b'G',b'T',b'C',b'G',b'A',b'G',b'T',b'T',b'C',b'C',b'A',b'C',b'G',b'T',b'G',b'A',b'A',b'A',b'C',b'A',b'A',b'A',b'A',b'T',b'C',b'C',b'A',b'C',b'T',b'T',b'T',b'A',b'T',b'A',b'A',b'T',b'C',b'G',b'A',b'A',b'T',b'G',b'A',b'C',b'G',b'A',b'T',b'T',b'C',b'G',b'T',b'A',b'T',b'T',b'A',b'A',b'A',b'G',b'A',b'T',b'A',b'A',b'T',b'G',b'G',b'C',b'A',b'T',b'T',b'G',b'G',b'T',b'A',b'T',b'T',b'C',b'C',b'T',b'A',b'T',b'C',b'A',b'A',b'T',b'A',b'A',b'A',b'G',b'T',b'C',b'G',b'A',b'T',b'A',b'A',b'G',b'A',b'T',b'A',b'T',b'T',b'C',b'G',b'A',b'C',b'C',b'G',b'A',b'T',b'T',b'C',b'T',b'A',b'T',b'C',b'G',b'T',b'G',b'T',b'A',b'G',b'A',b'T',b'A',b'A',b'G',b'G',b'C',b'A',b'C',b'G',b'T',b'A',b'C',b'G',b'C',b'G',b'T',b'A',b'A',b'A',b'A',b'T',b'G',b'G',b'G',b'T',b'G',b'G',b'T',b'A',b'C',b'T',b'G',b'G',b'A',b'T',b'T',b'A',b'G',b'G',b'A',b'C',b'T',b'A',b'G',b'C',b'C',b'A',b'T',b'T',b'T',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'G',b'A',b'T',b'T',b'G',b'T',b'G',b'G',b'A',b'A',b'G',b'C',b'G',b'C',b'A',b'C',b'A',b'A',b'T',b'G',b'G',b'T',b'C',b'G',b'T',b'A',b'T',b'T',b'T',b'G',b'G',b'G',b'C',b'A',b'A',b'A',b'C',b'A',b'G',b'T',b'G',b'T',b'A',b'G',b'A',b'A',b'G',b'G',b'T',b'C',b'A',b'A',b'G',b'G',b'T',b'A',b'C',b'A',b'T',b'C',b'T',b'A',b'T',b'C',b'T',b'T',b'T',b'A',b'T',b'C',b'A',b'C',b'A',b'C',b'T',b'T',b'C',b'C',b'A',b'T',b'G',b'T',b'G',b'A',b'A',b'G',b'T',b'C',b'A',b'T',b'T',b'G',b'A',b'A',b'G',b'A',b'C',b'G',b'G',b'T',b'G',b'A',b'T',b'T',b'G',b'G',b'G',b'A',b'T',b'G',b'A',b'A',b'T',b'A',b'A',b'T',b'A',b'A',b'G',b'G',b'A',b'G',b'C',b'A',b'T',b'A',b'T',b'T',b'A',b'A',b'A',b'T',b'C',b'T',b'G',b'T',b'C',b'A',b'T',b'T',b'T',b'T',b'A',b'G',b'C',b'A',b'C',b'T',b'A',b'C',b'T',b'C',b'G',b'T',b'C',b'T',b'T',b'G',b'A',b'T',b'G',b'A',b'G',b'T',b'G',b'T'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ k: 51, build_select: true, ..Default::default() });

        let k = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                sbwt.k()
            },
        };
        let threshold = 23;

        let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

        let translated = translate_ms_vec(&derand_ms, k, threshold);

        let refined = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.0000001_f64);

        let expected = vec!['M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','G','M','T','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M'];
        assert_eq!(refined, expected);
    }

    // Test with longer query and reference & default SBWT parameters
    #[test]
    fn fill_gaps_default_build_opts() {
        use crate::build;
        use crate::BuildOpts;
        use crate::index::query_sbwt;
        use crate::derandomize::derandomize_ms_vec;
        use crate::translate::translate_ms_vec;
        use super::fill_gaps;
        use sbwt::SbwtIndexVariant;

        let query: Vec<u8> = vec![b'G',b'A',b'A',b'T',b'T',b'C',b'T',b'T',b'A',b'A',b'T',b'T',b'T',b'T',b'T',b'G',b'T',b'C',b'C',b'G',b'T',b'T',b'T',b'A',b'A',b'A',b'A',b'A',b'T',b'C',b'T',b'G',b'G',b'C',b'T',b'A',b'G',b'T',b'A',b'A',b'C',b'G',b'A',b'A',b'C',b'T',b'A',b'T',b'T',b'T',b'T',b'T',b'A',b'C',b'T',b'T',b'A',b'A',b'C',b'A',b'T',b'T',b'T',b'A',b'A',b'T',b'A',b'C',b'T',b'A',b'A',b'G',b'C',b'A',b'A',b'C',b'A',b'G',b'T',b'T',b'T',b'T',b'T',b'G',b'A',b'A',b'C',b'G',b'A',b'A',b'G',b'T',b'G',b'A',b'G',b'T',b'T',b'T',b'A',b'G',b'C',b'G',b'A',b'A',b'T',b'T',b'T',b'G',b'C',b'A',b'G',b'C',b'G',b'A',b'A',b'T',b'T',b'C',b'T',b'T',b'A',b'A',b'T',b'T',b'T',b'T',b'T',b'A',b'T',b'C',b'T',b'G',b'T',b'T',b'A',b'A',b'G',b'A',b'A',b'A',b'T',b'C',b'T',b'G',b'G',b'C',b'T',b'A',b'G',b'T',b'A',b'A',b'C',b'G',b'A',b'A',b'C',b'T',b'A',b'T'];
        let reference: Vec<u8> = vec![b'A',b'G',b'C',b'C',b'G',b'A',b'G',b'C',b'A',b'A',b'A',b'T',b'C',b'T',b'C',b'G',b'C',b'T',b'G',b'T',b'G',b'T',b'T',b'T',b'G',b'A',b'G',b'T',b'G',b'A',b'A',b'A',b'C',b'G',b'A',b'G',b'T',b'T',b'T',b'A',b'G',b'C',b'G',b'A',b'A',b'T',b'T',b'T',b'G',b'C',b'A',b'G',b'T',b'G',b'A',b'A',b'T',b'T',b'C',b'T',b'T',b'A',b'A',b'T',b'T',b'T',b'T',b'T',b'A',b'T',b'C',b'T',b'G',b'T',b'T',b'A',b'A',b'G',b'A',b'A',b'A',b'T',b'C',b'T',b'G',b'G',b'C',b'T',b'A',b'G',b'T',b'A',b'A',b'C',b'G',b'A',b'A',b'C',b'T',b'A',b'T',b'T',b'T',b'T',b'C',b'A',b'A',b'A',b'T',b'T',b'A',b'A',b'A',b'T',b'A',b'T',b'T',b'T',b'G',b'A',b'A',b'G',b'A',b'G',b'A',b'G',b'G',b'T',b'A',b'A',b'A',b'A',b'A',b'A',b'T',b'G',b'T',b'T',b'T',b'C',b'T',b'T',b'A',b'T',b'G',b'A',b'T',b'A',b'G',b'A',b'T',b'A',b'A',b'C',b'T',b'A',b'C',b'G',b'A',b'C'];

        let (sbwt, lcs) = build(&[query], BuildOpts{ build_select: true, ..Default::default() });

        let (k, threshold) = match sbwt {
            SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
                (sbwt.k(), crate::derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, 0.0000001_f64))
            },
        };

        let noisy_ms = query_sbwt(&reference, &sbwt, &lcs);
        let derand_ms = derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);
        let translated = translate_ms_vec(&derand_ms, k, threshold);

        let refined = fill_gaps(&translated, &noisy_ms, &reference, &sbwt, threshold, 0.0000001_f64);

        let expected = vec!['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'C', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'];
        assert_eq!(refined, expected);
    }
}
