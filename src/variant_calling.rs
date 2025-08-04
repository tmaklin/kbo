//! Call all variants between a query and a reference.

use std::{cmp::min, ops::Range};
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

/// Describes a variant between the query and the reference.
#[derive(Debug, Eq, PartialEq)]
pub struct Variant {
    /// A position in the query that does not match the reference.
    pub query_pos: usize,

    /// The sequence of characters that are in the query but not the reference.
    /// If empty, the variant is a deletion.
    pub query_chars: Vec<u8>,

    /// The sequence of characters that are in the reference but not the query.
    /// If empty, the variant is an insertion.
    pub ref_chars: Vec<u8>,
}

/// Returned from [resolve_variant] if the variant can't be resolved.
///
/// Codes:
/// - 0: Generic failure.
/// - 1: Weird case, query_overlap == ref_overlap.
///
/// This error can be ignored if it's okay to leave the variant unresolved.
///
#[derive(Debug,Clone)]
pub struct ResolveVariantErr {
    code: usize,
    message: String,
}

impl std::fmt::Display for ResolveVariantErr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let err_msg = match self.code {
            1 => format!("{}: ", self.message),
            _ => String::new(),
        };
        write!(f, "{}Could not resolve variant.", err_msg)
    }
}

// Pads with dollars from the left if there is not full k-mer
fn get_kmer_ending_at(query: &[u8], end_pos: usize, k: usize) -> Vec<u8> {
    let mut query_kmer = Vec::<u8>::new();
    if end_pos >= k-1 {
        query_kmer.extend(&query[end_pos + 1 - k .. end_pos+1]);
    } else {
        let n_dollars = -(end_pos as isize - k as isize + 1);
        assert!(n_dollars > 0);
        query_kmer.resize(n_dollars as usize, b'$');
        query_kmer.extend(&query[0..end_pos+1]);
    }
    assert!(query_kmer.len() == k);
    query_kmer
}

fn longest_common_suffix(x: &[u8], y: &[u8]) -> usize {
    let mut len = 0_usize;
    for i in 0..min(x.len(), y.len()) {
        if x[x.len() - 1 - i] == y[y.len() - 1 - i] {
            // The if-check will not go to negatives because of the constraint on i
            len += 1;
        } else {
            break;
        }
    }
    len
}

fn get_rightmost_significant_peak(ms: &[(usize, Range<usize>)], significant_match_threshold: usize) -> Option<usize> {
    assert!(!ms.is_empty());
    for i in (0..ms.len()-1).rev() {
        let here = ms[i].0;
        let next = ms[i+1].0;
        if here >= significant_match_threshold && here > next {
            return Some(i);
        }
    }
    None
}

/// Resolve variation between a query _k_-mer and a reference _k_-mer.
///
/// Resolve a variant between the query kmer and the reference k-mer that occurs just before
/// the longest common suffix of the two k-mers. This is a low-level subroutine and there are 
/// a lot of assumptions on the parameters, so only call this if you know what you are doing. 
/// To call and resolve all variants in a longer query, use [call_variants] instead.
/// 
/// Returns the variant sequences in the query and in the reference.
///
/// The function may fail, in which case a [ResolveVariantErr] is returned. The
/// error may be ignored to leave the variant unresolved.
///
/// # Parameters:
/// * query_kmer: The query k-mer
/// * ref_kmer: The reference k-mer
/// * ms_vs_query: An array of elements (d_i, l_i..r_i) such that d_i is the length of the
///                longest substring of query_kmer ending at position i that matches the
///                reference SBWT, and l_i..r_i is the colexicographic range of that match.
/// * ms_vs_ref: the same as ms_vs_query, but with the roles of the query and the reference swapped.
/// * significant_match_threshold: length of a match that is considered statistically significant.
///                                Use e.g. [random_match_threshold](crate::derandomize::random_match_threshold) to determine this.
///
/// # Examples
/// ```rust
/// use kbo::variant_calling::resolve_variant;
/// use kbo::variant_calling::Variant;
/// use sbwt::BitPackedKmerSortingMem;
/// use sbwt::SbwtIndexBuilder;
/// use sbwt::StreamingIndex;
///
/// //                                                                 v here
/// let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
/// let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCTTCAACGTG";
///
/// let k = 20;
/// let threshold = 5;
///
/// let (sbwt_ref, lcs_ref) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[reference]);
/// let (sbwt_query, lcs_query) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[query]);
///
/// let index_ref = StreamingIndex::new(&sbwt_ref, lcs_ref.as_ref().unwrap());
/// let index_query = StreamingIndex::new(&sbwt_query, lcs_query.as_ref().unwrap());
///
/// let ms_vs_ref = index_ref.matching_statistics(query);
/// let ms_vs_query = index_query.matching_statistics(reference);
///
/// let variant = resolve_variant(query, reference, &ms_vs_query, &ms_vs_ref, threshold).ok().unwrap();
///
/// // `variant` is a tuple with contents: `([84], [65])`
///
/// # assert_eq!(variant.0, [84]);
/// # assert_eq!(variant.1, [65]);
/// ```
///
pub fn resolve_variant(
    query_kmer: &[u8], 
    ref_kmer: &[u8], 
    ms_vs_query: &[(usize, Range<usize>)], // Slice of length k
    ms_vs_ref: &[(usize, Range<usize>)], // Slice of length k
    significant_match_threshold: usize,
) -> Result<(Vec<u8>, Vec<u8>), ResolveVariantErr> {

    let k = query_kmer.len();
    assert!(ref_kmer.len() == k);
    assert!(ms_vs_query.len() == k);
    assert!(ms_vs_ref.len() == k);

    let common_suffix_len = longest_common_suffix(query_kmer, ref_kmer);
    assert!(common_suffix_len > 0);

    let query_ms_peak = get_rightmost_significant_peak(ms_vs_ref, significant_match_threshold);
    let ref_ms_peak = get_rightmost_significant_peak(ms_vs_query, significant_match_threshold);

    if let (Some(query_ms_peak), Some(ref_ms_peak)) = (query_ms_peak, ref_ms_peak) {
        let suffix_match_start = k - common_suffix_len;

        // Negative gap means overlap 
        let query_gap = suffix_match_start as isize - query_ms_peak as isize - 1;
        let ref_gap = suffix_match_start as isize - ref_ms_peak as isize - 1;

        if query_gap > 0 && ref_gap > 0 {
            return Ok(
                (
                    query_kmer[query_ms_peak+1..suffix_match_start].to_vec(), // Query chars
                    ref_kmer[ref_ms_peak+1..suffix_match_start].to_vec() // Ref chars
                )
            )
        } else {
            let query_overlap = -query_gap;
            let ref_overlap = -ref_gap;
            if query_overlap == ref_overlap {
                return Err(ResolveVariantErr{ code: 1, message: "query_overlap == ref_overlap".to_string() })
            }
            let variant_len = (query_overlap - ref_overlap).unsigned_abs();
            if query_overlap > ref_overlap {
                // Deletion in query
                return Ok(
                    (
                        vec![], // Query chars
                        ref_kmer[ref_ms_peak + 1 .. ref_ms_peak + 1 + variant_len].to_vec() // Ref chars
                    )
                );
            }  else {
                // Insertion in query
                return Ok(
                    (
                        query_kmer[query_ms_peak + 1 .. query_ms_peak + 1 + variant_len].to_vec(), // Query chars
                        vec![] // Ref chars
                    )
                );
            }
        }

    }
    
    Err(ResolveVariantErr{ code: 0, message: "".to_string() }) // Could not resolve variant
}

/// Call variants between a query and a reference sequence.
///
/// Resolves single and multi base substitutions, insertions, and deletions by
/// querying the ref and query sequences against each other and comparing the
/// matching statistic vectors and calling [resolve_variant].
///
/// Returns a vector containing the [variants](Variant), if any.
///
/// # Parameters:
/// * sbwt_ref: the SBWT index of the reference
/// * sbwt_lcs: the LCS array of the reference
/// * sbwt_query: the SBWT index of the query
/// * lcs_query: the LCS array of the query
/// * query: the query sequence as ASCII characters
/// * max_error_prob: the p-value for statisticially significant matches, e.g. 1e-8
///
/// # Examples
/// ```rust
/// use kbo::variant_calling::call_variants;
/// use kbo::variant_calling::Variant;
/// use sbwt::{BitPackedKmerSortingMem, SbwtIndexBuilder};
///
/// //                                 deleted characters    substituted        inserted
/// //                                        v                 v                v
/// let reference = b"TCGTGGATCGATACACGCTAGCAGGCTGACTCGATGGGATACTATGTGTTATAGCAATTCGGATCGATCGA";
/// let query =     b"TCGTGGATCGATACACGCTAGCCTGACTCGATGGGATACCATGTGTTATAGCAATTCCGGATCGATCGA";
///
/// let max_error_prob = 0.001;
/// let k = 20;
///
/// let (sbwt_ref, lcs_ref) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[reference]);
/// let (sbwt_query, lcs_query) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[query]);
///
/// let variants = call_variants(&sbwt_ref, lcs_ref.as_ref().unwrap(), &sbwt_query, lcs_query.as_ref().unwrap(), query, max_error_prob);
///
/// // `variants` is a vector with elements:
/// // - Variant { query_pos: 22, query_chars: [], ref_chars: [65, 71, 71] }
/// // - Variant { query_pos: 39, query_chars: [67], ref_chars: [84] }
/// // - Variant { query_pos: 57, query_chars: [67], ref_chars: [] }
///
/// # assert_eq!(variants[0], Variant{query_pos: 22, query_chars: vec![], ref_chars: b"AGG".to_vec()});
/// # assert_eq!(variants[1], Variant{query_pos: 39, query_chars: vec![b'C'], ref_chars: vec![b'T']});
/// # assert_eq!(variants[2], Variant{query_pos: 57, query_chars: vec![b'C'], ref_chars: vec![]});
/// # assert_eq!(variants.len(), 3);
/// ```
///
pub fn call_variants(
    sbwt_ref: &SbwtIndex<SubsetMatrix>,
    lcs_ref: &LcsArray,
    sbwt_query: &SbwtIndex<SubsetMatrix>,
    lcs_query: &LcsArray,
    query: &[u8],
    max_error_prob: f64,
) -> Vec<Variant> {

    assert_eq!(sbwt_ref.k(), sbwt_query.k());
    let k = sbwt_ref.k();
    let d = crate::derandomize::random_match_threshold(k, sbwt_ref.n_kmers(), 4, max_error_prob);

    let mut calls: Vec<Variant> = vec![];

    let index_ref = StreamingIndex::new(sbwt_ref, lcs_ref);
    let index_query = StreamingIndex::new(sbwt_query, lcs_query);
    let ms_vs_ref = index_ref.matching_statistics(query);

    for i in 1..query.len() {
        if ms_vs_ref[i].0 < ms_vs_ref[i-1].0 && ms_vs_ref[i-1].0 >= d && ms_vs_ref[i].0 < d {
            // Go to closest unique match position to the right
            for j in i+1..min(i+k+1, query.len()) {
                if ms_vs_ref[j].0 >= d && ms_vs_ref[j].1.len() == 1 {
                    let ref_colex = ms_vs_ref[j].1.start;

                    let query_kmer = get_kmer_ending_at(query, j, k);
                    let ref_kmer = sbwt_ref.access_kmer(ref_colex);

                    // MS vectors for k-mers (Possible future work: could we reuse slices of the existing MS vectors?)
                    let ms_vs_ref = index_ref.matching_statistics(&query_kmer);
                    let ms_vs_query = index_query.matching_statistics(&ref_kmer);

                    if let Ok(var) = resolve_variant(&query_kmer, &ref_kmer, &ms_vs_query, &ms_vs_ref, d) {
                        calls.push(Variant{query_chars: var.0, ref_chars: var.1, query_pos: i});
                    }
                    break;
                }
            }


        }
    }

    calls
}

/// Edits a nucleotide sequence based on some variants.
///
/// Takes a nucleotide sequence `seq` and edits its contents using the variants
/// in `variants`. If `variants` lists all differences between `seq` and the
/// reference used to call them, the edited sequence will equal the reference.
///
/// Panics if `variants` is not sorted by query_pos.
///
/// Returns a vector containing the edited nucleotide sequence.
///
/// # Examples
///
/// ```rust
/// use kbo::variant_calling::edit;
/// use kbo::variant_calling::Variant;
///
/// let sequence = b"TCGTGGAGGTC".to_vec();
/// let edited = b"TAGGCGTGTATCAAAGC".to_vec();
///
/// let variants = vec![
///                     Variant{query_pos: 1, query_chars: vec![], ref_chars: b"AGG".to_vec()},
///                     Variant{query_pos: 5, query_chars: vec![b'G'], ref_chars: vec![b'T']},
///                     Variant{query_pos: 7, query_chars: b"GGTC".to_vec(), ref_chars: b"TCAAAGC".to_vec()},
///                     ];
///
/// let got = edit(&sequence, &variants);
///
/// # assert_eq!(got, edited);
/// ```
pub fn edit(
    seq: &[u8],
    variants: &[Variant],
) -> Vec<u8> {
    let mut edited: Vec<u8> = seq.to_vec();
    let mut prev_pos = seq.len();
    variants.iter().rev().for_each(|variant| {
        assert!(variant.query_pos < prev_pos);
        prev_pos = variant.query_pos;

        let q = variant.query_pos;
        let q_nts: Vec<u8> = variant.query_chars.to_vec();
        let r_nts: Vec<u8> = variant.ref_chars.to_vec();
        if r_nts.is_empty() {
            let mut prev = edited[0..q].to_vec();
            prev.extend(edited[(q + q_nts.len())..edited.len()].iter());
            edited = prev;
        } else if q_nts.is_empty() {
            let mut prev = edited[0..q].to_vec();
            prev.extend(r_nts.iter());
            prev.extend(edited[q..edited.len()].iter());
            edited = prev;
        } else if r_nts.len() == q_nts.len() {
            r_nts.iter().enumerate().for_each(|(idx, nt)| {
                edited[q + idx] = *nt;
            });
        } else if q_nts.len() != r_nts.len() {
            let mut prev = edited[0..q].to_vec();
            prev.extend(r_nts.iter());
            prev.extend(edited[(q + q_nts.len())..edited.len()].iter());
            edited = prev;
        }
    });
    edited
}

/// Reverts edits made to a nucleotide sequence based on some variants.
///
/// Takes an edited sequence `edited` and removes the edits given in `variants`.
/// If `variants` lists all differences between the variant calling reference
/// and a query, the returned sequence will equal the original query.
///
/// Panics if `variants` is not sorted by query_pos.
///
/// Returns a vector containing the unedited nucleotide sequence.
///
/// # Examples
///
/// ```rust
/// use kbo::variant_calling::revert_edits;
/// use kbo::variant_calling::Variant;
///
/// let edited =   b"TAGGCGTGTATCAAAGC".to_vec();
/// let original = b"TCGTGGAGGTC".to_vec();
///
/// let variants = vec![
///                     Variant{query_pos: 1, query_chars: vec![], ref_chars: b"AGG".to_vec()},
///                     Variant{query_pos: 5, query_chars: vec![b'G'], ref_chars: vec![b'T']},
///                     Variant{query_pos: 7, query_chars: b"GGTC".to_vec(), ref_chars: b"TCAAAGC".to_vec()},
///                     ];
///
/// let got = revert_edits(&edited, &variants);
///
/// # assert_eq!(got, original);
/// ```
pub fn revert_edits(
    edited: &[u8],
    variants: &[Variant],
) -> Vec<u8> {
    let mut original: Vec<u8> = edited.to_vec();
    let mut prev_pos = 0;
    variants.iter().for_each(|variant| {
        assert!(variant.query_pos >= prev_pos);
        prev_pos = variant.query_pos;

        let q = variant.query_pos;
        let q_nts: Vec<u8> = variant.query_chars.to_vec();
        let r_nts: Vec<u8> = variant.ref_chars.to_vec();
        if r_nts.is_empty() {
            let mut prev = original[0..q].to_vec();
            prev.extend(q_nts.iter());
            prev.extend(original[q..original.len()].iter());
            original = prev;
        } else if q_nts.is_empty() {
            let mut prev = original[0..q].to_vec();
            prev.extend(original[(q + r_nts.len())..original.len()].iter());
            original = prev;
        } else if r_nts.len() == q_nts.len() {
            q_nts.iter().enumerate().for_each(|(idx, nt)| {
                original[q + idx] = *nt;
            });
        } else if q_nts.len() != r_nts.len() {
            let mut prev = original[0..q].to_vec();
            prev.extend(q_nts.iter());
            prev.extend(original[(q + r_nts.len())..original.len()].iter());
            original = prev;
        }
    });
    original
}

#[cfg(test)]
mod tests {

    use random::Source;
    use sbwt::{BitPackedKmerSortingMem, SbwtIndexBuilder};

    use super::*;

    fn run_variant_calling(query: &[u8], reference: &[u8], k: usize, p_value: f64) -> Vec<Variant> {
        let (sbwt_ref, lcs_ref) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[reference]);
        let (sbwt_query, lcs_query) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new().k(k).build_lcs(true).build_select_support(true).run_from_slices(&[query]);

        call_variants(&sbwt_ref, lcs_ref.as_ref().unwrap(), &sbwt_query, lcs_query.as_ref().unwrap(), query, p_value)
    }

    #[test]
    fn test_single_base_substitution() {
        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCTTCAACGTG";
        //                                                                 ^ here

        let variants = run_variant_calling(query, reference, 20, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 49, ref_chars: vec![b'A'], query_chars: vec![b'T']}]);
    }

    #[test]
    fn test_multi_base_substitution() {
        //                                             **
        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAGCGTCTATTGTACCAATCGGCATCAACGTG";
        //                                             ***

        let variants = run_variant_calling(query, reference, 30, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 29, ref_chars: b"AA".to_vec(), query_chars: b"GCG".to_vec()}]);
    }

    #[test]
    fn test_multi_base_insertion_non_overlap_case() {
        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATATCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAGCGTCTATTGTACCAATCGGCATCAACGTG";
        //                                             ***

        let variants = run_variant_calling(query, reference, 30, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 29, ref_chars: b"".to_vec(), query_chars: b"GCG".to_vec()}]);
    }

    #[test]
    fn test_multi_base_insertion_overlap_case() {
        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAAAAAATCTATTGTACCAATCGGCATCAACGTG";
        //                                               ****

        let variants = run_variant_calling(query, reference, 30, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 31, ref_chars: b"".to_vec(), query_chars: b"AAAA".to_vec()}]);
    }

    #[test]
    fn test_single_base_insertion_non_overlap_case() {
        // Case: Non-overlapping reference intervals

        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCAGTCAACGTG";
        //                                                                  ^ here

        let variants = run_variant_calling(query, reference, 20, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 50, ref_chars: vec![], query_chars: vec![b'G']}]);
    }

    #[test]
    fn test_single_base_insertion_overlap_case() {
        // Case: Overlapping reference intervals

        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCAATCAACGTG";
        //                                                                  ^ here

        let variants = run_variant_calling(query, reference, 20, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 50, ref_chars: vec![], query_chars: vec![b'A']}]);
    }

    #[test]
    fn test_single_base_deletion_non_overlap_case() {
        // Case: Non-overlapping query intervals

        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCAGTCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
        //                                                                  ^ here

        let variants = run_variant_calling(query, reference, 20, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 50, ref_chars: vec![b'G'], query_chars: vec![]}]);
    }

    #[test]
    fn test_single_base_deletion_overlap_case() {
        // Case: overlapping query intervals

        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATTCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";
        //                                                                  ^ here

        let variants = run_variant_calling(query, reference, 20, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 51, ref_chars: vec![b'T'], query_chars: vec![]}]);
    }

    #[test]
    fn test_multi_base_deletion_non_overlap_case() {
        //                                             ***
        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAGCGTCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATATCTATTGTACCAATCGGCATCAACGTG";

        let variants = run_variant_calling(query, reference, 30, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 29, query_chars: b"".to_vec(), ref_chars: b"GCG".to_vec()}]);
    }

    #[test]
    fn test_multi_base_deletion_overlap_case() {
        //                                               ****
        let reference = b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAAAAAATCTATTGTACCAATCGGCATCAACGTG";
        let query =     b"GCGGGGCTGTTGACGTTTGGGGTTGAATAAATCTATTGTACCAATCGGCATCAACGTG";

        let variants = run_variant_calling(query, reference, 30, 0.001);
        dbg!(&variants);

        assert_eq!(variants, vec![Variant{query_pos: 31, query_chars: b"".to_vec(), ref_chars: b"AAAA".to_vec()}]);
    }

    #[test]
    fn test_variants_in_same_query() {
        //                                 deleted character    substituted        inserted
        //                                        v                 v                v
        let reference = b"TCGTGGATCGATACACGCTAGCAGGCTGACTCGATGGGATACTATGTGTTATAGCAATTCGGATCGATCGA";
        let query =      b"TCGTGGATCGATACACGCTAGCAGCTGACTCGATGGGATACCATGTGTTATAGCAATTCCGGATCGATCGA";

        let variants = run_variant_calling(query, reference, 20, 0.001);
        dbg!(&variants);

        assert_eq!(variants[0], Variant{query_pos: 24, query_chars: vec![], ref_chars: vec![b'G']});
        assert_eq!(variants[1], Variant{query_pos: 41, query_chars: vec![b'C'], ref_chars: vec![b'T']});
        assert_eq!(variants[2], Variant{query_pos: 59, query_chars: vec![b'C'], ref_chars: vec![]});
        assert_eq!(variants.len(), 3);
    }

    fn random_nucleotide(rng: &mut random::Default) -> u8{
        match rng.read_u64() % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!("Impossible math")
        }
    }

    #[test]
    fn test_long_generated_testcase(){
        let mut rng = random::Default::new([123412, 121232]);
        let mut reference = Vec::<u8>::new();
        let mut query = Vec::<u8>::new();

        let n = 100_000;
        let variant_spacing = 25;
        let k = 63;
        let p_value = 1e-8;

        let mut true_variants = Vec::<Variant>::new();

        for i in 0..n {
            if i > variant_spacing && i < n - variant_spacing && i % variant_spacing == 0 {
                let mut query_variant_len = rng.read_u64() % 4;
                let mut ref_variant_len = rng.read_u64() % 4;
                while query_variant_len == 0 && ref_variant_len == 0 {
                    // Reroll
                    query_variant_len = rng.read_u64() % 4;
                    ref_variant_len = rng.read_u64() % 4;
                }

                let mut query_variant: Vec<u8> = (0..query_variant_len).map(|_| random_nucleotide(&mut rng)).collect();
                let ref_variant: Vec<u8> = (0..ref_variant_len).map(|_| random_nucleotide(&mut rng)).collect();

                while !query_variant.is_empty() && !ref_variant.is_empty() 
                    && (query_variant.first() == ref_variant.first() || query_variant.last() == ref_variant.last()){
                    // Reroll the first and the last character so that we have a mismatch
                    *query_variant.last_mut().unwrap() = random_nucleotide(&mut rng);
                    *query_variant.first_mut().unwrap() = random_nucleotide(&mut rng);
                }

                true_variants.push(Variant{query_pos: query.len(), query_chars: query_variant.clone(), ref_chars: ref_variant.clone()});

                reference.extend(ref_variant.iter());
                query.extend(query_variant.iter());

                // Figure out if we had a pure insertion or deletion
                let mut insertion_seq: Option<Vec<u8>> = None;
                if query_variant.is_empty() && !ref_variant.is_empty() {
                    insertion_seq = Some(ref_variant.clone());
                } else if !query_variant.is_empty() && ref_variant.is_empty() {
                    insertion_seq = Some(query_variant.clone());
                }

                if let Some(insertion_seq) = insertion_seq {
                    // Continue with a character that mismatches both the start
                    // and the end of the variant to avoid accidental matches to
                    // the border of the variant in the following genome. This would
                    // make the variant description incorrect.
                    let mut c = random_nucleotide(&mut rng);
                    while c == *insertion_seq.first().unwrap() || c == *insertion_seq.last().unwrap() {
                        c = random_nucleotide(&mut rng);
                    }
                    query.push(c);
                    reference.push(c);
                }

                
            } else {
                let c = random_nucleotide(&mut rng);
                query.push(c);
                reference.push(c);
            }
        }

        //eprintln!("{}", String::from_utf8_lossy(&reference));
        //eprintln!("{}", String::from_utf8_lossy(&query));

        let calls = run_variant_calling(&query, &reference, k, p_value);

        let n_calls = calls.len();
        let mut n_correct = 0_usize;
        for variant_idx in 0..min(calls.len(), true_variants.len()) {
            let true_var = &true_variants[variant_idx];
            let our_var = &calls[variant_idx];
            eprintln!("true: ({}: {} -> {}), our: ({}: {} -> {})", 
                true_var.query_pos, String::from_utf8_lossy(&true_var.ref_chars), String::from_utf8_lossy(&true_var.query_chars),
                our_var.query_pos, String::from_utf8_lossy(&our_var.ref_chars), String::from_utf8_lossy(&our_var.query_chars));
            n_correct += (*true_var == *our_var) as usize;
        }

        eprintln!("{} true variants", true_variants.len());
        eprintln!("{} calls", calls.len());
        eprintln!("{} correct out of {}", n_correct, n_calls);
        assert_eq!(n_calls, n_correct);
    }

    #[test]
    fn edit_multi_base_substitution() {
        use super::edit;
        use super::Variant;

        let sequence = b"TCGTGGAGGTC".to_vec();
        let edited =   b"TCGTGTCGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 5, query_chars: b"GA".to_vec(), ref_chars: b"TC".to_vec()},
        ];

        let got = edit(&sequence, &variants);

        assert_eq!(got, edited);
    }

    #[test]
    fn edit_multi_base_deletion() {
        use super::edit;
        use super::Variant;

        let sequence = b"TCGGAGGTC".to_vec();
        let edited =   b"TCGGTTAGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 4, query_chars: vec![], ref_chars: b"TT".to_vec()},
        ];

        let got = edit(&sequence, &variants);

        assert_eq!(got, edited);
    }

    #[test]
    fn edit_multi_base_insertion() {
        use super::edit;
        use super::Variant;

        let sequence = b"TCGGGGAGGTC".to_vec();
        let edited =   b"TCGGAGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 4, query_chars: b"GG".to_vec(), ref_chars: vec![]},
        ];

        let got = edit(&sequence, &variants);

        assert_eq!(got, edited);
    }

    #[test]
    fn edit_flanking_variant() {
        use super::edit;
        use super::Variant;

        let sequence = b"TCGGCGAGGTC".to_vec();
        let edited =   b"TCGGAGTGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 4, query_chars: b"CGA".to_vec(), ref_chars: b"AGT".to_vec()},
        ];

        let got = edit(&sequence, &variants);

        assert_eq!(got, edited);
    }

    #[test]
    fn edit_substituting_insertion_into_query() {
        use super::edit;
        use super::Variant;

        let sequence = b"TCGGAGCGTTGGTGGTGGTC".to_vec();
        let edited =   b"TCGGAGTGGTGGTGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 5, query_chars: b"GCGTT".to_vec(), ref_chars: b"GT".to_vec()},
        ];

        let got = edit(&sequence, &variants);

        assert_eq!(got, edited);
    }

    #[test]
    fn edit_substituting_insertion_into_ref() {
        use super::edit;
        use super::Variant;

        let sequence = b"TCGGAGTGGTGGTGGTC".to_vec();
        let edited =   b"TCGGAGCGTTGGTGGTGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 5, ref_chars: b"GCGTT".to_vec(), query_chars: b"GT".to_vec()},
        ];

        let got = edit(&sequence, &variants);

        assert_eq!(got, edited);
    }

    #[test]
    fn revert_edits_multi_base_substitution() {
        use super::revert_edits;
        use super::Variant;

        let original = b"TCGTGGAGGTC".to_vec();
        let edited =   b"TCGTGTCGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 5, query_chars: b"GA".to_vec(), ref_chars: b"TC".to_vec()},
        ];

        let got = revert_edits(&edited, &variants);

        assert_eq!(got, original);
    }

    #[test]
    fn revert_edits_multi_base_deletion() {
        use super::revert_edits;
        use super::Variant;

        let original = b"TCGGAGGTC".to_vec();
        let edited =   b"TCGGTTAGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 4, query_chars: vec![], ref_chars: b"TT".to_vec()},
        ];

        let got = revert_edits(&edited, &variants);

        assert_eq!(got, original);
    }

    #[test]
    fn revert_edits_multi_base_insertion() {
        use super::revert_edits;
        use super::Variant;

        let original = b"TCGGGGAGGTC".to_vec();
        let edited =   b"TCGGAGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 4, query_chars: b"GG".to_vec(), ref_chars: vec![]},
        ];

        let got = revert_edits(&edited, &variants);

        assert_eq!(got, original);
    }

    #[test]
    fn revert_edits_flanking_variant() {
        use super::revert_edits;
        use super::Variant;

        let original = b"TCGGCGAGGTC".to_vec();
        let edited =   b"TCGGAGTGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 4, query_chars: b"CGA".to_vec(), ref_chars: b"AGT".to_vec()},
        ];

        let got = revert_edits(&edited, &variants);

        assert_eq!(got, original);
    }

    #[test]
    fn revert_edits_substituting_insertion_into_query() {
        use super::revert_edits;
        use super::Variant;

        let original = b"TCGGAGCGTTGGTGGTGGTC".to_vec();
        let edited =   b"TCGGAGTGGTGGTGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 5, query_chars: b"GCGTT".to_vec(), ref_chars: b"GT".to_vec()},
        ];

        let got = revert_edits(&edited, &variants);

        assert_eq!(got, original);
    }

    #[test]
    fn revert_edits_substituting_insertion_into_ref() {
        use super::revert_edits;
        use super::Variant;

        let original = b"TCGGAGTGGTGGTGGTC".to_vec();
        let edited =   b"TCGGAGCGTTGGTGGTGGTC".to_vec();

        let variants = vec![
            Variant{query_pos: 5, ref_chars: b"GCGTT".to_vec(), query_chars: b"GT".to_vec()},
        ];

        let got = revert_edits(&edited, &variants);

        assert_eq!(got, original);
    }
}
