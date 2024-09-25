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

pub mod derandomize;
pub mod format;
pub mod index;
pub mod translate;

/// Builds an SBWT index from some fasta or fastq files.
///
/// Reads all sequence data in `seq_files` and builds an SBWT index
/// with the parameters and resources specified in `build_opts` (see
/// [index::BuildOpts] for details).
///
/// All files and sequence data in `seq_files` are merged into the
/// same index. It is not possible extract the individual sequences
/// from the index after it has been built; use [matches] with the TODO
/// options instead if you need to know which reference sequences the
/// alignments are for.
///
/// TODO Describe map syntax in lib.rs documentation.
///
/// Returns a tuple containing the built
/// [sbwt::SbwtIndexVariant](https://docs.rs/sbwt/latest/sbwt/enum.SbwtIndexVariant.html)
/// and
/// [sbwt::LcsArray](https://docs.rs/sbwt/latest/sbwt/struct.LcsArray.html).
///
/// Panics if a file in `seq_files` is not readable or a valid FASTX
/// file.
///
/// # Examples
/// ```rust
/// use sablast::build;
/// use sablast::index::BuildOpts;
///
/// let inputs: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
///
/// let (sbwt_index, lcs_array) = build(&inputs, BuildOpts::default());
/// ```
///
pub fn build(
    seq_data: &[Vec<u8>],
    build_opts: index::BuildOpts,
) -> (sbwt::SbwtIndexVariant, sbwt::LcsArray) {
    index::build_sbwt_from_vecs(seq_data, &Some(build_opts))
}

/// Matches a query fasta or fastq file against an SBWT index.
///
/// Queries the sequence data in `query_seq` against the SBWT index
/// `sbwt` and its LCS array `lcs` using [index::query_sbwt]. Then,
/// derandomizes the resulting _k_-bounded matching statistics vector
/// using [derandomize::derandomize_ms_vec] and translates the
/// matching statistics to a character representation of the alignment
/// using [translate::translate_ms_vec].
///
/// Returns a vector containing the character representation of the
/// alignment.
///
/// Panics if the query file is not readable or if it's not a valid
/// FASTX file.
///
/// # Output format
/// See the documentation for [translate].
///
/// # Example
/// ```rust
/// use sablast::build;
/// use sablast::matches;
/// use sablast::index::BuildOpts;
///
/// let reference: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let (sbwt, lcs) = build(&reference, BuildOpts{ k: 3, ..Default::default() });
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let ms_vectors = matches(&query, &sbwt, &lcs);
/// // `ms_vectors` has ['-','-','-','-','-','-','-','-','-','M','M','M','-','-']
/// # assert_eq!(ms_vectors, vec!['-','-','-','-','-','-','-','-','-','M','M','M','-','-']);
/// ```
///
pub fn matches(
    query_seq: &[u8],
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<char> {
    let (k, threshold) = match sbwt {
	SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
	    (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, 0.0000001_f64))
	},
    };

    let noisy_ms = index::query_sbwt(query_seq, sbwt, lcs);
    let derand_ms = derandomize::derandomize_ms_vec(&noisy_ms, k, threshold);

    translate::translate_ms_vec(&derand_ms, k, threshold)
}

/// Maps a query sequence against a reference sequence.
///
/// Maps the sequence data in `ref_seq` against the SBWT index
/// `query_sbwt` and `query_lcs` and converts the alignment to a
/// mapping relative to `ref_seq`.
///
/// Return the reference sequence with characters that are not present
/// in the query masked with a '-'.
///
/// # Examples
/// ```rust
/// use sablast::build;
/// use sablast::map;
/// use sablast::index::BuildOpts;
///
/// let query: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let (sbwt_query, lcs_query) = build(&query, BuildOpts{ k: 3, ..Default::default() });
///
/// let reference = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let alignment = map(&reference, &sbwt_query, &lcs_query);
/// // `ms_vectors` has [45,45,45,45,45,45,45,45,45,65,71,71,45,45]
/// # assert_eq!(alignment, vec![45,45,45,45,45,45,45,45,45,65,71,71,45,45]);
/// ```
///
pub fn map(
    ref_seq: &[u8],
    query_sbwt: &sbwt::SbwtIndexVariant,
    query_lcs: &sbwt::LcsArray,
) -> Vec<u8> {
    let aln = matches(ref_seq, &query_sbwt, &query_lcs);
    format::relative_to_ref(ref_seq, &aln)
}

/// Finds the _k_-mers from an SBWT index in a query fasta or fastq file.
///
/// Aligns the sequence data in `query_seq` against the SBWT index
/// `sbwt` and its LCS array `lcs` using [matches]. Then uses
/// [format::run_lengths] to extract the local alignments from the
/// matching statistics.
///
/// Returns a vector of tuples, where each element represents a local
/// alignment block and contains the following values:
/// 1. Start of local alignment block in query (1-based indexing).
/// 2. End of local alignment block in query.
/// 3. Number of matches in the block.
/// 4. Number of mismatches and 1-character insertions in the block.
///
/// # Examples
/// ```rust
/// use sablast::build;
/// use sablast::find;
/// use sablast::index::BuildOpts;
///
/// let reference: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let (sbwt, lcs) = build(&reference, BuildOpts{ k: 3, ..Default::default() });
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let local_alignments = find(&query, &sbwt, &lcs);
/// // `local_alignments` has [(10, 12, 3, 0)]
/// # assert_eq!(local_alignments, vec![(10, 12, 3, 0)]);
/// ```
///
pub fn find(
    query_seq: &[u8],
    sbwt: &sbwt::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Vec<(usize, usize, usize, usize)> {
    let aln = matches(query_seq, sbwt, lcs);
    format::run_lengths(&aln)
}

// TODO Implement refining the translated MS vectors
// 1. Search for Xs
// 2. Extract 2*k+1 region centered on the X.
// 3. Map region to query.
// 4. Resolve SNP vs insertion.
// 5. If SNP get the base.
pub fn refine_translation(
    ref_seq: &[u8],
    query_seq: &[u8],
    translation: &[char],
    raw: &[u8]
) -> Vec<u8> {
    let k: usize = 31;

    let mut refined = raw.to_vec().clone();

    let mut i = refined.len() - 1;
    while i > 0 {
	if translation[i] == 'X' {
	eprintln!("{}/{}", i, refined.len());
	    let l_try: i64 = i as i64 - k as i64;
	    let l_start: usize = if l_try < 0 { 0 } else { l_try as usize };
	    let r_end: usize = if i + 1 + k > refined.len() { refined.len() - 1 } else { i + 1 + k };

	    let left = Vec::from_iter(ref_seq[l_start..i].iter().cloned());
	    let right = Vec::from_iter(ref_seq[(i + 1)..r_end].iter().cloned());

	    let (sbwt_left, lcs_left) = index::build_sbwt_from_vecs(&[left], &Some(index::BuildOpts::default()));
	    let where_left = find(query_seq, &sbwt_left, &lcs_left);

	    let (sbwt_right, lcs_right) = index::build_sbwt_from_vecs(&[right], &Some(index::BuildOpts::default()));
	    let where_right = find(query_seq, &sbwt_right, &lcs_right);

	    // This assumes matches are unique
	    if where_right[0].0 - 2 == where_left[0].1 {
		refined[i] = query_seq[where_right[0].0 - 1];
	    }
	} else {
	    refined[i] = raw[i];
	}
	i -= 1;
    }

    refined
}

pub fn map_refine(
    ref_seq: &[u8],
    query_seq: &[u8],
    query_sbwt: &sbwt::SbwtIndexVariant,
    query_lcs: &sbwt::LcsArray,
) -> Vec<u8> {
    let aln = matches(ref_seq, &query_sbwt, &query_lcs);

    let raw = format::relative_to_ref(ref_seq, &aln);
    let refined = refine_translation(ref_seq, query_seq, &aln, &raw);

    refined
}
