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
//! Converting alignment representations into various output formats.

/// Run length encoding for an alignment segment
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RLE {
    /// Start position (1-based indexing)
    pub start: usize,
    /// End position (1-based indexing)
    pub end: usize,
    /// Number of matching bases ('M' or 'R')
    pub matches: usize,
    /// Number of mismatching bases (not 'M', 'R', or '-')
    pub mismatches: usize,
    /// Number of _k_-mer jumps (count of double 'R's)
    pub jumps: usize,
    /// Total number of missing bases ('-')
    pub gap_bases: usize,
    /// Number of consecutive '-' runs in segment regardless of length
    pub gap_opens: usize,
}

impl Default for RLE {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::format::RLE::default();
    /// opts.start = 0;
    /// opts.end = 0;
    /// opts.matches = 0;
    /// opts.mismatches = 0;
    /// opts.jumps = 0;
    /// opts.gap_bases = 0;
    /// opts.gap_opens = 0;
    /// # let expected = kbo::format::RLE::default();
    /// # assert_eq!(opts, expected);
    /// ```
    ///
    fn default() -> RLE {
        RLE {
            start: 0,
            end: 0,
            matches: 0,
            mismatches: 0,
            jumps: 0,
            gap_bases: 0,
            gap_opens: 0,
        }
    }
}

/// Extracts run length encodings from a translated alignment.
///
/// Traverses the character representation of the alignment stored in `aln` and
/// counts the consecutive run lengths of sections that align to the reference.
/// A base is counted as aligned if its character representation is not '-' or ' '.
///
/// This function can be used for both plain and refined translations.
///
/// Returns a vector where each element contains 4 elements with the contents:
/// - Start position of a run of consecutive characters (indexing starts from 1).
/// - End position of a run of consecutive characters (indexing starts from 1).
/// - Number matching bases in the run (character representation is 'M' or 'R').
/// - Number mismatching bases in the run (character representation is not 'M', '-', or ' ').
///
/// The run length is the sum of the matching and mismatching bases.
///
/// # Examples
/// ## Extract run lengths from a character representation
/// ```rust
/// use kbo::format::run_lengths;
/// use kbo::format::RLE;
///
/// // Parameters       : k = 3, threshold = 2
/// //
/// // Ref sequence     : A,A,A,G,A,A,C,C,A,-,T,C,A, -,-,G,G,G, C,G
/// // Query sequence   : C,A,A,G,-,-,C,C,A,C,T,C,A, T,T,G,G,G, T,C
/// // Input MS         : 0,1,2,3,    1,2,3,0,1,2,3,-1,0,1,2,3,-1,0
/// // Translation      : X,M,M,R, , ,R,M,M,X,M,M,M,  , ,M,M,M,  ,
///
/// // Expected output: [(1,  11, 9, 2),
/// //                   (14, 16, 3, 0)]
///
/// let input: Vec<char> = vec!['X','M','M','R','R','M','M','X','M','M','M','-','-','M','M','M','-','-'];
/// let run_lengths = run_lengths(&input);
/// # let expected = vec![RLE{start: 1, end: 11, matches: 9, mismatches: 2, jumps: 1, gap_bases: 0, gap_opens: 0},
/// #                     RLE{start: 14, end: 16, matches: 3, mismatches : 0, jumps : 0, gap_bases : 0, gap_opens : 0}];
/// # assert_eq!(run_lengths, expected);
/// ```
///
pub fn run_lengths(
    aln: &[char],
) -> Vec<RLE> {
    let mut encodings: Vec<RLE> = Vec::new();

    let mut i = 0;
    let mut match_start: bool = false;
    while i < aln.len() {
        match_start = (aln[i] != '-' && aln[i] != ' ') && !match_start;
        let mut jumps = 0;
        if match_start {
            let start = i;
            let mut matches: usize = 0;
            while i < aln.len() && (aln[i] != '-' && aln[i] != ' ') {
                matches += (aln[i] == 'M' || aln[i] == 'R') as usize;
                jumps += (aln[i] == 'R') as usize;
                i += 1;
            }
            let rle: RLE = RLE{
                start: start + 1,
                end: i,
                matches,
                mismatches: i - start - matches,
                jumps: jumps / 2,
                gap_bases: 0,
                gap_opens: 0,
            };
            encodings.push(rle);
            match_start = false;
        } else {
            i += 1;
        }
    }
    encodings
}

/// Format a refined translation relative to the reference sequence.
///
/// Jointly reads nucleotides from the reference sequence `ref_seq` and the
/// character representation of an alignment `alignment` to determine the
/// nucleotide sequence of the alignment relative to the reference.
///
/// Returns a vector containing alignment of the query against `ref_seq` in a
/// gapped alignment format.
///
/// Valid characters in the return format are:
/// - 'A', 'C', 'G', 'T': the nucleotide in the query sequence.
/// - '-': gap in the query.
///
/// If `alignment` is a refined translation, the nucleotides ACGT in the return
/// value may differ from the reference sequence and the gaps '-' represent
/// sections absent from the query.
///
/// If `alignment` is an unrefined translation, the nucleotides ACGT are always
/// the same as in the reference sequence. Gaps '-' may represent either
/// unresolved SNPs, insertions, or absent sections.
///
/// # Examples
/// ## Format an unrefined translation
/// Unrefined translations do not resolve SNPs or insertions in query, instead
/// replacing them with gaps '-'.
/// ```rust
/// use kbo::translate::translate_ms_vec;
/// use kbo::format::relative_to_ref;
///
/// // Parameters       : k = 4, threshold = 3
/// //
/// // Ref sequence     : T,T,G,A, T,T,G,G,C,T,G,G,G,C,A,G,A,G,C,T,G
/// // Query sequence   : T,T,G,A,     G,G,C,T,G,G,G,G,A,G,A,G,C,T,G
/// //
/// // Result MS vector : 1,2,3,4, 1,2,3,3,3,4,4,4,4,3,1,2,3,4,4,4,4
/// // Derandomized MS  : 1,2,3,4,-1,0,1,2,3,4,4,4,4,0,1,2,3,4,4,4,4
/// // Translation      : M,M,M,M, -,-,M,M,M,M,M,M,M,X,M,M,M,M,M,M,M
/// // Formatted        : T,T,G,A, -,-,G,G,C,T,G,G,G,-,A,G,A,G,C,T,G
///
/// let reference: Vec<u8> = vec![b'T',b'T',b'G',b'A',b'T',b'T',b'G',b'G',b'C',b'T',b'G',b'G',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
/// let derand_ms: Vec<i64> = vec![1,2,3,4,-1,0,1,2,3,4,4,4,4,0,1,2,3,4,4,4,4];
/// let translated = translate_ms_vec(&derand_ms, 4, 3);
///
/// let relative = relative_to_ref(&reference, &translated);
/// // `relative` has [T,T,G,A,-,-,G,G,C,T,G,G,G,-,A,G,A,G,C,T,G]
/// # let expected = vec![b'T',b'T',b'G',b'A',b'-',b'-',b'G',b'G',b'C',b'T',b'G',b'G',b'G',b'-',b'A',b'G',b'A',b'G',b'C',b'T',b'G'];
/// # assert_eq!(relative, expected);
/// ```
///
/// ## Format a refined translation
/// ```rust
/// use kbo::format::relative_to_ref;
///
/// // Ref sequence     : A,A,A,G,A,A,C,C,A,  T,C,A,G,G,G,C,G
/// // Query sequence   : C,A,A,G,-,-,C,C,A,C,T,C,A,G,G,G,-,-
/// // Translation      : C,M,M,R,-,-,R,M,M, ,M,M,M,M,M,M,-,-
///
/// let reference: Vec<u8> = vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G'];
/// let refined: Vec<char> = vec!['C','M','M','R','-','-','R','M','M','M','M','M','M','M','M','-','-'];
///
/// let relative = relative_to_ref(&reference, &refined);
/// # let expected = vec![b'C',b'A',b'A',b'G',b'-',b'-',b'C',b'C',b'A',b'T',b'C',b'A',b'G',b'G',b'G',b'-',b'-'];
/// // relative.iter().zip(reference.iter()).zip(expected.iter()).for_each(|((x1,x2),x3)| eprintln!("{}\t{}\t{}", *x1 as char, *x2 as char, *x3 as char));
/// # assert_eq!(relative, expected);
/// ```
///
pub fn relative_to_ref(
    ref_seq: &[u8],
    alignment: &[char],
) -> Vec<u8> {
    ref_seq.iter().zip(alignment.iter()).map(|x| {
        if *x.1 == 'M' || *x.1 == 'R' {
            *x.0
        } else if *x.1 == 'X' {
            // 'X' is an unresolved SNP
            b'-'
        } else if *x.1 != '-' {
            // Other possible values are resolved SNPs (A,C,G,T)
            *x.1 as u8
        } else {
            b'-'
        }
    }).collect()
}

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {
    #[test]
    fn run_lengths() {
        use crate::format::RLE;

        let expected: Vec<RLE> = vec![
            RLE{start: 6,
                end: 33,
                matches: 28,
                mismatches : 0,
                jumps : 0,
                gap_bases: 0,
                gap_opens: 0},
            RLE{start: 82,
                end: 207,
                matches: 126,
                mismatches: 0,
                jumps: 0,
                gap_bases: 0,
                gap_opens: 0},
            RLE{start: 373,
                end: 423,
                matches: 51,
                mismatches: 0,
                jumps: 0,
                gap_bases: 0,
                gap_opens: 0},
            RLE{start: 488,
                end: 512,
                matches: 25,
                mismatches: 0,
                jumps: 0,
                gap_bases: 0,
                gap_opens: 0}];
        let input = vec!['-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M','M'];
        let got = super::run_lengths(&input);
        assert_eq!(got, expected);
    }
}
