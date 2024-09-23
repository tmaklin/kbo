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
//! Derandomizing noisy _k_-bounded matching statistics.

/// Evaluates the CDF of _k_-bounded matching statistics random match distribution.
///
/// Computes the log-probability that a matching statistic with value
/// `t` or less that is the result of mapping a _k_-mer with
/// `alphabet_size` possible characters against an index containing
/// `n_kmers` _k_-mers was generated by chance.
///
/// TODO Add formula to log_rm_max_cdf documentation
///
/// Credit to Jarno N. Alanko for deriving the random match distribution.
///
/// # Examples
/// ```rust
/// use sablast::derandomize::log_rm_max_cdf;
///
/// let alphabet_size = 4;
/// let n_kmers = 20240921;
///
/// let res = log_rm_max_cdf(10, alphabet_size, n_kmers);
/// // `res` is -4.825812199808644
/// ```
///
pub fn log_rm_max_cdf(
    t: usize,
    alphabet_size: usize,
    n_kmers: usize,
) -> f64 {
    assert!(n_kmers > 0);
    assert!(alphabet_size > 0);

    n_kmers as f64 * (- ((1.0_f64.ln() - (alphabet_size as f64).ln()).exp()).powi(t as i32 + 1)).ln_1p()
}

/// Determines a lower bound for non-random _k_-bounded matching statistic values.
///
/// Computes the probabilities that the possible values for the
/// _k_-bounded matching statistics (MS) of a _k_-mer with size `k`
/// mapped against an index with `n_kmers` total _k_-mers and
/// `alphabet_size` possible values at each character are random
/// matches. Computation terminates when the MS value that produces a
/// random match probability below `max_error_prob` is found and
/// returned.
///
/// If no MS value passes the check, the function returns `k` instead.
///
/// # Examples
/// ```rust
/// use sablast::derandomize::random_match_threshold;
///
/// let k = 31;
/// let n_kmers = 20240921;
/// let alphabet_size = 4;
/// let max_error_prob = 0.01_f64;
///
/// let threshold = random_match_threshold(k, n_kmers, alphabet_size, max_error_prob);
/// // `threshold` is 15
/// ```
pub fn random_match_threshold(
    k: usize,
    n_kmers: usize,
    alphabet_size: usize,
    max_error_prob: f64,
) -> usize {
    assert!(k > 0);
    assert!(n_kmers > 0);
    assert!(alphabet_size > 0);
    assert!(max_error_prob <= 1_f64);
    assert!(max_error_prob > 0_f64);

    for i in 1..k {
	if log_rm_max_cdf(i, alphabet_size, n_kmers) > (-max_error_prob).ln_1p() {
	    return i;
	}
    }
    k
}

/// Derandomizes a single noisy _k_-bounded matching statistic.
///
/// Derandomizes the `current_noisy_ms` matching statistic (MS) based
/// on the `next_derand_ms` value obtained from the output of this
/// function for the next noisy MS when read left-to-right, the
/// _k_-mer size `k`, and the `threshold` which specifies a lower
/// bound to consider the MS a non-random match.
///
/// Positive values of the output i64 value mean that i64 characters
/// from the beginning of the k-mer match the reference, ie. same as
/// the MS, while negative values denote distance from the last
/// character in the last _k_-mer that produced a match.
///
/// # Examples
/// ## Noisy MS has only matches
/// ```rust
/// use sablast::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 3, threshold = 2
/// //
/// // Noisy MS         : 1,2,3,3,3
/// // Derandomized MS  : 1,2,3,3,3
/// // Testing this pos :     |
///
/// let derand_ms = derandomize_ms_val(3, 3, 2, 3);
/// // `derand_ms` is 3
/// ```
///
/// ## Noisy MS has only noise
/// ```rust
/// use sablast::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 3, threshold = 2
/// //
/// // Noisy MS         :  0, 0, 2, 1,0
/// // Derandomized MS  : -4,-3,-2,-1,0
/// // Testing this pos :        |
///
/// let derand_ms = derandomize_ms_val(2, -1, 2, 3);
/// // `derand_ms` is -2
/// ```
///
/// ## Noisy MS is at beginning of a full _k_-mer match
/// ```rust
/// use sablast::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 3, threshold = 2
/// //
/// // Noisy MS         : 1,2,3, 1,2
/// // Derandomized MS  : 1,2,3,-1,0
/// // Testing this pos :     |
///
/// let derand_ms = derandomize_ms_val(3, -1, 2, 3);
/// // `derand_ms` is 3
/// ```
///
/// ## Noisy MS is at beginning of a partial _k_-mer match
/// ```rust
/// use sablast::derandomize::derandomize_ms_val;
///
/// // Parameters       : k = 4, threshold = 2
/// //
/// // Noisy MS         : 1,2,3,-1,0,1,2,3,4,4
/// // Derandomized MS  : 1,2,3,-1,0,1,2,3,4,4
/// // Testing this pos :     |
///
/// let derand_ms = derandomize_ms_val(3, -1, 2, 4);
/// // `derand_ms` is 3
/// ```
///
pub fn derandomize_ms_val(
    curr_noisy_ms: usize,
    next_derand_ms: i64,
    threshold: usize,
    k: usize,
) -> i64 {
    assert!(k > 0);
    assert!(threshold > 1);
    assert!(curr_noisy_ms <= k);
    assert!(next_derand_ms <= k as i64);

    // Default is to decrease MS by 1.
    let mut run: i64 = next_derand_ms - 1;

    if curr_noisy_ms == k {
	// Beginning of a full k-mer match
	run = k as i64;
    }

    if curr_noisy_ms > threshold && next_derand_ms < curr_noisy_ms as i64 {
	// Beginning of a partial k-mer match
	// Only useful if threshold > 1 and k > 3
	run = curr_noisy_ms as i64;
    }

    run
}

/// Derandomizes a sequence of noisy _k_-bounded matching statistics.
///
/// Iterates over a sequence of noisy _k_bounded matching statistics
/// `ms` in reverse to identify values that are the result of random
/// matching between _k_-mers of size `k` and an index that the lower
/// bound `threshold` was calculated for.
///
/// # Examples
/// ```rust
/// use sablast::derandomize::derandomize_ms_vec;
///
/// let k = 3;
/// let threshold = 2;
/// let noisy_ms = vec![1,2,2,3,2,2,3,2,1,2,3,1,1,1,2,3,1,2];
///
/// let derand_ms = derandomize_ms_vec(&noisy_ms, k, threshold);
/// // `derand_ms` has [0,1,2,3,1,2,3,0,1,2,3,-1,0,1,2,3,-1,0]
/// ```
///
pub fn derandomize_ms_vec(
    noisy_ms: &[usize],
    k: usize,
    threshold: usize,
) -> Vec<i64> {
    assert!(k > 0);
    assert!(threshold > 1);
    assert!(noisy_ms.len() > 2);

    let len = noisy_ms.len();
    let mut derand_ms: Vec<i64> = vec![0; len];

    // Traverse the matching statistics in reverse.
    derand_ms[len - 1] = if noisy_ms[len - 1] > threshold { noisy_ms[len - 1]} else { 0 } as i64;
    for i in 2..len {
	derand_ms[len - i] = derandomize_ms_val(noisy_ms[len - i], derand_ms[len - i + 1], threshold, k);
    }

    derand_ms
}

////////////////////////////////////////////////////////////////////////////////
// Tests
//
#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn log_rm_max_cdf() {
	let expected = vec![-1306319.1078024083,-318761.2492719044,-79220.9269610741,-19776.1823255263,-4942.2344281681,-1235.4454790664,-308.8543003470,-77.2131332649,-19.3032557026,-4.8258121998,-1.2064529421,-0.3016132288,-0.0754033068,-0.0188508267,-0.0047127067,-0.0011781767,-0.0002945442,-0.0000736360,-0.0000184090,-0.0000046023,-0.0000011506,-0.0000002876,-0.0000000719,-0.0000000180,-0.0000000045,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000];
	let alphabet_size = 4;
	let n_kmers = 20240921;
	let k = 1..32;
	k.for_each(|t| assert_approx_eq!(super::log_rm_max_cdf(t, alphabet_size, n_kmers), expected[t - 1], 1e-8f64));
    }

    #[test]
    fn random_match_threshold() {
	let expected = vec![15,18,22,25,28];
	let alphabet_size = 4;
	let n_kmers = 20240921;
	let k = 31;
	let factor = 1..6;
	factor.for_each(|i| assert_eq!(super::random_match_threshold(k, n_kmers, alphabet_size, (0.01_f64).powf(i as f64)), expected[i - 1]));
    }

    #[test]
    fn derandomize_ms_val_full_match() {
	// Parameters       : k = 3, threshold = 2
	//
	// Noisy MS         : 1,2,3,3,3
	// Derandomized MS  : 1,2,3,3,3
	// Testing this pos :     |

	let expected = 3;
	let got = super::derandomize_ms_val(3, 3, 2, 3);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_val_only_noise() {
	// Parameters       : k = 3, threshold = 2
	//
	// Noisy MS         :  0, 0, 2, 1,0
	// Derandomized MS  : -4,-3,-2,-1,0
	// Testing this pos :        |

	let expected = -2;
	let got = super::derandomize_ms_val(2, -1, 2, 3);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_val_beginning_of_full_match() {
	// Parameters       : k = 3, threshold = 2
	//
	// Noisy MS         : 1,2,3, 1,2
	// Derandomized MS  : 1,2,3,-1,0
	// Testing this pos :     |

	let expected = 3;
	let got = super::derandomize_ms_val(3, -1, 2, 3);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_val_beginning_of_partial_match() {
	// Parameters       : k = 4, threshold = 2
	//
	// Noisy MS         : 1,2,3,-1,0,1,2,3,4,4
	// Derandomized MS  : 1,2,3,-1,0,1,2,3,4,4
	// Testing this pos :     |

	let expected = 3;
	let got = super::derandomize_ms_val(3, -1, 2, 4);

	assert_eq!(got, expected);
    }

    #[test]
    fn derandomize_ms_vec() {
	let noisy_ms = vec![1,2,2,3,2,2,3,2,1,2,3,1,1,1,2,3,1,2];
	let expected = vec![0,1,2,3,1,2,3,0,1,2,3,-1,0,1,2,3,-1,0];
	let got = super::derandomize_ms_vec(&noisy_ms, 3, 2);

	assert_eq!(got, expected);
    }
}
