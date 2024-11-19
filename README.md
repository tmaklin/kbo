# kbo
Spectral Burrows-Wheeler transform accelerated local alignment search.

kbo is an approximate local aligner based on converting [_k_-bounded matching
statistics](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1)
into a character representation of the underlying alignment sequence.

Documentation is available at [https://docs.rs/kbo](https://docs.rs/kbo).

## Installation
kbo is distributed as three separate Rust packages:
- [kbo](https://github.com/tmaklin/kbo) contains a Rust library implementing the core algorithm (this repository).
- [kbo-cli](https://github.com/tmaklin/kbo-cli) provides a command-line interface for `kbo find` and `kbo map`.
- [kbo-gui](https://github.com/tmaklin/kbo-gui) is a work in progress WebAssembly graphical user interface for running kbo in the browser.

Precompiled binaries for kbo-cli are available from the [Releases page](https://github.com/tmaklin/kbo-cli/releases).

## About
kbo supports two main operations:

- `kbo find` matches the _k_-mers in a query sequence with the
  reference and reports the local alignment segments found within the
  reference. Find is useful for problems that can be solved with
  [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
- `kbo map` maps the query sequence against a reference
  sequence, and reports the nucleotide sequence of the alignment relative to
  the reference. Map solves the same problem as
  [snippy](https://github.com/tseemann/snippy) and [ska
  map](https://docs.rs/ska/latest/ska/#ska-map).

kbo uses the [Spectral Burrows-Wheeler
Transform](https://docs.rs/sbwt/latest/sbwt/) data structure that allows
efficient _k_-mer matching between a target and a query sequence and
fast retrieval of the _k_-bounded matching statistic for each _k_-mer match.

## License
kbo is dual-licensed under the [MIT](LICENSE-MIT) and [Apache 2.0](LICENSE-APACHE) licenses.
