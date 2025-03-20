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

//! kbo is an approximate local aligner based on converting [_k_-bounded matching
//! statistics](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1)
//! into a character representation of the underlying alignment sequence.
//!
//! Currently, kbo supports two main operations:
//!
//! - `kbo find` [matches](matches()) the _k_-mers in a query sequence with the
//!   reference and reports the local alignment segments found within the
//!   reference. Find is useful for problems that can be solved with
//!   [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
//! - `kbo map` [maps](map()) the query sequence against a reference
//!   sequence, and reports the nucleotide sequence of the alignment relative to
//!   the reference. Map solves the same problem as
//!   [snippy](https://github.com/tseemann/snippy) and [ska
//!   map](https://docs.rs/ska/latest/ska/#ska-map).
//!
//! kbo uses the [Spectral Burrows-Wheeler
//! Transform](https://docs.rs/sbwt/latest/sbwt/) data structure that allows
//! efficient _k_-mer matching between a target and a query sequence and
//! fast retrieval of the _k_-bounded matching statistic for each _k_-mer match.
//!
//! # Installing the kbo executable
//! See installation instructions at [GitHub](https://github.com/tmaklin/kbo-cli).
//!
//! # Usage
//!
//! kbo can be run directly on fasta files without an initial indexing step.
//! Prebuilt indexes are supported via `kbo build` but are only
//! relevant in `kbo find` analyses where the reference _k_-mers can be
//! concatenated into a single contig.
//!
//! kbo can read inputs compressed in the DEFLATE format (gzip, zlib, etc.).
//! bzip2 and xz support can be enabled by adding the "bzip2" and "xz" feature
//! flags to [needletail](https://docs.rs/needletail) in the kbo Cargo.toml.
//!
//! ## kbo find
//!
//! To set up the example, download the fasta sequence of the [_Escherichia
//! coli_ Nissle
//! 1917](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000714595.1/) genome
//! from the NCBI and the [pks
//! island](https://raw.githubusercontent.com/tmaklin/clbtype/refs/heads/main/db/db.fasta)
//! gene sequences from GitHub. Example output was generated with versions
//! ASM71459v1 and rev 43bbd36.
//!
//! ### Find gene sequence locations
//! In the directory containing the input files, run
//! ```text
//! kbo find --reference db.fasta GCF_000714595.1_ASM71459v1_genomic.fna
//! ```
//! <details>
//! <summary>
//! This will produce the output (click to expand)
//! </summary>
//!
//! |query|ref|q.start|q.end|strand|length|mismatches|query.contig|ref.contig|
//! |-----|---|-------|-----|------|------|----------|------------|----------|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|226708|227226|+|519|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|2289596|2290543|+|949|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3152039|3161660|+|9623|1|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3161701|3164301|+|2601|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3164311|3165180|+|870|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3165210|3165458|+|249|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3165462|3167857|+|2397|1|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3167905|3172701|+|4797|2|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3172751|3175783|+|3033|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3175827|3182327|+|6501|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3182338|3190258|+|7922|1|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3190320|3196124|+|5807|1|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3196155|3198614|+|2460|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3198627|3200856|+|2231|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3200891|3201403|+|513|1|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|4502887|4503405|+|519|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5145962|5147210|+|1249|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5147272|5149449|+|2179|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5351015|5351533|+|519|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5352280|5352503|+|224|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5354674|5356713|+|2040|1|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5381795|5381945|+|151|0|db.fasta|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//!
//! </summary>
//! </details>
//!
//! ### Find gene sequence locations with names
//! If you need to know which gene in db.fasta the matches are for, add the `--detailed` toggle:
//! ```text
//! kbo find --detailed --reference db.fasta GCF_000714595.1_ASM71459v1_genomic.fna
//! ```
//!
//! <details>
//! <summary>
//! This replaces the query.contig column with the name of the contig (click to expand)
//! </summary>
//!
//! |query|ref|q.start|q.end|strand|length|mismatches|query.contig|ref.contig|
//! |-----|---|-------|-----|------|------|----------|------------|----------|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|226708|227226|+|519|0|clbS-like_4ce09a|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|2289596|2289808|+|213|0|clbR locus_tag=ECOK1_RS11410 product="colibactin biosynthesis LuxR family transcriptional regulator ClbR" protein_id=WP_000357141.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|2289809|2290543|+|735|0|clbA locus_tag=ECOK1_RS11415 product="colibactin biosynthesis phosphopantetheinyl transferase ClbA" protein_id=WP_001217110.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3152039|3161660|+|9623|1|clbB locus_tag=ECOK1_RS11405 product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbB" protein_id=WP_001518711.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3161701|3164301|+|2601|0|clbC locus_tag=ECOK1_RS11400 product="colibactin polyketide synthase ClbC" protein_id=WP_001297908.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3164311|3165180|+|870|0|clbD locus_tag=ECOK1_RS11395 product="colibactin biosynthesis dehydrogenase ClbD" protein_id=WP_000982270.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3165210|3165458|+|249|0|clbE locus_tag=ECOK1_RS11390 product="colibactin biosynthesis aminomalonyl-acyl carrier protein ClbE" protein_id=WP_001297917.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3165462|3166592|+|1131|0|clbF locus_tag=ECOK1_RS11385 product="colibactin biosynthesis dehydrogenase ClbF" protein_id=WP_000337350.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3166589|3167857|+|1269|1|clbG locus_tag=ECOK1_RS11380 product="colibactin biosynthesis acyltransferase ClbG" protein_id=WP_000159201.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3167905|3172701|+|4797|2|clbH locus_tag=ECOK1_RS11375 product="colibactin non-ribosomal peptide synthetase ClbH" protein_id=WP_001304254.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3172751|3175783|+|3033|0|clbI locus_tag=ECOK1_RS11370 product="colibactin polyketide synthase ClbI" protein_id=WP_000829570.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3175827|3182327|+|6501|0|clbJ locus_tag=ECOK1_RS11365 product="colibactin non-ribosomal peptide synthetase ClbJ" protein_id=WP_001468003.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3180417|3181703|+|1287|2|clbK locus_tag=ECOK1_RS11360 product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK" protein_id=WP_000222467.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3182338|3188802|+|6465|2|clbK locus_tag=ECOK1_RS11360 product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK" protein_id=WP_000222467.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3186070|3187356|+|1287|1|clbJ locus_tag=ECOK1_RS11365 product="colibactin non-ribosomal peptide synthetase ClbJ" protein_id=WP_001468003.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3188795|3190258|+|1464|0|clbL locus_tag=ECOK1_RS11355 product="colibactin biosynthesis amidase ClbL" protein_id=WP_001297937.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3190320|3191759|+|1440|0|clbM locus_tag=ECOK1_RS11350 product="precolibactin export MATE transporter ClbM" protein_id=WP_000217768.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3191756|3196124|+|4370|1|clbN locus_tag=ECOK1_RS11345 product="colibactin non-ribosomal peptide synthetase ClbN" protein_id=WP_001327259.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3196155|3198614|+|2460|0|clbO locus_tag=ECOK1_RS11340 product="colibactin polyketide synthase ClbO" protein_id=WP_001029878.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3198627|3200141|+|1515|0|clbP locus_tag=ECOK1_RS11335 product="precolibactin peptidase ClbP" protein_id=WP_002430641.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3200134|3200856|+|723|0|clbQ locus_tag=ECOK1_RS11330 product="colibactin biosynthesis thioesterase ClbQ" protein_id=WP_000065646.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|3200891|3201403|+|513|1|clbS locus_tag=ECOK1_RS11325 product="colibactin self-protection protein ClbS" protein_id=WP_000290498.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|4502887|4503405|+|519|0|clbS-like_4ce09a|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5145962|5147210|+|1249|0|clbL locus_tag=ECOK1_RS11355 product="colibactin biosynthesis amidase ClbL" protein_id=WP_001297937.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5147272|5148479|+|1208|0|clbM locus_tag=ECOK1_RS11350 product="precolibactin export MATE transporter ClbM" protein_id=WP_000217768.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5148478|5149449|+|972|0|clbN locus_tag=ECOK1_RS11345 product="colibactin non-ribosomal peptide synthetase ClbN" protein_id=WP_001327259.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5351015|5351533|+|519|0|clbS-like_4ce09a|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5352280|5352503|+|224|0|clbS-like_4ce09a|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5354674|5356713|+|2040|1|clbN locus_tag=ECOK1_RS11345 product="colibactin non-ribosomal peptide synthetase ClbN" protein_id=WP_001327259.1|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//! |GCF_000714595.1_ASM71459v1_genomic.fna|db.fasta|5381795|5381945|+|151|0|clbS-like_4ce09a|NZ_CP007799.1 Escherichia coli Nissle 1917 chromosome, complete genome|
//!
//! </details>
//!
//! Note that the current implementation `--detailed` significantly slows down
//! the algorithm. Future versions of kbo may address this by incorporating
//! colors in the index structure.
//!
//! ### Find containment of gene sequences in assembly
//! Alternatively, if you are only interested in whether the contigs in `db.fasta` are present in the assembly, run
//! ```text
//! kbo find --reference GCF_000714595.1_ASM71459v1_genomic.fna db.fasta
//! ```
//!
//! <details>
//! <summary>
//! which will return (click to expand)
//! </summary>
//!
//! |query|ref|q.start|q.end|strand|length|mismatches|query.contig|ref.contig|
//! |-----|---|-------|-----|------|------|----------|------------|----------|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|513|+|513|1|GCF_000714595.1_ASM71459v1_genomic.fna|clbS\|locus_tag=ECOK1_RS11325\|product="colibactin self-protection protein ClbS"\|protein_id=WP_000290498.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|723|+|723|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbQ\|locus_tag=ECOK1_RS11330\|product="colibactin biosynthesis thioesterase ClbQ"\|protein_id=WP_000065646.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1515|+|1515|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbP\|locus_tag=ECOK1_RS11335\|product="precolibactin peptidase ClbP"\|protein_id=WP_002430641.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|2460|+|2460|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbO\|locus_tag=ECOK1_RS11340\|product="colibactin polyketide synthase ClbO"\|protein_id=WP_001029878.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|4368|+|4369|1|GCF_000714595.1_ASM71459v1_genomic.fna|clbN\|locus_tag=ECOK1_RS11345\|product="colibactin non-ribosomal peptide synthetase ClbN"\|protein_id=WP_001327259.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1208|+|1208|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbM\|locus_tag=ECOK1_RS11350\|product="precolibactin export MATE transporter ClbM"\|protein_id=WP_000217768.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1440|+|1440|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbM\|locus_tag=ECOK1_RS11350\|product="precolibactin export MATE transporter ClbM"\|protein_id=WP_000217768.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1464|+|1464|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbL\|locus_tag=ECOK1_RS11355\|product="colibactin biosynthesis amidase ClbL"\|protein_id=WP_001297937.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|6465|+|6465|2|GCF_000714595.1_ASM71459v1_genomic.fna|clbK\|locus_tag=ECOK1_RS11360\|product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbK"\|protein_id=WP_000222467.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|6501|+|6501|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbJ\|locus_tag=ECOK1_RS11365\|product="colibactin non-ribosomal peptide synthetase ClbJ"\|protein_id=WP_001468003.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|3033|+|3033|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbI\|locus_tag=ECOK1_RS11370\|product="colibactin polyketide synthase ClbI"\|protein_id=WP_000829570.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|4797|+|4797|2|GCF_000714595.1_ASM71459v1_genomic.fna|clbH\|locus_tag=ECOK1_RS11375\|product="colibactin non-ribosomal peptide synthetase ClbH"\|protein_id=WP_001304254.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1269|+|1269|1|GCF_000714595.1_ASM71459v1_genomic.fna|clbG\|locus_tag=ECOK1_RS11380\|product="colibactin biosynthesis acyltransferase ClbG"\|protein_id=WP_000159201.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|1131|+|1131|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbF\|locus_tag=ECOK1_RS11385\|product="colibactin biosynthesis dehydrogenase ClbF"\|protein_id=WP_000337350.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|249|+|249|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbE\|locus_tag=ECOK1_RS11390\|product="colibactin biosynthesis aminomalonyl-acyl carrier protein ClbE"\|protein_id=WP_001297917.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|870|+|870|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbD\|locus_tag=ECOK1_RS11395\|product="colibactin biosynthesis dehydrogenase ClbD"\|protein_id=WP_000982270.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|2601|+|2601|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbC\|locus_tag=ECOK1_RS11400\|product="colibactin polyketide synthase ClbC"\|protein_id=WP_001297908.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|9621|+|9622|1|GCF_000714595.1_ASM71459v1_genomic.fna|clbB\|locus_tag=ECOK1_RS11405\|product="colibactin hybrid non-ribosomal peptide synthetase/type I polyketide synthase ClbB"\|protein_id=WP_001518711.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|213|+|213|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbR\|locus_tag=ECOK1_RS11410\|product="colibactin biosynthesis LuxR family transcriptional regulator ClbR"\|protein_id=WP_000357141.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|735|+|735|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbA\|locus_tag=ECOK1_RS11415\|product="colibactin biosynthesis phosphopantetheinyl transferase ClbA"\|protein_id=WP_001217110.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|519|+|519|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbS-like_4ce09a|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1|519|+|519|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbS-like_4ce09a|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|216|1464|+|1249|0|GCF_000714595.1_ASM71459v1_genomic.fna|clbL\|locus_tag=ECOK1_RS11355\|product="colibactin biosynthesis amidase ClbL"\|protein_id=WP_001297937.1|
//! |db.fasta|GCF_000714595.1_ASM71459v1_genomic.fna|1156|4167|+|3013|1|GCF_000714595.1_ASM71459v1_genomic.fna|clbN\|locus_tag=ECOK1_RS11345\|product="colibactin non-ribosomal peptide synthetase ClbN"\|protein_id=WP_001327259.1|
//!
//! </details>
//!
//! ## kbo map
//!
//! kbo map can be used to align a query sequence against a reference sequence.
//! This is useful in for example generating a reference-based alignment of
//! multiple related genomes against a good reference assembly.
//!
//! To run this example, download the genome sequence of the [_E. coli_
//! UTI89](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013265.1/) strain
//! from the NCBI (ASM1326v1).
//!
//! ### Reference-based alignment
//! Run
//! ```text
//! kbo map --reference GCF_000714595.1_ASM71459v1_genomic.fna GCF_000013265.1_ASM1326v1_genomic.fna > result.aln
//! ```
//!
//! which will write the alignment sequence to `result.aln`. Note that kbo map
//! always writes to stdout.
//!
//! If you have multiple sequences you need to align, either supply them as
//! arguments to `kbo map` or process them using gnu parallel:
//!
//! ```text
//! parallel -j 'kbo map --reference GCF_000714595.1_ASM71459v1_genomic.fna {}' < query_paths.txt > result.aln
//! ```
//!
//! kbo map also accepts the `--threads` argument to parallelise either the
//! index construction (in the case of a single query), or run in parallel over
//! the input files (multiple queries).
//!

#![warn(missing_docs,
        missing_debug_implementations, missing_copy_implementations,
        trivial_casts, trivial_numeric_casts,
        unsafe_code,
        unstable_features,
        unused_import_braces, unused_qualifications)]

use sbwt::SbwtIndexVariant;

pub mod derandomize;
pub mod format;
pub mod index;
pub mod translate;
pub mod variant_calling;

/// Options and parameters for [find]
#[non_exhaustive]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FindOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,
    /// Maximum length of a single gap segment before splitting an alignment.
    pub max_gap_len: usize,
}

impl Default for FindOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::FindOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// opts.max_gap_len = 0;
    /// # let expected = kbo::FindOpts::default();
    /// # assert_eq!(opts, expected);
    /// ```
    ///
    fn default() -> FindOpts {
        FindOpts {
            max_error_prob: 0.0000001,
            max_gap_len: 0,
        }
    }
}

/// Options and parameters for [matches](matches())
#[non_exhaustive]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MatchOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,
}

impl Default for MatchOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::MatchOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// # let expected = kbo::MatchOpts::default();
    /// # assert_eq!(opts, expected);
    /// ```
    ///
    fn default() -> MatchOpts {
        MatchOpts {
            max_error_prob: 0.0000001,
        }
    }
}

/// Options and parameters for [map]
#[non_exhaustive]
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MapOpts {
    /// Prefix match lengths with probability higher than `max_error_prob` to
    /// happen at random are considered noise.
    pub max_error_prob: f64,
}

impl Default for MapOpts {
    /// Default to these values:
    /// ```rust
    /// let mut opts = kbo::MapOpts::default();
    /// opts.max_error_prob = 0.0000001;
    /// # let expected = kbo::MapOpts::default();
    /// # assert_eq!(opts, expected);
    /// ```
    ///
    fn default() -> MapOpts {
        MapOpts {
            max_error_prob: 0.0000001,
        }
    }
}

/// Builds an SBWT index from some fasta or fastq files.
///
/// Reads all sequence data in `seq_files` and builds an SBWT index
/// with the parameters and resources specified in `build_opts` (see
/// [index::BuildOpts] for details).
///
/// Prebuilt indexes can currently only be used with kbo find.
///
/// All files and sequence data in `seq_files` are merged into the
/// same index. It is not possible extract the individual sequences
/// from the index after it has been built; run `kbo map -r
/// <query_file> <seq_files>` if you need to know which reference
/// sequences the alignments are for.
///
/// Returns a tuple containing the built
/// [SbwtIndexVariant](https://docs.rs/sbwt/latest/sbwt/enum.SbwtIndexVariant.html)
/// and
/// [sbwt::LcsArray](https://docs.rs/sbwt/latest/sbwt/struct.LcsArray.html).
///
/// Panics if a file in `seq_files` is not readable or a valid FASTX
/// file.
///
/// # Examples
/// ```rust
/// use kbo::build;
/// use kbo::index::BuildOpts;
///
/// let inputs: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
///
/// let opts = BuildOpts::default();
/// let (sbwt_index, lcs_array) = build(&inputs, opts);
/// ```
///
pub fn build(
    seq_data: &[Vec<u8>],
    build_opts: index::BuildOpts,
) -> (SbwtIndexVariant, sbwt::LcsArray) {
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
/// use kbo::build;
/// use kbo::matches;
/// use kbo::index::BuildOpts;
/// use kbo::MatchOpts;
///
/// let reference: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// let (sbwt, lcs) = build(&reference, opts);
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let ms_vectors = matches(&query, &sbwt, &lcs, MatchOpts::default());
/// // `ms_vectors` has ['-','-','-','-','-','-','-','-','-','M','M','M','-','-']
/// # assert_eq!(ms_vectors, vec!['-','-','-','-','-','-','-','-','-','M','M','M','-','-']);
/// ```
///
pub fn matches(
    query_seq: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
    match_opts: MatchOpts,
) -> Vec<char> {
    let (k, threshold) = match sbwt {
        SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
            (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, match_opts.max_error_prob))
        },
    };

    let noisy_ms: Vec<usize> = index::query_sbwt(query_seq, sbwt, lcs).iter().map(|x| x.0).collect();
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
/// use kbo::build;
/// use kbo::map;
/// use kbo::index::BuildOpts;
/// use kbo::MapOpts;
///
/// let query: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// opts.build_select = true;
/// let (sbwt_query, lcs_query) = build(&query, opts);
///
/// let reference = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let alignment = map(&reference, &sbwt_query, &lcs_query, MapOpts::default());
/// // `ms_vectors` has [45,45,45,45,45,45,45,45,45,65,71,71,45,45]
/// # assert_eq!(alignment, vec![45,45,45,45,45,45,45,45,45,65,71,71,45,45]);
/// ```
///
pub fn map(
    ref_seq: &[u8],
    query_sbwt: &SbwtIndexVariant,
    query_lcs: &sbwt::LcsArray,
    map_opts: MapOpts,
) -> Vec<u8> {
    let (k, threshold) = match query_sbwt {
        SbwtIndexVariant::SubsetMatrix(ref sbwt) => {
            (sbwt.k(), derandomize::random_match_threshold(sbwt.k(), sbwt.n_kmers(), 4_usize, map_opts.max_error_prob))
        },
    };

    let noisy_ms = index::query_sbwt(ref_seq, query_sbwt, query_lcs);
    let derand_ms = derandomize::derandomize_ms_vec(&noisy_ms.iter().map(|x| x.0).collect::<Vec<usize>>(), k, threshold);

    let translation = translate::translate_ms_vec(&derand_ms, k, threshold);
    let refined = translate::refine_translation(&translation, &noisy_ms, query_sbwt, threshold);

    format::relative_to_ref(ref_seq, &refined)
}

/// Finds the _k_-mers from an SBWT index in a query fasta or fastq file.
///
/// Aligns the sequence data in `query_seq` against the SBWT index
/// `sbwt` and its LCS array `lcs` using [matches][matches()]. Then uses
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
///
/// TODO Add better examples to find()
///
/// ```rust
/// use kbo::build;
/// use kbo::find;
/// use kbo::FindOpts;
/// use kbo::index::BuildOpts;
/// use kbo::format::RLE;
///
/// let reference: Vec<Vec<u8>> = vec![vec![b'A',b'A',b'A',b'G',b'A',b'A',b'C',b'C',b'A',b'-',b'T',b'C',b'A',b'G',b'G',b'G',b'C',b'G']];
/// let mut opts = BuildOpts::default();
/// opts.k = 3;
/// let (sbwt, lcs) = build(&reference, opts);
///
/// let query = vec![b'G',b'T',b'G',b'A',b'C',b'T',b'A',b'T',b'G',b'A',b'G',b'G',b'A',b'T'];
///
/// let local_alignments = find(&query, &sbwt, &lcs, FindOpts::default());
/// // `local_alignments` has [(10, 12, 3, 0)]
/// # assert_eq!(local_alignments, vec![RLE{start: 10, end: 12, matches: 3, mismatches: 0, jumps: 0, gap_bases: 0, gap_opens: 0}]);
/// ```
///
pub fn find(
    query_seq: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
    find_opts: FindOpts,
) -> Vec<format::RLE> {
    let match_opts = MatchOpts { max_error_prob: find_opts.max_error_prob };
    let aln = matches(query_seq, sbwt, lcs, match_opts);
    if find_opts.max_gap_len > 0 {
        format::run_lengths_gapped(&aln, find_opts.max_gap_len)
    } else {
        format::run_lengths(&aln)
    }
}
