---
title: "PlainMap: a simple mapping pipeline for ancient and modern DNA"
tags:
  - Bash
  - bioinformatics
  - read mapping
  - ancient DNA
  - FASTQ
authors:
  - name: Michael V. Wesbury
    affiliation: 1
    corresponding: true
affiliations:
  - name: Department of Health Technology, Section for Bioinformatics, Technical University of Denmark, Kongens Lyngby, Denmark
    index: 1
date: 07 July 2026
bibliography: paper.bib
---

# Summary

PlainMap is a lightweight and restartable read mapping pipeline for ancient and modern DNA sequencing data. Implemented as a single Bash script, the software focuses specifically on preprocessing, mapping, duplicate removal, and BAM generation while avoiding workflow engines and complex configuration systems. PlainMap supports heterogeneous sequencing datasets, including mixed single-end and paired-end libraries, multiple sequencing platforms, and ambiguous SRA-exported FASTQ formats. The pipeline additionally provides deterministic chunking, restartable execution, conservative validation procedures, and optional pilot subsampling of raw fragments for empirical comparison of mapping strategies.

![PlainMap workflow overview. Heterogeneous FASTQ inputs are classified by header content before optional pre-trimming pilot subsampling and adaptive chunking, preprocessing with fastp, mode-specific BWA mapping, and conversion into duplicate-filtered BAM files and summary statistics. Optional outputs include mapDamage analyses and FASTQ files containing unmapped reads.](figure1.png)

# Statement of Need

Read mapping is a foundational step in next-generation sequencing analysis, but practical implementation frequently requires multiple preprocessing and optimization steps, including adapter trimming, read merging, filtering, duplicate removal, and mapper parameter tuning. These decisions become particularly important when working with ancient and historical DNA, where short fragment lengths, post-mortem damage, and divergence from available reference genomes can influence mapping performance and make optimal alignment strategies difficult to predict a priori. In addition, many sequencing projects involve heterogeneous datasets containing mixed short-read sequencing platforms, inconsistent FASTQ naming conventions, and combinations of single-end and paired-end libraries. Together, these factors often necessitate empirical testing of alternative mapping approaches and substantial preprocessing before mapping can begin. At the same time, increasing sequencing throughput has led to larger datasets and longer mapping runtimes, making restartable execution increasingly important for efficient exploratory analyses.

Several established frameworks support read mapping for ancient and modern DNA, including PALEOMIX [@schubert2014paleomix], nf-core/eager [@fellows_yates2021eager], GENERODE [@kutschera2022generode], and Mapache [@neuenschwander2023mapache]. These frameworks provide extensive functionality and scalability, frequently integrating downstream analyses such as contamination estimation, genotyping, or population genetic inference. While well suited to large-scale and standardized analyses, they typically rely on workflow engines and structured configuration systems that introduce additional layers of abstraction between users and pipeline execution. For researchers performing exploratory analyses, testing alternative mapping strategies, or working with heterogeneous datasets, this complexity can complicate debugging, obscure failure modes, and increase the overhead associated with rapid methodological iteration. This creates a need for lightweight mapping workflows that minimize configuration overhead, support heterogeneous datasets, facilitate rapid testing of alternative mapping strategies, and remain compatible with downstream project-specific analytical pipelines.

# Software Design

PlainMap was developed to provide a robust mapping workflow that minimizes configuration complexity while supporting heterogeneous sequencing datasets and restartable execution. The pipeline adopts a deliberately narrow scope focused on preprocessing, read mapping, duplicate removal, and BAM generation. Adapter trimming, quality filtering, polyG trimming, and read merging are performed using fastp [@chen2018fastp]. Mapping is conducted using the Burrows-Wheeler Aligner (BWA) [@li2009bwa], with alternative alignment strategies implemented for modern and degraded DNA datasets. Modern mode uses `bwa mem` for both single-end and paired-end data. Ancient-fast mode uses `bwa aln` followed by `bwa samse`, mapping merged and single-end-like fragments as single-end reads. Ancient-all mode uses `bwa aln` followed by `bwa sampe` for paired-end reads and `bwa samse` for single-end-like reads, thereby retaining both paired-end and single-end information.

Ancient-fast was developed for highly degraded datasets where authentic endogenous fragments are expected to be sufficiently short that most paired-end reads merge during preprocessing. By excluding unmerged paired-end reads from alignment, this mode reduces computational requirements and may reduce the contribution of longer fragments that are more likely to originate from contamination or exogenous sources. Ancient-all retains all available reads and may therefore recover additional endogenous fragments in datasets where substantial numbers of authentic paired-end reads remain unmerged. Both ancient-fast and ancient-all use alignment parameters intended for short and damaged DNA fragments [@schubert2012ancientmapping]. BAM processing, sorting, duplicate removal, indexing, and summary statistic generation are performed using SAMtools [@li2009samtools]. Optional ancient DNA damage profiling is implemented through integration with mapDamage2 [@jonsson2013mapdamage].

PlainMap performs FASTQ header-based pairing inference rather than relying on filename conventions, allowing mixed single-end and paired-end datasets from multiple sequencing platforms and public repositories to be processed within a single workflow. Conservative validation procedures are applied throughout execution, including gzip integrity checks, pairing consistency checks, and BAM validation, with the pipeline terminating when ambiguous or inconsistent inputs are detected.

To support exploratory analyses and large sequencing datasets, PlainMap incorporates deterministic FASTQ chunking, checkpoint-based restartability, and optional pilot subsampling of raw fragments. These features enable rapid comparison of mapping strategies on identical input data while allowing interrupted analyses to resume without recomputing completed work.

# Research Impact Statement

PlainMap was developed to support the analysis of ancient, historical, and modern sequencing datasets, particularly those containing heterogeneous combinations of sequencing platforms, library types, and FASTQ formats. The software is intended for researchers who require a lightweight and robust mapping workflow that can be deployed with minimal configuration and readily integrated into project-specific downstream analyses.

PlainMap is particularly well suited to exploratory sequencing projects in which mapping strategies, reference genomes, or library characteristics are uncertain and require empirical evaluation. The inclusion of pilot subsampling, multiple alignment modes, and restartable execution enables rapid testing of alternative approaches while reducing the computational cost of interrupted or repeated analyses.

By focusing exclusively on preprocessing, mapping, and BAM generation, PlainMap provides a lightweight alternative to larger workflow-based frameworks for users who require a lightweight mapping workflow while retaining flexibility over downstream analytical decisions.

# AI Usage Disclosure

Generative AI tools were used throughout the design, implementation, documentation, testing, and manuscript preparation of PlainMap. AI-assisted interactions contributed suggestions for software architecture, implementation details, debugging, edge-case identification, error handling, code refinement, and technical writing.

The author retained full responsibility for all software design decisions, validation procedures, and final code. AI-generated suggestions were treated as advisory and were critically assessed before adoption. All code, analyses, and manuscript text were manually reviewed, tested, revised, and verified by the author prior to inclusion in the final software release and manuscript.
