---
title: 'Capra: Scalable HEALPix-native pipeline for reconstructing all-sky maps from event-counting spacecraft data'
tags:
  - Mathematica
  - Python
  - astronomy
  - heliophysics
  - energetic neutral atoms
  - HEALPix
  - image reconstruction
  - sky maps
authors:
  - name: Nikola Bukowiecka
    orcid: "0000-0003-1129-5310"
    affiliation: Department of Physics, University of Rhode Island, Kingston, RI, USA
date: 06/07/26
bibliography: paper.bib
---

# Summary

Capra is a software pipeline for reconstructing all-sky intensity maps from event-counting spacecraft observations. It was developed for energetic neutral atom (ENA) measurements from the Interstellar Boundary Explorer (IBEX) mission [@mccomas2009ibex; @funsten2009ibexhi] and is designed with future Interstellar Mapping and Acceleration Probe (IMAP) applications in mind [@mccomas2018imap]. ENA maps provide global views of the outer heliosphere, but the raw observations are not collected on a regular image grid. Instead, the spacecraft scans the sky in repeated arcs, producing binned counts, exposure times, and background estimates that must be converted into physically interpretable sky maps.

Capra performs this reconstruction directly on the sphere using the HEALPix tessellation [@gorski2005healpix]. For each energy channel, it produces maps of counts, exposure, background-subtracted signal, signal rate, and ENA intensity. The pipeline supports two reconstruction products: a baseline boresight-assigned map, which assigns each scan-angle bin to the corresponding HEALPix direction, and an optional smoothed map, in which counts, exposure, and signal are redistributed over neighboring HEALPix pixels using normalized kernel weights. The current implementation includes a Gaussian smoothing option and is structured so that alternative kernels can be introduced by the user.

Capra also includes postprocessing and comparison tools. These include detector-efficiency correction, relative normalization, globally distributed flux (GDF) normalization, transformations between HEALPix and rectangular comparison grids, coordinate-frame transformations, HEALPix tessellation changes, plotting utilities, and diagnostic checks for map consistency. The software is intended for transparent, reproducible ENA map-making and for comparison between reconstruction products such as elementary IBEX maps, Capra HEALPix maps, and THESEUS maps [@osthus2024theseus].

# Statement of need

All-sky ENA imaging is central to observational studies of the heliosphere because ENAs carry information from distant plasma regions back to near-Earth spacecraft without being tied to magnetic field lines. However, reconstructing an ENA sky map from event-counting data is not simply a plotting problem. The observations are collected with nonuniform sky coverage, finite exposure, background contributions, instrument response, spacecraft spin geometry, and repeated scan arcs. A useful reconstruction pipeline must therefore preserve the relationship between counts, exposure, signal, rate, and intensity while also producing maps that can be compared across methods, energy channels, coordinate frames, and angular resolutions.

Capra addresses this need by providing a HEALPix-native reconstruction workflow for IBEX-Hi-like orbit-arc data. Its target users are heliophysics researchers working with ENA sky maps, mission-data analysts comparing map-making approaches, and researchers preparing high-resolution workflows for IMAP. The software is designed to make the reconstruction assumptions explicit: the baseline product follows the measurement geometry as directly as possible, while smoothing is an optional user-controlled post-reconstruction choice rather than an implicit part of the data model.

The pipeline was developed as part of a broader scientific study applying Capra to IBEX-Hi Map2018B orbit-arc data and comparing the resulting maps with elementary IBEX maps and THESEUS products. That companion study, written with Daniel B. Reisenfeld and Maciej Bzowski, uses Capra to test morphology, normalization, grid transformations, resolution scaling, and performance. The present JOSS paper describes the software itself.

# State of the field

The standard IBEX map-making workflow and the THESEUS reconstruction provide important reference points for IBEX ENA imaging [@mccomas2009ibex; @funsten2009ibexhi; @osthus2024theseus]. Standard IBEX products are commonly represented on rectangular sky grids and include smoothing choices tied to the production of visually continuous maps. THESEUS uses a statistical reconstruction approach to infer sharper and smoother maps from the underlying observations. These methods are valuable, but they are not optimized for the specific role Capra fills: a transparent, HEALPix-native, modular reconstruction pipeline that keeps the bookkeeping of counts, exposure, signal, rate, and intensity explicit through reconstruction, normalization, transformation, and diagnostic stages.

Capra was built rather than added to an existing package because the research goal required a pipeline that combines mission-specific orbit-arc geometry with equal-area spherical pixelization, reconstruction diagnostics, and direct comparison products. General spherical-pixelization packages such as HEALPix and healpy provide the underlying sky discretization and coordinate tools [@gorski2005healpix; @zonca2019healpy], but they do not define the IBEX/IMAP-specific event-counting reconstruction, exposure handling, background subtraction, ENA intensity conversion, or GDF-normalization workflow. Capra builds on this ecosystem by supplying the mission-analysis layer needed to turn binned spacecraft observations into physically interpretable ENA maps.

# Software design

Capra is organized around the physical quantities that must be preserved during map reconstruction. For each orbit arc and scan-angle bin, the software reads the measured look direction, counts, exposure time, and background rate. It then identifies the relevant HEALPix pixels on the sphere, assigns or redistributes the bin contribution, and accumulates partial maps into all-sky products. Signal is computed from the difference between measured counts and exposure-scaled background. Signal rate and intensity are then derived from the accumulated signal and exposure, using the appropriate geometric and energy-channel factors.

A central design choice is the separation between reconstruction and interpretation. The baseline boresight-assigned product is the most direct representation of the binned measurement geometry. Optional smoothing is implemented as an explicit kernel-weighting step, with normalized discrete weights that preserve the bin contribution while redistributing it across nearby HEALPix pixels. This makes the smoothing choice inspectable and replaceable, rather than hidden inside plotting or downstream postprocessing.

A second design choice is to treat extensive and intensive quantities differently during grid transformations. Counts, exposure, and signal are summed when maps are transformed or aggregated, while rate and intensity are averaged as density-like quantities. This distinction is essential for preserving the physical meaning of exposure-corrected quantities. Capra therefore provides separate postprocessing routines for relative normalization, GDF normalization, rectangular-grid comparison products, coordinate-frame transformations, tessellation changes, and diagnostic checks.

The implementation uses Wolfram Language/Mathematica notebooks and packages for the core reconstruction and postprocessing workflow, with Python plotting utilities for HEALPix and rectangular-grid visualization. The reconstruction is parallelized over independent map-generation tasks, allowing high-resolution products to be generated efficiently while keeping the algorithmic structure close to the mathematical description.

![Capra workflow. The pipeline takes orbit-arc event-counting data, reconstructs HEALPix-native count, exposure, signal, rate, and intensity maps, and provides normalization, coordinate-transformation, tessellation-change, gridding, diagnostic, and comparison products.](capra_workflow.png)

# Research impact statement

Capra has been used in a companion research study of IBEX-Hi ENA map reconstruction, co-authored by Nikola Bukowiecka, Daniel B. Reisenfeld, and Maciej Bzowski. In that application, Capra was run on IBEX-Hi orbit-arc data from the second half of 2018, producing ENA intensity maps across multiple energy channels and enabling direct comparison with elementary IBEX maps and THESEUS reconstructions.

The study shows that Capra recovers the main large-scale ENA structures, including the Ribbon and broad nose-centered enhancement, while reducing sampling-driven mottling in low-exposure regions.

The same study uses Capra’s diagnostic tools to test preservation of counts, signal, rate, and intensity-like quantities under reconstruction and transformation. It also demonstrates HEALPix-native coordinate transformations, tessellation changes, and projection onto rectangular comparison grids. Resolution-scaling tests show convergence for the representative dataset by `Nside = 64`, and performance tests show approximately linear strong scaling over the tested range of parallel workers.

These results make Capra useful both as a research tool for current IBEX analyses and as a prototype workflow for future high-resolution IMAP ENA map-making. Its main near-term impact is to provide a reproducible, open-source, methodologically transparent reconstruction pipeline that can be inspected, modified, and compared against existing ENA map-making approaches.

The companion manuscript reports preservation tests, resolution-scaling tests, and parallel performance benchmarks; these results are summarized here only to document demonstrated software use, not to replace the full science paper.

# AI usage disclosure

Generative AI tools were used to assist with editing documentation. The scientific content, algorithmic choices, equations, implementation and validation tests were performed by the authors. AI-generated text was treated as editorial assistance and was not accepted without author review.

# Acknowledgements

The author thanks Daniel B. Reisenfeld and Maciej Bzowski for their collaboration on the companion scientific manuscript applying Capra to IBEX-Hi ENA map reconstruction, and for scientific discussions that helped shape the validation and comparison workflow used to demonstrate the software. The author also thanks the IBEX and IMAP communities for the mission context and data-analysis foundations that motivated this software.

The author acknowledges the hospitality of the Center for Computational Astrophysics-Flatiron Institute. The author also acknowledges support from the URI Institute for AI & Computational Research. The computations were performed on the UMass-URI UNITY high-performance computing (HPC) cluster hosted at the Massachusetts Green HPC Center.

# References
