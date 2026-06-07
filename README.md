# Capra

Capra is a HEALPix-native reconstruction pipeline for producing all-sky intensity maps from event-counting spacecraft data. It was developed for IBEX-Hi orbit-arc data and designed with future IMAP applications in mind.

The pipeline reconstructs count, exposure, background-subtracted signal, rate, and
intensity maps directly on the sphere using the HEALPix tessellation. It provides
both a baseline boresight-assigned product and an optional smoothed product, where
counts, exposure, and signal are redistributed over neighboring HEALPix pixels using
a user-selected kernel.

Capra is intended for transparent, reproducible ENA map-making and comparison
between reconstruction products such as elementary IBEX maps and THESEUS maps.

DOI 10.5281/zenodo.20584912

## Main features

- HEALPix-native all-sky reconstruction
- Baseline boresight-assigned maps
- Optional Gaussian smoothing kernel
- Count, exposure, signal, rate, and intensity products
- Detector-efficiency and calibration-factor support
- Relative and GDF-normalized intensity maps
- Transformation between HEALPix and rectangular comparison grids
- Coordinate-frame transformations
- HEALPix tessellation changes
- Validation diagnostics for signal, rate, flux, and map consistency
- Parallelized map generation

## Basic workflow

### 1. Prepare IBEX-Hi orbit-arc input data

The pipeline expects binned orbit-arc data containing, for each orbit arc and
scan-angle bin:

- ecliptic longitude
- ecliptic latitude
- counts
- exposure time
- background rate

The same orbit-arc data can be used to construct elementary IBEX comparison maps
and Capra HEALPix maps.

### 2. Choose the reconstruction settings (Run the code.nb)

Typical settings include:

- HEALPix resolution, e.g. `Nside = 16`, `64`, or `128`
- ESA energy step
- smoothing option:
  - no smoothing / boresight assignment
  - Gaussian smoothing
- collimator or kernel radius
- output directory

The baseline reconstruction assigns each scan-angle bin to the corresponding
HEALPix direction. The Gaussian-smoothed reconstruction redistributes each bin over
nearby HEALPix pixels using normalized discrete kernel weights.

### 3. Run the Capra reconstruction (Run the code.nb)

For each ESA step, Capra constructs native HEALPix maps of:

- counts
- exposure
- background-subtracted signal
- signal rate
- ENA intensity / flux

The basic intensity relation is

```math
F_j =
\frac{S_j}{E_j \, G_i \, C_{e,i}},
````

where `S_j` is the background-subtracted signal, `E_j` is exposure, `G_i` is the
geometric factor, and `C_{e,i}` is the central energy for ESA step `i`.

For the 2018 maps, the reconstructed signal rate or derived intensity maps may also
be divided by the detector-efficiency scale factor:

```math
\epsilon_\mathrm{eff}=0.96.
```

### 4. Postprocess the maps (Postprocess the code.nb)

Capra includes postprocessing routines for producing comparison products.

#### Relative normalization

Relative maps divide each intensity map by its solid-angle-weighted mean:

```math
F_\mathrm{rel} =
\frac{F}{\langle F \rangle_\Omega}.
```

These maps are useful for comparing morphology independent of absolute intensity
scale.

#### GDF normalization

GDF-normalized maps divide each intensity map by an off-Ribbon reference baseline,
estimated using an exposure-weighted median over a quiet off-Ribbon mask:

```math
F_\mathrm{GDF} =
\frac{F}{s_\mathrm{ESA}}.
```

These maps are useful for comparing Ribbon-to-background contrast.

### 5. Transform maps for comparison or visualization

Native HEALPix maps can be transformed to rectangular grids for comparison with
IBEX and THESEUS products. The transformation treats quantities according to their
physical type:

* intensity and rate: averaged as density-like quantities
* counts, exposure, and signal: summed as extensive quantities

The transformed maps can then be plotted with standard Python tools.

### 6. Run validation diagnostics

The diagnostic routines check consistency relations such as:

```math
R_j = \frac{S_j}{E_j},
```

and

```math
F_j =
\frac{S_j}{E_j G_i C_{e,i}}.
```

They also compare global signal, reconstructed signal, exposure-weighted rates,
and flux consistency factors. These checks are intended to verify that map
construction, smoothing, and grid transformations preserve the expected physical
relationships.

## Typical output products

For each ESA step, the pipeline can produce:

* native HEALPix count maps
* native HEALPix exposure maps
* native HEALPix signal maps
* native HEALPix rate maps
* native HEALPix intensity maps
* relative-normalized maps
* GDF-normalized maps
* rectangular-grid comparison maps
* diagnostic tables and plots

## Intended use

Capra is designed for ENA map reconstruction and method comparison. In its current
form it is primarily tested on IBEX-Hi data, but the HEALPix-native structure is
intended to make the approach adaptable to future IMAP ENA map-making workflows.

## Notes

The elementary IBEX maps used for comparison in this repository are not official
ISOC-rendered maps. They are direct maps made from the same orbit-arc data used as
input to the Capra pipeline and do not include the standard ISOC smoothing.


## Contributors
Nikola Bukowiecka led the analysis, software development and writing of the manuscript (in preparation).
Daniel Reisenfeld from Los Alamos National Laboratory acquired the observations and contributed relevant scientific expertise to the project and wiritng of the manuscript (in preparation).
Maciej Bzowski from Polish Space Research Centre, Polish Academy of Sciences contributed with geometricPackage.wl, and with early science advice to the project. This study was supported by the Polish Ministry for Education and Science under contract MEiN/2021/2/DIR. 

