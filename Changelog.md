# Changelog

## Unreleased

### Changed

- Rename the following constraints (M. Kirk)
  - `B->pi::FormFactors[parametric,LCSR+LQCD]@LvDM:2021A` -> `LMvD:2021A`
  - `B_s->D_s^*::A_1[s_max]@HPQCD:2019A` -> `HPQCD:2019B`
  - `B->K^*::A_12@FNAL+MILC:2013A` -> `HPQCD:2013B`
  - `B->K^*::A_1@FNAL+MILC:2013A` -> `HPQCD:2013B`
  - `B->K^*::V@FNAL+MILC:2013A` -> `HPQCD:2013B`
  - `B_s->K^*::A_12@FNAL+MILC:2013A` -> `HPQCD:2013B`
  - `B_s->K^*::A_1@FNAL+MILC:2013A` -> `HPQCD:2013B`
  - `B_s->K^*::V@FNAL+MILC:2013A` -> `HPQCD:2013B`
  - `D^+->Kbar^0e^+nu::BR[0.0,q2max]@PDG:2023A` -> `PDG:2024A`
  - `0->pipi::Abs{f_+}^2@CMD3:2023C` -> `CMD3:2023A`
  - `B->K^*::FormFactors[A_0,A_1,A_2,V,T_1,T_2,T_23]@GRvDV:2021A` -> `GRvDV:2022A`


### Added

- Add functionality to figure framework to allow export to multiple files at once (M. Kirk)
- Export `wilson-polynomials` to python (M. Reboud)
- Add missing references to bibliography (M. Kirk, E. McPartland)


### Removed

- Purge `integrate1d` (M. Reboud, F. Herren)


### Fixed

- Correct references list (M. Kirk, E. McPartland)



## [v1.0.19] - 2025-12-02

### Changed
### Added

- Start building wheels for Python version 3.14 (M. Kirk)
- Allow for log-scale plot in the new figure framework (M. Kirk)
- Include the ``eos-data``, ``eos-figure``, and ``eos-list-observables`` scripts in the Python wheels (M. Kirk)
- Add code to draw posterior-predictive samples in the 1D and 2D KDE (D. van Dyk)
- Add a differential dBR/dkperp distribution for $\bar{B}_s^0\to D_s^+\ell^-\bar\nu$ decays (D. van Dyk)

### Removed

### Fixed

- Fix a bug that prohibits building with gcc version 15 or newer (M. E. Smith)
- Fix a bug that selects a non-existent font as default in the figure framework (M. Kirk)
- Fix a [bug](https://github.com/eos/eos/issues/1115) in the Python bindings that broke the documentation (D. van Dyk)
- Fix various issues with LaTeX representations of observables (M. Kirk, D. van Dyk)



## [v1.0.18] - 2025-10-10

### Changed

- Improve the functionality of ``eos-analysis report`` (D. van Dyk)
- Improve the functionality of corner figures (M. Kirk, D. van Dyk)

### Added

- Add parameters to estimate power correction in class-I B decays (S. Meiser)
- Add (/ correct) documentation of the WET Wilson coefficient bases (D. van Dyk)
- Add tau decay observables (M. Kirk, D. van Dyk)
- Add pseudo observables for vacuum to K pi form factors (M. Kirk)

### Removed

### Fixed

- Fix a bug in the logging callback when running within a Python environment (D. van Dyk)
- Fix a [bug](https://github.com/eos/eos/issues/1078) that prohibited building against Boost version 1.89 or later (D. van Dyk)
- Fix a [bug](https://github.com/eos/eos/issues/1090) in the validation of analysis files (M. Kirk)
- Fix a [bug](https://github.com/eos/eos/issues/1091) in the (now deprecated) ``corner-plot`` task (M. Kirk)
- Fix a [bug](https://github.com/eos/eos/issues/1108) in the argument parsing for ``eos-analysis draw-figure`` (M. Kirk)



[v1.0.19]: https://github.com/eos/eos/releases/tag/v1.0.19
[v1.0.18]: https://github.com/eos/eos/releases/tag/v1.0.18
