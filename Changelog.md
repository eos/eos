# Changelog

## Unreleased

### Changed

- Rename the following signal PDFs (D. van Dyk)
  - `B->gammalnu::d^2Gamma/dEgamma/dcos(theta_l)` -> `B_u->gammalnu::P(E_gamma,cos(theta_l))`
  - `B_u->enumumu::d^5Gamma` -> `B_u->enumumu::P(q2,k2,z_gamma,z_w,phi)`
  - `B_u->munuee::d^5Gamma` -> `B_u->munuee::P(q2,k2,z_gamma,z_w,phi)`
  - `B->pilnu::dGamma/dq2` -> `B->pilnu::P(q2)`
  - `B->pilnu::d^2Gamma/dq2/dcos(theta_l)` -> `B->pilnu::P(q2,cos(theta_l))`
  - `B->Dlnu::dGamma/dq2` -> `B->Dlnu::P(q2)`
  - `B->Dlnu::d^2Gamma/dq2/dcos(theta_l)` -> `B->Dlnu::P(q2,cos(theta_l))`
  - `B->D^*lnu::dBR` -> `B->D^*lnu::P(q2)`
  - `B->D^*lnu::d^4Gamma` -> `B->D^*lnu::P(q2,cos(theta_l),cos(theta_D),phi)`
  - `B_s->K^*lnu::dBR` -> `B_s->K^*lnu::P(q2)`
  - `B_s->K^*lnu::d^4Gamma` -> `B_s->K^*lnu::P(q2,cos(theta_l),cos(theta_K),phi)`
  - `B->pipimunu::d^3Gamma@QCDF` -> `B^+->pi^+pi^-lnu::PDF(q2,k2,cos(theta_pi))`
  - `B->Knunu::dGamma/dq2` -> `B->Knunu::P(q2)`
  - `B->Kll::d^2Gamma@LargeRecoil` -> `B->Kll::P(q2,cos(theta_l))`
  - `B->K^*nunu::dGamma/dq2` -> `B->K^*nunu::P(q2)`
  - `B->K^*ll::d^2Gamma@LargeRecoil` -> `B->K^*ll::P(q2,cos(theta_l))`
- Make ``pypmc`` an optional Python dependency and adjust the documentation (D. van Dyk)
- Split the ``eos.data`` package into multiple modules, one per class (D. van Dyk)
- Allow the characters `|`, `(`, `)`, and `*` in ``qnp::OptionValue`` (D. van Dyk)
- Switch ``Options`` to use ``qnp::OptionValue`` in lieu of ``std::string`` for its option values; hard-coded option values are now expressed using the ``_ov`` user-defined literal (D. van Dyk)
- Use ``qnp::OptionValue`` keys in the the factory method ``Model::make`` (D. van Dyk)

### Added

- Add G2026 form factor parameterisation (N. Gubernari)
- Add B(s) -> K(*) constraints from BGMT:2025A (N. Gubernari)
- Add BSZ:2025A constraints in the Traditional Basis (N. Gubernari)
- Add B->pi f0 and f+ constraints from FLAG:2024A (N. Gubernari)
- Add Lambda_c->pll observables (D. Suelmann)
- Add the following signal PDFs (D. van Dyk)
  - `B_s->Klnu::P(q2)`
  - `B_s->Klnu::P(q2,cos(theta_l))`
  - `B_s->D_slnu::P(q2)`
  - `B_s->D_slnu::P(q2,cos(theta_l))`
  - `B->etalnu::P(q2)`
  - `B->etalnu::P(q2,cos(theta_l))`
  - `B->eta'lnu::P(q2)`
  - `B->eta'lnu::P(q2,cos(theta_l))`
  - `B_s->D_s^*lnu::P(q2)`
  - `B_s->D_s^*lnu::P(q2,cos(theta_l),cos(theta_D_s),phi)`
  - `B^+->pi^+pi^-lnu::PDF(q2,k2,cos(theta_pi))`
- Expand the unit test coverage of the Python interface (D. van Dyk)
- Add a workflow to determine code coverage with test cases (D. van Dyk)
- Add comparison operators for ``qnp::OptionValue`` (D. van Dyk)
- Add an output stream operator for ``qnp::OptionValue`` (D. van Dyk)
- Document the user-facing classes and methods of the ``eos`` Python modules, and auto-generate the analysis file format reference from the ``eos.analysis_file_description`` classes (D. van Dyk)
- Add a ``vertical`` plot item that draws a vertical line at a fixed position on the x axis, e.g. to mark a kinematic threshold (D. van Dyk)
- Add an optional ``size`` field to ``grid`` figures to set the figure size explicitly (D. van Dyk)
- Add an optional ``watermark_plot`` field to ``grid`` figures to stamp the watermark on a single panel, addressed either by a flattened index or by a ``(row, col)`` pair (D. van Dyk)
- Allow an explicit anchor (``xy``) and a configurable ``offset`` for the figure ``watermark`` (D. van Dyk)
- Add ``tight_layout`` and ``shared_axes`` options to ``grid`` figures (D. van Dyk)
- Provide a sensible default ``legend()`` entry for most figure item types (D. van Dyk)
- Make the tick label format configurable via a ``format`` field on the ``xticks`` and ``yticks`` of a plot (D. van Dyk)
- Add the unbinned log-likelihood block ``LogLikelihoodBlock::Unbinned1D`` (Python: ``eos.LogLikelihoodBlock.Unbinned1D``), which convolves a signal PDF with a resolution function on a grid via a discrete Fourier transform (D. van Dyk)
- Add classes for discrete Fourier transforms, ``eos::dft::Container`` and ``eos::dft::Plan`` (D. van Dyk)
- Add ``SignalPDF::evaluate_linear()`` to evaluate a signal PDF on the linear scale (D. van Dyk)
- Allow registering a custom ``SignalPDF`` at run time via ``SignalPDFs::insert`` (Python: ``eos.SignalPDFs.insert``) (D. van Dyk)
- Add an ``expression`` plot item that draws an arbitrary mathematical expression of the x-axis variable (D. van Dyk)
- Add a ``constraint2D`` plot item that draws the 2D contours of two correlated observables from a single constraint (D. van Dyk)
- Add a ``point`` plot item that draws a single marker at a fixed position (D. van Dyk)
- Add a ``contours2D`` plot item that draws 2D probability contours from a histogram of pre-existing samples (D. van Dyk)
- Add a ``scaling_factor`` field to the ``xticks`` and ``yticks`` of a plot to rescale the displayed tick values (D. van Dyk)

### Deprecated

- Deprecate the curtailed Gaussian prior description, i.e. a ``gauss``/``gaussian`` prior with ``min``/``max`` keys; use ``type: gaussian`` without ``min``/``max`` instead (D. van Dyk)
- Deprecate constraint-based prior descriptions without a ``type`` key; use ``type: constraint`` instead (D. van Dyk)

### Removed

- Remove the following univariate signal PDFs, since the corresponding multivariate PDFs render them superfluous (D. van Dyk)
  - `Lambda_b->Lambda_clnu::dBR`
  - `Lambda_b->Lambda_c(2625)lnu::dGamma`

- Remove the following signal PDFs, since they are unused and not worth maintaining (D. van Dyk)
  - `B->pimu1nu::d^2Gamma`
  - `B->pimu3nu::d^5Gamma`
  - `B->Dmu1nu::d^2Gamma`
  - `B->Dmu3nu::d^5Gamma`
  - `B->Kll::d^2Gamma@LowRecoil`
  - `B->K^*ll::d^2Gamma@LowRecoil`

### Fixed

- Fix some linewidth and linestyles arguments in figure framework (D. Suelmann)
- Fix undefined name errors in Python code (D. van Dyk)
- Fix usage of mutable default arguments in Python code (D. van Dyk)
- Fix exception handling in Python code (D. van Dyk)
- Remove dead code and unused imports in Python code (D. van Dyk)
- Fix exception handling in ``li22logA1`` (D. van Dyk)
- Fix use of ``std::abs`` in inclusive B->X_sll decays (D. van Dyk)
- Remove unused debugging code (D. van Dyk)
- Remove call to ``exit()`` in ``cubature`` code (D. van Dyk)
- Use ``context.nullcontext`` instead of ``context.suppress`` (D. van Dyk)
- Fix inconsistent use of ``s`` rather than ``q2`` across B->Plnu, Lambda_b->Lambda_clnu, Lambda_b->Lambda_c(2595)lnu, and Lambda_b->Lambda_c(2625)lnu (D. van Dyk)
- Fix inconsistent use of ``s`` rather than ``q2`` across D->Plnu, and Lambda_c->(S=1/2)^+lnu (D. van Dyk)
- Fix inconsistent use of ``s`` rather than ``q2`` across B->Kll, B->K^*ll, B_s->phill, Lambda_b->Lambdall, Lambda_b->Lambda(1520)ll, and B->X_sll (D. van Dyk)
- Fix bugs in figure framework (D. Suelmann)
- Fix bug in (non-)packaging of constraint files (D. van Dyk)
- Fix the default ``Item.legend`` to return a list rather than a tuple (D. van Dyk)
- Fix the type annotation of ``MaskComponent.description`` (D. van Dyk)
- Do not draw an empty legend (D. van Dyk)
- Fix a potential out-of-range bug when computing the contour level for ``eos.figure.TwoDimensionalKernelDensityEstimateItem`` (D. van Dyk)


## [v1.0.20] - 2026-04-28

### Changed

- Rename the following constraints (M. Kirk, M. Reboud)
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
  - `B_(s)->D_(s)^*::FormFactors@HPQCD:2023A` -> `B_(s)->D_(s)^*::FormFactors[HQETBasis]@HPQCD:2023A`
  - `B_(s)->D_(s)^*::FormFactors[TraditionalBasis]@HPQCD:2023A` -> `B_(s)->D_(s)^*::FormFactors@HPQCD:2023A`
- Impose the exact relation `f+(0) = f0(0)` in the BGL parameterisation of P -> P form factors (M. Reboud)
- Bump the Ubuntu image to resolute and the PyPI container (D. van Dyk)
- Reorganize the structure of constraints (F. Herren)


### Added

- Add functionality to figure framework to allow export to multiple files at once (M. Kirk)
- Export `wilson-polynomials` to python (M. Reboud, D. van Dyk)
- Add missing references to bibliography and implement a commit hook to check for missing references (M. Kirk, E. McPartland)
- Implement a new figure item to plot constraints residues (M. Reboud)
- Update the examples and the tests of the figure framework (M. Reboud)
- Implements an observable option to specify which option needs to be plotted in an uncertainty band item (M. Reboud)
- Implements B_c -> J/psi ell nu decays (C. Bolognani, M. Reboud)
- Add basis for WCs with Delta c = Delta u = 1 transitions (D. Sue)
- Extend the WC RGE and test it (D. Sue)
- Implement L_c -> neutron ell nu decays (C. Bolognani)


### Removed

- Purge `integrate1d` (M. Reboud, F. Herren)
- Remove the matplotlib backend choice (M. Kirk)


### Fixed

- Correct references list (M. Kirk, E. McPartland)
- Fix bugs in the figure framework (M. Reboud)
- Fix the constraint `B_(s)->D_(s)^*::FormFactors@HPQCD:2023A` (M. Reboud)
- Fix the bound saturations in the BGL parameterisation (M. Reboud)
- Fix the font used by matplotlib (M. Kirk)
- Fix a race condition in TreadPool (D. van Dyk)


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



[v1.0.20]: https://github.com/eos/eos/releases/tag/v1.0.20
[v1.0.19]: https://github.com/eos/eos/releases/tag/v1.0.19
[v1.0.18]: https://github.com/eos/eos/releases/tag/v1.0.18
