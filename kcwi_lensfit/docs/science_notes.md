# Science notes: why every default is what it is

The science goal is **redshifts** (and, where the data allow, **stellar kinematics**) of lens
systems, with emphasis on **compound lenses** (several galaxies/source planes along one sightline).
This file justifies the physical choices; every one is a setting you can change.

## 1. Data handling

**Wavelength frame.** Air and vacuum wavelengths differ by 2.8×10⁻⁴ (≈85 km/s, Δz≈3×10⁻⁴ at any z).
Keck DRP blue cubes use `CTYPE3=AWAV`, red cubes `WAVE` (vacuum). The loader reads `CTYPE3`, converts
everything to **air** (Ciddor/Morton, the pPXF convention), and the fitting code converts the
*wavelength array* (not a constant factor) into the frame of the chosen template library: E-MILES is
air-native, FSPS/GALAXEV/XSL vacuum-native. Line lists are stored in vacuum and converted at use.

**Velocity frame.** DRP cubes already carry `VCORR` (heliocentric). The default does *not* re-apply
it (double correction shifts z by ~1×10⁻⁴); the log prints the keyword so you can decide.

**Footprint.** DRP stacks store “no data” as exact zeros (or NaN). Both become NaN so a zero-valued
border can never be selected as “quiet sky” (which would give zero noise and absurd S/N).

**Units.** `BUNIT` per arcsec² is multiplied by the spaxel area so sums over spaxels are true fluxes.

**Registration.** Blue and red are two cameras behind one slicer: the same (x, y) is the same sky by
design, so `identity` is a valid choice, but the default is `wcs` (follow the header WCS; verify it in the
stage-2 overlay and switch to `identity` or `xcorr` if the zero-points disagree). Header WCS zero-points
are set independently and can be arcseconds off; image cross-correlation is a weak test because a blue star-forming arc and a red
quiescent core do not trace the same structure. Both are *computed and displayed* for comparison.
Sky-coordinate region files are, at the user's request, mapped with each cube's own WCS.

## 2. Noise

Measured from a source-free box: background = median, noise = 1.4826·MAD per wavelength. Resampling
the image-slicer data onto a square grid **correlates neighbouring spaxels**, so the sum over *N*
spaxels has σ_N = σ₁·N^α with α between 0.5 (independent) and 1 (fully correlated). α is fitted
from block sums of 1×1…k×k spaxels in the sky box and extrapolated to the aperture size. Ignoring it
makes errors of summed apertures too small, often by 1.5–3×. The error of the background estimate
itself (perfectly correlated across the aperture) is added. Source photon noise is not modelled;
pPXF errors are multiplied by √χ²_red when χ²_red > 1 and are bootstrapped.

## 3. Extraction and deblending

At KCWI seeing (~1″) the deflector halo, the arcs and (in compound lenses) a second lens overlap.
Subtracting a scaled copy of one aperture's spectrum from another needs the contaminant's
*spatial* profile; a single circular Gaussian (the quick approach) fails for extended lens light.
The default instead fits, once, a spatial model of every object in a white-light image of the band
where the deflector dominates (**PSF-convolved elliptical Sérsic** for galaxies, the **seeing disc**
for point-like images, the **aperture mask itself** for arcs whose true shape is unknown), then solves

  data(x, λ) = Σₖ cₖ(λ)·Pₖ(x)

by linear least squares at every λ. The covariance of the fit is propagated: two blended objects
cannot both be measured as precisely as either alone, and the error bars say so. With one object this
reduces to optimal (profile-weighted) extraction. Diagnostics: χ² of the profile fit, template
condition number, template-degeneracy warnings and the residual image. Plain sums (optionally
intensity-weighted, i.e. the co-add-then-fit aperture of Cappellari+2006) remain available.

**Images of one source** share a `group` label and are inverse-variance combined: lensing is
achromatic, so every image has the same redshift (time-delay variability changes the spectral
shape, not z). Different source planes of a compound lens get different groups.

## 4. Redshifts of sources: matched-filter line search

Continuum is removed with a running median (150 Å) and the noise is rescaled to the local empirical
scatter (a bad variance cannot create detections). At each trial z every catalogue **multiplet** (e.g.
[O III] 4960/5008 with the exact 1:2.98 ratio; Mg II, C IV, Si IV doublets; [O II]) has a best-fit
amplitude A and error σ_A; the signed S/N s = A/σ_A must be positive for emission and negative for
absorption. The score

  S(z) = Σ_covered lines [ max(s_signed, 0)² − ½ ]

subtracts the null expectation (½) *per covered line*, so redshifts where more lines happen to fall
in the wavelength coverage are not favoured over redshifts with fewer. (The naive “sum of S/N” score
systematically prefers low-z aliases with many lines in range over a real single-line high-z
detection, and mis-ranks the true answer.) The search needs no initial guess and covers emission and
absorption, because lensed Lyman-break-type sources are often identified only through interstellar
absorption (Si II, C II, Fe II, Al II, Mg II …).

**Line classes and systemic redshift.** Nebular emission and stellar absorption trace the systemic
redshift; interstellar/outflow absorption is typically blueshifted by 50–300 km/s; resonant
emission (Lyα, Mg II, C IV) is typically redshifted by 100–500 km/s. The joint fit therefore
gives each class its own velocity offset relative to the reference class. If only outflow-affected
lines are present the redshift is *defined by that class*, a **systemic estimate** and a
**systematic error** (≈150 km/s for ISM absorption, ≈250 km/s for resonant emission) are attached.

**A single line is never “secure”.** `min_lines_secure` (default 2) multiplets must be detected.
Otherwise the result is *tentative* and the report tabulates every possible identification of the
feature ([O II], Mg II, Hβ, [O III], Hα, …), the redshift each implies, and how many *other* lines
should then be visible versus actually detected.

## 5. Redshifts and σ of deflectors: pPXF

* **Scan**: χ²(z) with σ fixed (a free σ lets a noisy fit “cheat” by smearing the template over
  noise, flattening χ²(z)); each trial is a non-negative linear fit of high-pass-filtered templates
  (continuum shape, reddening and any blue/red scale mismatch drop out), so the whole range is
  scanned in seconds; the best minima are refined with full pPXF.
* **Joint blue+red** spectrum, log-rebinned **once** (separate rebinning breaks the constant
  `velscale` requirement), one z and σ, a free red/blue scale, **multiplicative** polynomials by
  default (they cannot erase or fake line depth, unlike additive ones).
* **Resolution matching.** `sps_lib` can only broaden templates. For z ≳ 0.5 the rest-frame KCWI
  resolution (~1.5 Å) is *sharper* than E-MILES (2.5 Å); the galaxy is then degraded to the template
  floor so σ is not inflated (a synthetic test without this recovered 389 km/s for a true
  235 km/s). The noise is deliberately **not** reduced after smoothing (correlated pixels carry no
  extra information).
* **Rest window.** Default 3540–8950 Å: where E-MILES has its native 2.5 Å resolution. In the
  rest-UV the library is ~5 Å, so including it makes σ depend on the template. Redshifts do not care.
* **Source emission in the deflector window.** A z=1.3 source's Mg II/[O II]/He I/[Ne III] fall in
  the deflector's absorption region and bias σ if left in. Once any plane's z is known its emission
  lines are masked in the other planes' *continuum* fits (`mask_other_planes_lines`). Redshift searches are never masked, and a region found at the same redshift as another (within `same_plane_tol_kms`, default 1500 km/s) keeps that plane's lines, so ungrouped regions that may be the same source are not blinded; the report flags such pairs as “SAME PLANE?”. Balmer lines are *kept* for
  passive galaxies (they are stellar absorption that constrains σ); only forbidden lines are masked.
* **Errors.** Bootstrap (default: noise Monte Carlo around the best-fit model; optional wild/Rademacher resampling of the residuals with *locally* scaled scatter);
  **systematics** from a sweep over multiplicative-polynomial degrees and, optionally, template
  libraries (each in its own native frame). Polynomial/library choices are usually the dominant σ
  systematic, which is why they are reported, not hidden.
* **Reliability.** σ is flagged unreliable if S/N per Å < 10, if σ is below half a pixel, or if
  only empirical (5 Å) templates were used. The redshift is still reported.
* **Cross-checks.** Blue-only vs red-only (independent detectors/calibrations) must agree on z, and
  an independent stellar-absorption-line fit (Ca H&K, G band, Mg b, Na D, Ca triplet) is compared.

On σ_e: quote an *aperture* σ with its radius (e.g. R_eff) and apply your own aperture-correction
convention; use `weighting: intensity` with a circular aperture for the I-weighted aperture used in
the literature.

## 6. Lens-system consistency (reported, never imposed)

* Every source must lie **behind the primary (lowest-z) deflector**; otherwise it is an interloper or
  the line identification is wrong. Candidates at or below the deflector's redshift are marked
  impossible in the candidate list and skipped by the automatic choice.
* Images in one group must agree within their errors; if not they are different sources or a
  spectrum is contaminated.
* A source within 2000 km/s of the deflector is probably a companion, not a background galaxy.
* The report prints the **redshift ladder**; ≥3 distinct planes with deflector + sources identifies
  a compound (multi-plane) lens.

## 7. Known simplifications

Single Gaussian seeing independent of wavelength; instrumental FWHM constant per arm (nominal
R(grating, slicer) unless overridden — a measured value from arcs/sky lines is better and sets the σ
zero point); no telluric correction; rest-frame equivalent widths and fluxes are not calibrated.
