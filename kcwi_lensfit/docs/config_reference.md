# Configuration reference

Generated from `kcwilensfit/schema.py` (`kcwilensfit fields` prints the same text). Settings marked *advanced* are hidden behind “Show advanced settings” in the GUI; their defaults suit most data.


## 1. Data

### `target_name` — Target name

*default:* `'lens_system'`

Label used in output file names and plot titles (e.g. 'AGEL012441').

### `blue_cube` — Blue (KCWI) intensity cube

*default:* `''`

Reduced KCWI data cube (FITS, 3-D: wavelength x y x x). Leave empty to use the red arm only. A cube from the Keck DRP '_icubes' (intensity) product is expected; NaNs or all-zero spaxels are treated as 'no data'.

### `red_cube` — Red (KCRM) intensity cube

*default:* `''`

Reduced KCRM data cube. Leave empty to use the blue arm only. Both arms are fitted jointly (one redshift, one velocity dispersion) when both are given.

### `blue_var_cube` — Blue variance cube (optional)  *(advanced)*

*default:* `''`

Variance cube ('_vcubes') matching the blue flux cube. If given it is used as the noise model; otherwise the noise is measured from a source-free sky box (Noise & sky tab).

### `red_var_cube` — Red variance cube (optional)  *(advanced)*

*default:* `''`

Variance cube matching the red flux cube (see above).

### `output_dir` — Output directory

*default:* `'kcwilensfit_output'`

Everything the pipeline writes (spectra, plots, tables, the settings actually used) goes here.


## 1b. Instrument (same tab as Data)

### `blue_wave_frame` — Blue wavelength frame  *(advanced)*

*default:* `'auto'` · *choices:* auto, air, vacuum

Whether the blue cube's wavelength axis is in AIR or VACUUM. 'auto' reads CTYPE3 (AWAV = air, WAVE = vacuum) and assumes air if absent. Getting this wrong shifts every line by ~85 km/s (dz ~ 3e-4) and goes straight into the redshift.

### `red_wave_frame` — Red wavelength frame  *(advanced)*

*default:* `'auto'` · *choices:* auto, air, vacuum

Same for the red cube. 'auto' reads CTYPE3 and assumes VACUUM if absent (the KCRM DRP convention).

### `velocity_correction` — Velocity-frame correction

*default:* `'none'` · *choices:* none, header, manual

'none': use wavelengths as they are in the file (the Keck DRP normally already writes heliocentric-corrected cubes; check the VCORR/VCORRTYP keywords printed in the log). 'header': ADD the header VCORR on top (only if the file is NOT already corrected). 'manual': add the value set under Advanced. Double-correcting shifts z by ~1e-4.

### `manual_vcorr_kms` — Manual velocity correction [km/s]  *(advanced)*

*default:* `0.0`

Used only if the correction mode is 'manual'. Positive values make redshifts larger.

### `blue_wave_min` — Blue min wavelength [A]  *(advanced)*

*default:* `0.0`

Trim the blue cube to this observed wavelength. 0 = automatic (header WAVGOOD0 if present, else the first valid channel). Edges of the bandpass are low-throughput and noisy.

### `blue_wave_max` — Blue max wavelength [A]  *(advanced)*

*default:* `0.0`

As above for the red end of the blue cube. 0 = automatic. The dichroic cross-over region near 5600 A is usually the least reliable part of the blue arm.

### `red_wave_min` — Red min wavelength [A]  *(advanced)*

*default:* `0.0`

Trim the red cube (0 = header WAVGOOD0 if present, else first valid channel).

### `red_wave_max` — Red max wavelength [A]  *(advanced)*

*default:* `0.0`

Trim the red cube (0 = header WAVGOOD1 if present, else last valid channel).

### `blue_grating` — Blue grating  *(advanced)*

*default:* `'auto'` · *choices:* auto, BL, BM, BH1, BH2, BH3

KCWI grating. 'auto' reads the header / file name; if neither identifies it, BL is assumed.

### `red_grating` — Red grating  *(advanced)*

*default:* `'auto'` · *choices:* auto, RL, RM1, RM2, RH1, RH2, RH3, RH4

KCRM grating ('auto' as above; RL assumed if undetected).

### `slicer` — Image slicer  *(advanced)*

*default:* `'auto'` · *choices:* auto, Large, Medium, Small

IFU slicer (Large/Medium/Small). Sets the resolving power R together with the grating. 'auto' reads IFUNAM; Medium assumed otherwise.

### `blue_fwhm_A` — Blue instrumental FWHM [A]

*default:* `0.0`

Instrumental line-spread FWHM of the blue arm at its central wavelength. LEAVE EMPTY to use the nominal value lambda_c/R from the grating+slicer table. Enter a measured value (arc lamp / sky lines) if you have one -- it directly sets the velocity-dispersion zero point (sigma_inst = FWHM/2.355).

### `red_fwhm_A` — Red instrumental FWHM [A]

*default:* `0.0`

As above for the red arm (empty = nominal table value).

### `convert_sb_to_flux` — Convert surface brightness to flux per spaxel  *(advanced)*

*default:* `True`

If the cube BUNIT is per arcsec^2 (the DRP 'icubes'), multiply by the spaxel area so extracted spectra are true fluxes (1e-16 erg/s/cm2/A) and summing spaxels is physically meaningful. Has no effect on redshifts.

### `arm_registration` — Blue/red spatial registration

*default:* `'wcs'` · *choices:* identity, wcs, xcorr, manual

How a spaxel in the blue cube maps to the red cube. 'identity': same (x,y) index -- correct when both arms were resampled to one grid from the same IFU exposure (the KCWI+KCRM design: one telescope focus, one image slicer, split by a dichroic). 'wcs': follow the header WCS. 'xcorr': measure from white-light images. 'manual': use dx, dy (under Advanced). Header WCS zero-points can be off by arcsec, so verify in the registration check (Data tab) before trusting 'wcs'.

### `arm_dx` — Manual dx [pix]  *(advanced)*

*default:* `0.0`

Red x = blue x + dx (only for 'manual').

### `arm_dy` — Manual dy [pix]  *(advanced)*

*default:* `0.0`

Red y = blue y + dy (only for 'manual').


## 2. Noise & sky

### `sky_box_mode` — Sky-region selection

*default:* `'auto'` · *choices:* auto, manual

'auto' picks the faintest, flattest source-free box inside the valid field of view. 'manual' uses the box given below. The box supplies (i) the residual background level to subtract and (ii) the empirical noise -- it must contain NO source or lens light.

### `sky_box_size` — Auto sky-box size [pix]

*default:* `12`

Side of the square box searched for in 'auto' mode. Bigger = better statistics but harder to keep source-free in a small field.

### `sky_box_blue` — Manual sky box, blue (y0,y1,x0,x1)

*default:* `''`

Pixel limits, 0-indexed, end-exclusive, as 'y0,y1,x0,x1'. Ignored in 'auto' mode.

### `sky_box_red` — Manual sky box, red (y0,y1,x0,x1)

*default:* `''`

Leave empty to reuse the blue box for the red arm (valid with identity registration).

### `subtract_background` — Subtract residual background

*default:* `True`

Subtract the sky-box median spectrum from every spaxel. Removes small sky-subtraction residuals common to the whole field. Turn OFF if the box might contain diffuse lens or cluster light (it would then over-subtract real signal).

### `noise_source` — Noise model  *(advanced)*

*default:* `'auto'` · *choices:* auto, sky_box, variance_cube

'sky_box': per-wavelength scatter of the sky box (assumes spatially uniform noise). 'variance_cube': use the supplied variance cubes. 'auto': variance cubes if provided, else sky box. Neither contains source photon noise, so bright apertures are slightly under-estimated -- pPXF errors are rescaled by sqrt(chi2/dof) to compensate.

### `correct_noise_correlation` — Correct for correlated spaxel noise  *(advanced)*

*default:* `True`

Resampling the KCWI slices onto a square grid correlates neighbouring spaxels, so the noise of a SUM over N spaxels grows faster than sqrt(N). The pipeline measures this from block sums in the sky box (sigma_N = sigma_1 N^alpha, alpha in [0.5,1]) and inflates aperture errors accordingly. Skipping it makes every error too optimistic.


## 3. Apertures

### `sky_regions_use_wcs` — Map sky-coordinate regions with each cube's WCS

*default:* `True`

Regions loaded from a DS9 file in fk5/icrs (WCS) coordinates are converted to pixels separately for the blue and the red cube using EACH cube's own header WCS (centre, size and rotation), and the registered shapes actually used are drawn on every plot. Turn OFF to treat the regions as defined once on the blue grid and shift them with the blue/red registration (Data tab) instead (safer if the header WCS zero-points disagree).


## 4. Extraction

### `extraction_method` — Extraction method

*default:* `'deblend'` · *choices:* deblend, sum

'deblend': simultaneous spatial-profile (PSF-convolved Sersic / point / arc-mask) fit of ALL apertures at every wavelength -- spectral deblending, the right tool when a deflector and its arcs (or two lens galaxies) overlap on the sky. Errors include the blending covariance. 'sum': plain (optionally flux-weighted) aperture sums with no contamination correction -- adequate only for well-separated, isolated objects.

### `weighting` — Weighting for 'sum' extraction  *(advanced)*

*default:* `'uniform'` · *choices:* uniform, intensity

'uniform': plain sum. 'intensity': weight spaxels by the white-light flux (the S/N-optimal, co-add-then-fit aperture of Cappellari+2006; suppresses faint contaminating arc light in a deflector aperture).

### `seeing_fwhm_arcsec` — Seeing / PSF FWHM [arcsec]

*default:* `1.0`

Gaussian PSF FWHM used to build the spatial profiles for deblending. KCWI is seeing-limited (typically 0.7-1.3"); take it from the header or a field star. The Sersic fit is degenerate with the PSF, so a wrong value mainly biases the profile r_e.

### `profile_band_A` — Band used to fit profiles [A]  *(advanced)*

*default:* `''`

Observed wavelength range 'lo-hi' of the image used to fit the spatial profiles. Empty = the redder half of the red arm (or blue arm if only blue is given), where the deflector dominates and the blue star-forming arcs are faintest.

### `deblend_pad_pix` — Deblend window padding [pix]  *(advanced)*

*default:* `8`

Margin added around the aperture bounding box for the spatial fit region. The profile wings outside the apertures are what pin down the underlying lens light.

### `deblend_fit_background` — Fit a flat background in deblend  *(advanced)*

*default:* `False`

Add a free constant to the deblend model at every wavelength. Usually unnecessary after the sky-box subtraction and it can be degenerate with extended lens light.

### `combine_groups` — Combine images of the same source

*default:* `True`

Apertures that share a 'group' label (multiple images of ONE lensed source) are inverse-variance combined into an extra spectrum '<group>' for a higher-S/N redshift. Valid because lensing is achromatic -- all images share one redshift -- but only label apertures with the same group when you are confident they are images of the same source.

### `manual_subtractions` — Manual subtraction of one region from another  *(advanced)*

*default:* `''`

Remove light you KNOW belongs to another region, one rule per line: 'target: contaminator [scale]'. Example: 'imgA: deflector 0.3' subtracts 0.3 x the deflector's spectrum from imgA's; 'imgA: deflector auto' (or no scale) estimates the fraction from the contaminator's fitted spatial profile at the target aperture. Works with both extraction methods (with 'deblend' the overlap is already modelled, so use it only for residuals). The subtracted spectrum's noise is added to the target's variance. Uses the contaminator's spectrum as extracted, i.e. before any other manual subtraction.

### `reject_spikes` — Reject single-pixel spikes  *(advanced)*

*default:* `True`

Mask isolated one-pixel positive outliers (cosmic-ray / hot-pixel residuals) in the extracted spectra. Real emission lines are resolved by the instrument (>= ~2 pixels wide), so a lone bright pixel cannot be a line.

### `display_smooth_A` — Display smoothing [A]  *(advanced)*

*default:* `0.0`

Gaussian smoothing (sigma, Angstrom) applied to plots ONLY. Fits always use unsmoothed data: smoothing correlates neighbouring pixels and fakes S/N.


## 5. Redshift & kinematics

### `fit_arms` — Arms to fit

*default:* `'both'` · *choices:* both, blue, red

'both': blue+red fitted jointly with a single redshift/sigma and one free relative flux scale for red (only the ratio between arms is identifiable). Otherwise a single arm.

### `fit_red_scale` — Fit relative red/blue flux scale

*default:* `True`

(Only when both arms are fitted.) Fit one free multiplicative factor for the red arm relative to the blue (iterated in closed form). Corrects flux-calibration differences between the arms.

### `crosscheck_arms` — Blue-only vs red-only consistency

*default:* `True`

(Only when both arms are fitted.) Re-fit each arm alone: independent determinations from the two cameras (different detectors, calibrations and systematics) must agree on z.

### `mask_airglow` — Mask night-sky emission lines  *(advanced)*

*default:* `True`

[O I] 5577/6300/6364 and Na D airglow: residuals here are atmospheric, not astrophysical at any redshift.

### `mask_oh_forest` — Mask the OH forest (> 8750 A)  *(advanced)*

*default:* `True`

The crowded red OH Meinel bands, at the very edge of the usable KCRM range.

### `mask_light_pollution` — Mask Hg/Na light-pollution lines  *(advanced)*

*default:* `True`

Street-lamp lines; weak on Maunakea but cheap insurance.

### `mask_telluric` — Mask telluric absorption bands  *(advanced)*

*default:* `True`

O2 B (6860-6890), H2O (7150-7350) and O2 A (7580-7700) A. No telluric correction is applied, so these bands are excluded from fits.

### `extra_mask_windows` — Extra masked windows [A]

*default:* `''`

Additional observed-frame windows to exclude, one 'lo-hi' per line (e.g. a cosmic-ray residual or a bad column).

### `mask_other_planes_lines` — Mask other planes' lines in each fit

*default:* `True`

Compound-lens essential: once a redshift is secured for ANY aperture, its emission lines are masked in the continuum fits of the other apertures (e.g. a z=1.3 source's [O II]/Mg II emission lands inside the deflector's absorption window and biases sigma if left in). Never applied to a redshift SEARCH (a region may itself lie at another plane's redshift), and never to a region that turns out to be at the same redshift (see the same-plane tolerance).

### `other_line_mask_kms` — Line-mask half-width [km/s]  *(advanced)*

*default:* `600.0`

Half-width of the mask around each masked emission line (covers line width + velocity uncertainty of the redshift).

### `same_plane_tol_kms` — Same-plane tolerance [km/s]  *(advanced)*

*default:* `1500.0`

When fitting a region at (or near) redshift z, other regions whose redshift lies within this velocity of z are treated as the SAME source plane and their emission lines are NOT masked (they would otherwise erase this region's own lines). Also marks ungrouped regions with agreeing redshifts as 'possibly one source' in the report.

### `sps_library` — Stellar template library

*default:* `'emiles'` · *choices:* emiles, fsps, galaxev, xsl, kc96_starburst, kc96_galaxy, kc96_all, custom

Template library for continuum fits. 'emiles' (R~2.5 A FWHM, air) is the best match to KCWI resolution for old stellar populations; 'fsps'/'galaxev'/'xsl' are alternatives; 'kc96_*' are empirical galaxy spectra (starburst/spiral/elliptical, 5 A) useful for blue star-forming sources; 'custom' uses your own pPXF-format .npz. Downloaded on first use.

### `source_sps_library` — Template library for sources

*default:* `'kc96_starburst'` · *choices:* emiles, fsps, galaxev, xsl, kc96_starburst, kc96_galaxy, kc96_all, custom

Library used for continuum fits of 'source'/'other' apertures. Lensed sources are usually blue star-forming galaxies whose spectra are better matched by empirical starburst templates than by old-population SPS models.

### `custom_template_path` — Custom template file (.npz)  *(advanced)*

*default:* `''`

pPXF sps_lib-format .npz (templates, lam, fwhm, ages, metals, masses). Used if the library is 'custom'.

### `custom_template_frame` — Custom template frame  *(advanced)*

*default:* `'air'` · *choices:* air, vacuum

Wavelength frame of the custom file.

### `age_range` — Template age range [Gyr]  *(advanced)*

*default:* `'0.5-14.0'`

Restrict SPS templates to ages 'lo-hi' Gyr (ignored for the empirical kc96 libraries). Rarely needs changing: the default drops very young stars, which are degenerate with the continuum polynomial for a passive deflector; widen it only for a star-forming deflector.

### `metal_range` — Template [M/H] range  *(advanced)*

*default:* `'-1.0-0.3'`

Restrict SPS templates to [M/H] 'lo-hi' (ignored for the empirical kc96 libraries). Rarely needs changing: massive lens galaxies are metal-rich and this removes only very sub-solar templates that add degeneracy.

### `moments` — Kinematic moments  *(advanced)*

*default:* `'2'` · *choices:* 2, 4

2 = (V, sigma). 4 adds Gauss-Hermite h3, h4 (needs S/N >~ 30 per pixel to be useful).

### `mdegree` — Multiplicative polynomial degree  *(advanced)*

*default:* `8`

Order of the multiplicative polynomial. It absorbs the blue/red flux-calibration mismatch and slow reddening/continuum mismatch WITHOUT distorting line depths (unlike additive polynomials), so it is the safer choice when sigma matters. Raise it for a very wide wavelength range.

### `degree` — Additive polynomial degree  *(advanced)*

*default:* `-1`

Order of the additive polynomial (-1 = none). Additive terms can fake or hide line depth; use sparingly and test the effect with the robustness sweep.

### `sigma_init` — Initial sigma [km/s]  *(advanced)*

*default:* `200.0`

Starting velocity dispersion for the non-linear fit.

### `vel_bound_kms` — Velocity bound [km/s]  *(advanced)*

*default:* `1500.0`

|V| is restricted to +/- this about the chosen redshift guess. A result sitting on the bound is flagged as unconverged.

### `mask_nebular` — Mask nebular emission in continuum fits  *(advanced)*

*default:* `'forbidden'` · *choices:* forbidden, all, none

Lines at the aperture's OWN redshift that stellar templates cannot reproduce. 'forbidden': [O II], [Ne III], [O III], [N II], [S II] (Balmer lines are kept because for a passive galaxy they are stellar ABSORPTION that constrains sigma). 'all': also Balmer. 'none'.

### `fit_rest_range` — Rest-frame fitting window for deflectors [A]

*default:* `'3540-8950'`

Rest-frame wavelength range used in stellar-continuum fits of 'deflector' apertures (empty = everything). The default is the range where E-MILES has its native ~2.5 A resolution; in the rest-UV the libraries are much coarser than KCWI, so including it forces the galaxy to be degraded and makes sigma depend on the template. Redshifts are insensitive to this choice; velocity dispersions are not.

### `deflector_z_range` — Deflector z search range

*default:* `'0.05-1.2'`

Redshift range 'lo-hi' scanned for aperture role 'deflector'. Strong-lens deflectors are massive early-type galaxies, mostly 0.1 < z < 1; Ca H&K must lie in the data (z < 1.27 with KCRM to 8940 A).

### `source_z_range` — Source z search range

*default:* `'0.2-3.5'`

Redshift range scanned for aperture role 'source'. Physically a lensed source must be BEHIND the lens, so results at or below a measured deflector redshift are flagged as impossible.

### `z_step_continuum` — z step, continuum scan  *(advanced)*

*default:* `0.002`

Spacing of the chi2(z) scan (velocity dispersion held fixed, continuum removed by high-pass filtering). Must be smaller than the width of the narrowest strong features (sigma ~ 250 km/s -> dz ~ 0.0015 at z~0.7); the best minima are then refined with full pPXF fits.

### `z_step_lines` — z step, line search  *(advanced)*

*default:* `0.0005`

Spacing of the matched-filter line scan. Lines are ~150-300 km/s wide, so steps of dz ~5e-4 (150 km/s) sample them adequately.

### `scan_sigma_kms` — sigma(s) held fixed in the scan [km/s]

*default:* `'150, 350'`

Velocity dispersion(s) assumed while scanning for the redshift. Give ONE value for a single scan, or several (default '150, 350') to run a scan for each: a candidate present at every value is stable against the assumption, one that appears only for some is flagged. Holding sigma fixed stops a noisy fit from smearing the template over the noise; a wrong value mainly reduces contrast and rarely shifts the minimum. Massive lens galaxies span roughly 150-450 km/s.

### `line_sets` — Line search uses  *(advanced)*

*default:* `'both'` · *choices:* both, emission, absorption

Which kinds of features the blind line search scores: 'emission' (nebular + resonant), 'absorption' (stellar + ISM) or 'both'. Source arcs of star-forming galaxies often show only one or the other.

### `line_sigma_kms` — Intrinsic line width [km/s]  *(advanced)*

*default:* `100.0`

Assumed intrinsic Gaussian sigma of searched lines, added in quadrature to the instrumental width. 100 km/s suits compact star-forming galaxies; use 250+ for AGN or massive galaxies.

### `line_continuum_window_A` — Continuum filter width [A]  *(advanced)*

*default:* `150.0`

Width of the running-median filter that removes the continuum before the matched-filter search. Must exceed ~3x the widest feature of interest.

### `detect_snr` — Per-line detection threshold [sigma]  *(advanced)*

*default:* `3.0`

A line counts as 'detected' when its fitted amplitude is at least this many sigma.

### `n_null_sims` — Noise simulations for significance  *(advanced)*

*default:* `20`

Number of sign-flipped noise realisations of the whole redshift scan used to measure how large the best score gets in pure noise (look-elsewhere correction: a scan of ~10^4 redshifts x ~45 lines produces chance 3-sigma peaks). A candidate must beat the 95th percentile of these maxima to be called significant. 0 disables the test.

### `min_lines_secure` — Lines needed for a secure z

*default:* `2`

Independent line features that must be detected before a redshift is called 'secure'. A single line is never proof of a redshift (an isolated emission feature can be [O II], Mg II, H-beta, [O III], H-alpha...); the report lists the competing identifications.

### `n_bootstrap` — Bootstrap realisations

*default:* `100`

Number of re-fits used to measure the redshift and sigma uncertainties (0 = use pPXF's formal error only, which is usually optimistic). 100 gives ~7% precision on the error; use 500+ for publication numbers.

### `bootstrap_method` — Bootstrap method  *(advanced)*

*default:* `'mc'` · *choices:* wild, mc

'mc' (default): re-draw Gaussian noise from the noise array around the best-fit model and re-fit -- the error then follows the measured noise level. 'wild': multiply the fit residuals by random +/-1 (Rademacher), scaled by their local scatter -- also captures template mismatch, but shrinks if the noise array is over-estimated.

### `sweep_degrees` — Polynomial-degree sweep  *(advanced)*

*default:* `'6,8,10,12'`

Comma list of multiplicative degrees re-fit after the main fit; the scatter in V and sigma is reported as a SYSTEMATIC error (empty = off). Continuum modelling is usually the largest systematic in a sigma measurement.

### `sweep_libraries` — Template-library sweep  *(advanced)*

*default:* `''`

Comma list of libraries (e.g. 'emiles,fsps') re-fit with the same settings; the scatter in V and sigma is added to the systematic budget (empty = off). Each library is compared in its own native wavelength frame.

### `snr_min_kinematics` — Min S/N per A for sigma  *(advanced)*

*default:* `10.0`

Below this median S/N per Angstrom the velocity dispersion is reported but flagged unreliable (the redshift can still be good; sigma needs far more S/N).

### `auto_accept` — Auto-select best candidate (headless)  *(advanced)*

*default:* `True`

In non-interactive runs, adopt the best-scoring candidate redshift that is not geometrically impossible for its role. In the GUI you choose from the candidate list.


## 6. Output

### `save_plots` — Save diagnostic plots

*default:* `True`

Write PNG figures for every stage.

### `save_fits_bundle` — Save FITS spectra bundle  *(advanced)*

*default:* `True`

Write all extracted spectra into one multi-extension FITS file.

### `save_ascii_spectra` — Save ASCII spectra  *(advanced)*

*default:* `True`

Write one text/ECSV file per aperture and arm.
