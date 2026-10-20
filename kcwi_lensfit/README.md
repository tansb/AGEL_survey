# kcwilensfit — redshifts and kinematics of strong lenses from KCWI/KCRM cubes

An end-to-end, compound-lens-aware pipeline for **Keck/KCWI (blue) + KCRM (red)** integral-field
data of gravitational lenses. It takes reduced cubes and produces, for every deflector and every
lensed source image,

* a **secure/tentative/failed redshift** with statistical and systematic errors,
* the deflector's **stellar velocity dispersion** (with a systematic budget and a reliability flag),
* a lensing-geometry **consistency report** (is each source behind the lens? do all images of a
  source agree? how many redshift planes does a compound lens have?),

through a simple **tabbed GUI** (every input pre-filled with a documented default and described on
screen) or a one-line command-line run driven by the same JSON settings file.

```
 1 Data & instrument → 2 Noise & sky → 3 Apertures → 4 Extraction → 5 Redshift & kinematics → 6 Output
```

## Install

```bash
git clone <this repo> && cd kcwi_lensfit
conda env create -f environment.yml && conda activate kcwilensfit      # or: pip install -e .
kcwilensfit demo --run          # synthetic compound lens: checks your install end to end (~2 min)
```

Requirements: Python ≥ 3.9, numpy, scipy, astropy, matplotlib, [pPXF](https://pypi.org/project/ppxf/),
and Tkinter for the GUI (bundled with python.org / conda Python; `sudo apt install python3-tk` on Debian).
Stellar-population template libraries (E-MILES, FSPS, …) are downloaded once on first use into
`~/.cache/kcwilensfit` (or set `KCWILENSFIT_CACHE`); the Kinney–Calzetti empirical galaxy templates
ship with the package, so redshift-only use of blue star-forming sources works offline.

## Use

```bash
kcwilensfit gui                       # open the GUI (optionally: kcwilensfit gui my_config.json)
kcwilensfit run my_config.json        # headless: all stages, writes ./kcwilensfit_output
kcwilensfit init-config my.json       # a config with every default
kcwilensfit fields                    # every setting, default and description in the terminal
```

### GUI walk-through

| tab | what you do | what you see |
|---|---|---|
| 1 Data & instrument | choose blue/red cubes (either is optional), output folder, wavelength frame (air/vacuum), velocity correction, usable range, grating/slicer → FWHM (leave empty for the nominal value), blue/red registration → **Load cubes** | white-light images; registration overlay (red contours on blue, blue vs red aperture outlines) |
| 2 Noise & sky | auto/manual source-free sky box → **Estimate sky & noise** | background, noise and the noise-correlation law (σ_sum ∝ N^α) |
| 3 Apertures | table of deflectors/sources/others; load a **DS9 region file (pixel or WCS/fk5)**; auto-detect; click-to-place | the apertures **exactly as registered on each cube** |
| 4 Extraction | deblend vs sum, seeing, grouping of images → **Run extraction** | spectra, the spatial model & residuals |
| 5 Redshift & kinematics | **Scan** → pick a candidate (or type z) → **Fit**; or **Auto-fit ALL**; **Clear fit (this region)** discards one region's fit before saving | χ²(z)/line-score curves, pPXF fit, line fits, bootstrap, redshift ladder, results table with flags |
| 6 Output | review the summary, output directory, aperture and spectra+z plots, then **Save all outputs** | summary table, plots and full `report.md` preview *before* anything is written |
| | **Clear all (new target)** (bottom bar / File menu) resets everything for the next target | |

### Keeping the forms short

Each tab shows only the settings you would realistically change (31 of 85 in total); everything
else sits behind **“Show N advanced settings”** with a sensible default (e.g. template age and
metallicity cuts, polynomial orders, mask classes, bootstrap flavour, line-search internals).
All of them are still read from the config file and documented in `docs/config_reference.md`.
File/folder dialogs open in the folder you last used, including for the output directory.

### Regions in WCS coordinates

A DS9 region file saved in **fk5/icrs** is converted to pixels *separately for the blue and the
red cube using each cube's own header WCS* (centre, size and rotation), and the shapes actually
used are drawn on every preview. Image/physical-coordinate files are read as 1-indexed DS9 pixels of
the blue cube and shifted to the red cube with the registration of the Data tab. Turn the behaviour off
with *“Map sky-coordinate regions with each cube's WCS”* if you distrust the header zero-points
(they can be off by arcseconds; the registration overlay shows whether they agree).

### Outputs (`output_dir/`)

`report.md` (human-readable summary) · `redshifts.csv` / `results.json` · `spectra.fits` and
`spectra/*.ecsv` (wave, flux, variance, mask for every aperture and arm) · `figures/*.png` ·
`config_used.json` (re-run exactly with `kcwilensfit run output/config_used.json`).

## What it does, and why (short version; see `docs/science_notes.md`)

* **Air vs vacuum is tracked explicitly** (blue = air, red = vacuum in Keck DRP cubes); the wrong
  frame costs ~85 km/s.
* **Noise is measured, not assumed**: sky-box scatter or variance cubes, plus the covariance created
  by resampling the slicer data (errors of summed apertures grow as N^α, α>0.5).
* **Spectral deblending** (default): lens galaxy (PSF-convolved Sérsic), point-like images and arcs
  are fitted simultaneously at every wavelength, so overlapping deflector/arc/second-lens light is
  separated with honest, covariance-aware errors. Plain aperture sums remain available.
* **Redshift scan**: continuum χ²(z) is scanned for *several* fixed dispersions (default 150 and
  350 km/s); candidates present at every σ are marked stable, σ-dependent ones are flagged.
* **Deflectors**: pPXF joint blue+red fit (one z and σ, free red/blue scale, multiplicative
  polynomials, frame-exact wavelengths, resolution matched in the physically correct direction),
  rest window restricted to where the templates are sharp, bootstrap errors, polynomial/library sweep
  as a systematic, blue-only vs red-only and stellar-absorption-line consistency checks.
* **Sources**: null-subtracted matched-filter search over emission *and* absorption multiplets
  (no initial guess needed), joint line fit with separate offsets for nebular, ISM-absorption and
  resonant-emission classes, an explicit systemic-velocity systematic when only outflow-affected
  lines are seen, and a **competing-identification table** for single-line detections.
* **Compound lenses**: any number of apertures; images of one source share a *group* and are
  combined (achromatic lensing); once one plane's redshift is known its emission lines are masked
  in the other planes' fits (e.g. a z=1.3 source's [O II] inside the deflector's window); the report
  lists the redshift ladder and flags geometry violations.

## Honest limitations

* No telluric correction (bands are masked); no flux calibration beyond what the DRP provides.
* Source photon noise is not in the noise model (pPXF errors are rescaled by √χ²/dof when χ² > 1).
* Velocity dispersions need S/N ≳ 10 Å⁻¹ and an aperture free of contaminating arc light; quote
  σ_e only with the aperture and aperture-correction conventions of your paper (see notes).
* Region angles from sky-coordinate files assume a (near) north-up cube.
* Validated on synthetic compound-lens cubes with known truth (`tests/`, `kcwilensfit demo`) and run
  on real AGEL cubes; please report problems with a copy of `config_used.json` and the log.

## Tests

```bash
pip install -e .[test] && pytest -q            # fast unit tests
pytest -q -m slow                              # full synthetic compound-lens recovery
```

## Citation / credits

Please cite pPXF (Cappellari & Emsellem 2004; Cappellari 2017, 2023), the E-MILES library
(Vazdekis et al. 2016) if used, and Kinney et al. (1996) for the empirical templates.
