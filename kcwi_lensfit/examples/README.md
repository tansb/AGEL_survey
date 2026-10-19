# Examples

* `example_config.json` – a template with every setting at its default; fill in `blue_cube`,
  `red_cube` and the `apertures`, then `kcwilensfit run example_config.json` (or open it in the GUI).
* `lens_regions_fk5.reg` – a region file in WCS coordinates: load it in tab 4 (“Load DS9 region
  file…”); each cube's own WCS is used to place it. Give regions a `# text={name}` label to name
  the apertures (labels containing “defl”/“lens” become deflectors).
* `kcwilensfit demo --run` – a fully synthetic compound lens (deflector at z=0.62, source planes at
  z=1.30 and z=2.09) with a known answer.
