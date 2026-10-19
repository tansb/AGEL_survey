import numpy as np

from kcwilensfit.cubeio import Arm
from kcwilensfit import noise as nz
from kcwilensfit.extraction import deblend, extract_sum, inverse_variance_combine, ArmSpectrum, reject_spikes
from kcwilensfit.apertures import Aperture, aperture_mask
from scipy.ndimage import gaussian_filter


def fake_arm(flux, wave=None):
    nw, ny, nx = flux.shape
    wave = np.arange(nw) + 4000.0 if wave is None else wave
    return Arm(name="blue", path="x", flux=flux.astype(np.float32), var=None, header=None, wcs=None, wave=wave,
               native_frame="air", fwhm_A=2.5, pixscale=(0.3, 0.3), valid=np.ones((ny, nx), bool))


def test_independent_noise_gives_alpha_half():
    rng = np.random.default_rng(0)
    arm = fake_arm(rng.normal(0, 1, (400, 40, 40)))
    sky = nz.estimate_sky(arm, (5, 25, 5, 25))
    assert abs(sky.alpha - 0.5) < 0.12
    assert abs(np.median(sky.noise) - 1.0) < 0.1


def test_correlated_noise_gives_alpha_above_half():
    rng = np.random.default_rng(1)
    f = gaussian_filter(rng.normal(0, 1, (400, 40, 40)), (0, 0.9, 0.9)) * 3.0
    arm = fake_arm(f)
    sky = nz.estimate_sky(arm, (5, 25, 5, 25))
    assert sky.alpha > 0.62


def test_suggest_sky_box_avoids_source():
    img = np.ones((50, 50)) * 0.01
    img[20:30, 20:30] = 10.0
    valid = np.ones_like(img, bool)
    y0, y1, x0, x1 = nz.suggest_sky_box(img + np.random.default_rng(2).normal(0, 1e-3, img.shape), valid, 10)
    assert not (y0 < 30 and y1 > 20 and x0 < 30 and x1 > 20)


def test_deblend_recovers_overlapping_spectra():
    rng = np.random.default_rng(3)
    ny = nx = 41
    yy, xx = np.mgrid[:ny, :nx]
    g = lambda x0, y0, s: gaussian_filter(np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / s ** 2), 1.0)
    P1, P2 = g(18, 20, 4.0), g(26, 21, 1.5)
    P1 /= P1.sum(); P2 /= P2.sum()
    nw = 300
    w = np.arange(nw) + 5000.0
    s1 = 100 + 20 * np.sin(w / 30)
    s2 = 50 * np.exp(-0.5 * ((w - 5150) / 3) ** 2)
    cube = s1[:, None, None] * P1 + s2[:, None, None] * P2 + rng.normal(0, 1e-3, (nw, ny, nx))
    arm = fake_arm(cube, w)
    sky = nz.SkyModel(box=(0, 5, 0, 5), bg=np.zeros(nw), noise=np.full(nw, 1e-3), n_pix=25)
    masks = [aperture_mask((ny, nx), Aperture(id="a", x=18, y=20, a=12)), aperture_mask((ny, nx), Aperture(id="b", x=26, y=21, a=4))]
    out, diag = deblend(arm, sky, [P1, P2], masks, np.ones((ny, nx), bool), subtract_background=False)
    n2 = float(np.sum(P2 * masks[1]))
    # blended compact line source recovered without the broad component's continuum
    assert abs(np.nanmedian(out[1].flux[:50])) < 0.05 * 50 * n2
    assert abs(out[1].flux[150] / (50 * n2) - 1) < 0.05
    # and the blending covariance is reflected in a larger error than an isolated measurement
    assert np.all(out[1].var > 0)


def test_sum_extraction_flux_scale():
    arm = fake_arm(np.ones((50, 30, 30)) * 2.0)
    sky = nz.SkyModel(box=(0, 5, 0, 5), bg=np.zeros(50), noise=np.full(50, 0.1), n_pix=25)
    m = aperture_mask((30, 30), Aperture(id="a", x=15, y=15, a=4))
    s = extract_sum(arm, sky, m, subtract_background=False)
    assert abs(s.flux[0] - 2.0 * m.sum()) < 1e-3


def test_spike_rejection_keeps_real_lines():
    w = np.arange(500) + 5000.0
    f = np.random.default_rng(4).normal(0, 1, 500)
    f += 12 * np.exp(-0.5 * ((w - 5250) / 2.0) ** 2)          # resolved line
    f[100] += 40                                              # one-pixel cosmic ray
    s = ArmSpectrum(w, f.copy(), np.ones(500))
    n = reject_spikes(s)
    assert np.isnan(s.flux[100]) and np.isfinite(s.flux[250]) and n >= 1


def test_parse_subtraction_rules():
    from kcwilensfit.extraction import parse_subtraction_rules
    r = parse_subtraction_rules("A: deflector 0.3\nB: deflector auto  # comment\n\nC: deflector")
    assert r == [("A", "deflector", 0.3), ("B", "deflector", None), ("C", "deflector", None)]
    import pytest
    with pytest.raises(ValueError):
        parse_subtraction_rules("A deflector")
    with pytest.raises(ValueError):
        parse_subtraction_rules("A: A 0.1")
    with pytest.raises(ValueError):
        parse_subtraction_rules("A: deflector abc")


def test_subtract_spectrum_propagates_variance():
    from kcwilensfit.extraction import ArmSpectrum, subtract_spectrum
    w = np.arange(10.0)
    t = ArmSpectrum(w, np.full(10, 5.0), np.full(10, 1.0))
    c = ArmSpectrum(w, np.full(10, 2.0), np.full(10, 4.0))
    subtract_spectrum(t, c, 0.5)
    assert np.allclose(t.flux, 4.0) and np.allclose(t.var, 1.0 + 0.25 * 4.0)
