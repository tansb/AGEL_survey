from types import SimpleNamespace as N

import numpy as np

from kcwilensfit import plotting as P
from kcwilensfit.extraction import ArmSpectrum
from kcwilensfit.pipeline import Pipeline


def _arm(n=1900, seed=0):
    rng = np.random.default_rng(seed)
    w = np.linspace(3600, 5500, n)
    f = rng.normal(0, 1, n) + 8 * np.exp(-0.5 * ((w - 4500) / 3) ** 2)
    return w, f, ArmSpectrum(w, f, np.ones(n), 1.0)


def test_clean_blanks_spike_and_windows_but_keeps_line():
    w, f, a = _arm()
    a.flux[1850] += 500
    d, bad = P._clean_for_display(a, [(4000.0, 4010.0)])
    assert bad[1850]
    assert bad[(w >= 4000) & (w <= 4010)].all()
    assert np.isfinite(d.flux[np.argmin(abs(w - 4500))])
    assert np.isfinite(a.flux[1850])          # original untouched


def test_single_arm_xlim_ignores_windows_in_missing_arm():
    _, _, a = _arm()
    pipe = N(spectra={"A": N(arms={"blue": a}, role="source", meta={})}, results={})
    fig = P.fig_spectra(pipe, ["A"], windows=[(7580.0, 7700.0)])
    assert fig.axes[0].get_xlim()[1] < 5600


def test_combined_members_hidden_by_default():
    _, _, a = _arm()
    mk = lambda **m: N(arms={"blue": a}, role="source", meta=m)
    pipe = N(spectra={"A": mk(), "B": mk(), "G": mk(members=["A", "B"])}, results={})
    assert P.combined_member_ids(pipe) == {"A", "B"}
    assert len(P.fig_spectra(pipe).axes) == 1


def test_clear_fit_removes_only_that_result():
    p = Pipeline(log=lambda m: None)
    p.results = {"A": N(), "B": N()}
    assert p.clear_fit("A") is True
    assert list(p.results) == ["B"]
    assert p.clear_fit("A") is False


def test_quick_summary_lists_fits_and_unfitted(tmp_path):
    from kcwilensfit.report import quick_summary
    p = Pipeline(log=lambda m: None)
    _, _, a = _arm()
    p.spectra = {"A": N(arms={"blue": a}, role="source", meta={}), "B": N(arms={"blue": a}, role="source", meta={})}
    r = N(role="source", status="secure", z=1.3, z_err=1e-4, sigma=float("nan"), sigma_err=float("nan"),
          kinematics_ok=True, flags=[])
    p.results = {"A": r}
    txt = quick_summary(p, str(tmp_path))
    assert str(tmp_path) in txt and "1.30000" in txt and "B" in txt.split("No recorded fit")[1]


def test_fwhm_fields_blank_by_default_in_gui_spec():
    from kcwilensfit.schema import FIELD_BY_KEY, STAGES
    assert FIELD_BY_KEY["blue_fwhm_A"].blank_default and FIELD_BY_KEY["red_fwhm_A"].blank_default
    assert FIELD_BY_KEY["arm_registration"].stage == "inputs" and "instrument" not in dict(STAGES)


def test_unblanked_spectra_plot_keeps_spike_and_window_data():
    w, f, a = _arm()
    a.flux[1850] += 500
    pipe = N(spectra={"A": N(arms={"blue": a}, role="x", meta={})}, results={})
    fig = P.fig_spectra(pipe, ["A"], windows=[(4000.0, 4010.0)], despike=False, blank_windows=False)
    ydata = np.concatenate([ln.get_ydata() for ln in fig.axes[0].lines if len(ln.get_ydata()) == w.size])
    assert np.nanmax(ydata) > 400
    ln = next(l for l in fig.axes[0].lines if l.get_label() == "blue")
    assert np.isfinite(ln.get_ydata()[(w >= 4000) & (w <= 4010)]).all()
