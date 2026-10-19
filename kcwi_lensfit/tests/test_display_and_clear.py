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
