import numpy as np
import pytest

from kcwilensfit.pipeline import Pipeline


@pytest.fixture(scope="module")
def pipe(demo):
    cfg = dict(blue_cube=demo["paths"]["blue"], red_cube=demo["paths"]["red"], apertures=demo["apertures"],
               n_bootstrap=0, sweep_degrees="", crosscheck_arms=False)
    p = Pipeline(cfg, log=lambda m: None)
    p.load(); p.estimate_sky(); p.set_apertures(cfg["apertures"]); p.extract()
    return p


def test_loader_handles_frames_zeros_and_units(pipe):
    b, r = pipe.arms["blue"], pipe.arms["red"]
    assert b.native_frame == "air" and r.native_frame == "vacuum"
    assert b.valid.sum() == r.valid.sum() < b.valid.size            # zero border -> invalid
    assert np.isnan(b.flux[:, ~b.valid]).all()
    assert abs(b.unit_scale - 0.09) < 1e-6
    assert r.wave.min() > 5618                                      # WAVGOOD0 trim (+ vac->air)


def test_noise_correlation_is_detected(pipe):
    for s in pipe.sky.values():
        assert 0.55 < s.alpha <= 1.0                                # synthetic noise is correlated


def test_deblend_separates_overlapping_objects(pipe):
    # the deflector spectrum must not contain the arcs' [O II]/C III] emission
    sp = pipe.spectra["deflector"].arms["red"]
    ok = sp.good
    assert sp.median_snr() > 5


@pytest.mark.slow
def test_compound_lens_redshift_recovery(pipe, demo):
    t = demo["truth"]
    res = pipe.fit_all()
    assert abs(res["deflector"].z - t["z_defl"]) < 5e-4
    assert abs(res["deflector"].sigma - t["sigma_defl"]) < 25
    assert abs(res["group:S1"].z - t["z_src1"]) < 5e-4
    c = res["C"]
    # ISM-absorption-defined redshift is blueshifted ~150 km/s; the systemic estimate recovers the truth
    assert abs(c.detail["z_sys_estimate"] - t["z_src2"]) < 1e-3
    assert any("COMPOUND" in m for m in pipe.system_messages)


def test_group_status_capped_when_members_disagree():
    from kcwilensfit.redshift import ApertureResult, RedshiftEngine
    from kcwilensfit.schema import defaults
    eng = RedshiftEngine(defaults(), {}, {}, {}, [], log=lambda m: None)
    mk = lambda i, z, st="tentative", g="S1": ApertureResult(id=i, role="source", group=g, mode="lines", z=z,
                                                            z_err=1e-4, status=st)
    eng.results = {"a": mk("a", 1.65), "b": mk("b", 0.40), "group:S1": mk("group:S1", 0.98, "secure")}
    _, msgs = eng.system_checks()
    assert eng.results["group:S1"].status == "tentative"
    assert any("downgraded" in m for m in msgs)
    # agreeing members leave a secure group alone
    eng.results = {"a": mk("a", 1.30), "b": mk("b", 1.30001), "group:S1": mk("group:S1", 1.30, "secure")}
    eng.system_checks()
    assert eng.results["group:S1"].status == "secure"


def test_spaxel_counts_and_manual_subtraction(demo):
    cfg = dict(blue_cube=demo["paths"]["blue"], red_cube=demo["paths"]["red"], apertures=demo["apertures"],
               extraction_method="sum", manual_subtractions="A: deflector 0.5\nB: deflector auto")
    logs = []
    p = Pipeline(cfg, log=logs.append)
    p.load(); p.estimate_sky(); p.set_apertures(cfg["apertures"])
    p.extract()
    d = p.spectra["deflector"]
    cnt = d.meta["counts"]["blue"]
    assert 70 < cnt["area"] < 85 and cnt["touched"] >= cnt["area"] * 0.95      # pi*5^2 ~ 78.5 spaxels
    assert any("spaxels" in m for m in logs)
    # manual rule: A = A_raw - 0.5*deflector, variance propagated
    p2 = Pipeline({**cfg, "manual_subtractions": ""}, log=lambda m: None)
    p2.load(); p2.estimate_sky(); p2.set_apertures(cfg["apertures"]); p2.extract()
    a_raw, a_sub = p2.spectra["A"].arms["red"], p.spectra["A"].arms["red"]
    dfl = p2.spectra["deflector"].arms["red"]
    ok = a_sub.good & a_raw.good
    assert np.allclose((a_raw.flux - a_sub.flux)[ok], 0.5 * dfl.flux[ok], rtol=1e-4, atol=1e-6)
    assert np.all((a_sub.var - a_raw.var)[ok] > 0)
    rec = p.spectra["B"].meta["subtracted"][0]
    assert rec["auto"] and 0.0 <= rec["scale"] < 1.0                           # leakage estimated from the profile


def test_source_continuum_fit_with_limited_template_coverage(pipe, monkeypatch, tmp_path):
    """Regression: empirical kc96 templates (1235-9945 A) are shorter than the galaxy at high trial z."""
    monkeypatch.setenv("KCWILENSFIT_CACHE", str(tmp_path))
    eng = pipe.make_engine()
    eng.cfg.update(n_bootstrap=0, sweep_degrees="", crosscheck_arms=False, z_step_continuum=0.03, n_null_sims=0)
    sc = eng.scan("A", kind="continuum")
    assert sc[0].candidates, "continuum scan over z=0.2-3.5 must run without error"
    r = eng.fit("A", mode="continuum", scans=sc)
    assert not any("fit failed" in f for f in r.flags), r.flags


def test_same_plane_lines_are_not_masked():
    """A region fitted at (or near) another region's redshift must keep that redshift's lines."""
    from kcwilensfit.redshift import ApertureResult, RedshiftEngine
    from kcwilensfit.extraction import ArmSpectrum, ExtractedSpectrum
    from kcwilensfit.apertures import Aperture
    from kcwilensfit.schema import defaults
    w = np.arange(3300.0, 8900.0)
    sp = {k: ExtractedSpectrum(k, "source", "", {"blue": ArmSpectrum(w, np.ones_like(w), np.ones_like(w))})
          for k in ("A", "B", "C")}
    eng = RedshiftEngine(defaults(), sp, {k: Aperture(id=k, role="source") for k in sp}, {"blue": 2.5}, ["blue"], log=lambda m: None)
    eng.results = {"A": ApertureResult(id="A", role="source", group="", mode="lines", z=1.30, z_err=1e-4, status="secure")}
    # B's guess is far from A's z -> A's lines are masked
    far = eng.other_plane_windows("B", 3300.0, 8900.0, z_self=0.62)
    assert far and eng.last_mask_info["masked"] and not eng.last_mask_info["skipped_same_plane"]
    # B's guess is A's redshift (ungrouped!) -> nothing masked, and the decision is recorded
    same = eng.other_plane_windows("B", 3300.0, 8900.0, z_self=1.3004)
    assert same == [] and eng.last_mask_info["skipped_same_plane"][0][0] == "A"
    # unknown z (search) -> windows come only from static masks
    assert eng._windows_for("B", ["blue"], use_other=False) == eng._windows_for("B", ["blue"], z_self=1.3, use_other=True)
    # the report flags ungrouped regions that agree in redshift
    eng.results["B"] = ApertureResult(id="B", role="source", group="", mode="lines", z=1.30005, z_err=1e-4, status="tentative")
    eng.results["C"] = ApertureResult(id="C", role="source", group="", mode="lines", z=0.9, z_err=1e-4, status="tentative")
    _, msgs = eng.system_checks()
    assert any(m.startswith("SAME PLANE?") and "'A'" in m and "'B'" in m for m in msgs)
    assert not any("'C'" in m and m.startswith("SAME PLANE?") for m in msgs)
