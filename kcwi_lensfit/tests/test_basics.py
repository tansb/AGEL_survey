import numpy as np
import pytest

from kcwilensfit import schema
from kcwilensfit.lines import LINES_BY_NAME, expected_lines, observed_wavelengths_air
from kcwilensfit.wavelength import air_to_vac, vac_to_air, z_to_velocity, velocity_to_z


def test_air_vac_roundtrip_and_size():
    lam = np.array([3500.0, 5000.0, 8000.0])
    assert np.allclose(vac_to_air(air_to_vac(lam)), lam, atol=1e-3)
    # ~ 1.4 A at 5000 A (n-1 ~ 2.8e-4)
    assert 1.2 < air_to_vac(5000.0) - 5000.0 < 1.6


def test_velocity_redshift_roundtrip():
    z = 0.62
    assert abs(velocity_to_z(z_to_velocity(0.63, z), z) - 0.63) < 1e-12


def test_line_observed_wavelength_known_value():
    # [O III] 5008.24 (vac) at z=1 -> vac 10016.5 -> air ~ 10013.7
    lam = observed_wavelengths_air(LINES_BY_NAME["[O III] 5007"], 1.0)
    assert abs(lam[1] - 10013.7) < 1.0
    # H-alpha at z=0.2 falls in the red arm
    assert expected_lines(0.2, 5600, 8900, ["H-alpha"])


def test_schema_defaults_validate_and_roundtrip(tmp_path):
    cfg = schema.defaults()
    cfg["blue_cube"] = ""
    probs = schema.validate(cfg)
    assert any("at least one" in p for p in probs)
    p = tmp_path / "c.json"
    schema.save_config(cfg, str(p))
    cfg2 = schema.load_config(str(p))
    assert cfg2["mdegree"] == cfg["mdegree"]
    assert schema.parse_range("-1.0-0.3") == (-1.0, 0.3)
    assert schema.parse_range("0.5-14") == (0.5, 14.0)
    assert schema.parse_range("0.1,1.2") == (0.1, 1.2)


def test_every_field_has_help_and_default():
    for f in schema.FIELDS:
        assert f.help and len(f.help) > 20, f.key
        assert f.stage in dict(schema.STAGES)
        if f.kind == "choice":
            assert str(f.default) in f.choices, f.key


def test_output_dir_falls_back_when_not_writable(tmp_path, monkeypatch):
    import os
    from kcwilensfit.report import writable_output_dir
    ro = tmp_path / "ro"
    ro.mkdir()
    os.chmod(ro, 0o500)
    monkeypatch.setenv("HOME", str(tmp_path))
    msgs = []
    try:
        out = writable_output_dir(str(ro), "AGEL", msgs.append)
    finally:
        os.chmod(ro, 0o700)
    if os.access(ro, os.W_OK):          # running as a user that ignores permissions (e.g. root): nothing to test
        return
    assert out != str(ro) and out.startswith(str(tmp_path)) and msgs
    assert writable_output_dir(str(tmp_path / "new"), "x", msgs.append) == str(tmp_path / "new")


def test_every_empirical_library_builds_and_loads(tmp_path, monkeypatch):
    """Regression: building a kc96 library once crashed (np.arange got a type as its step)."""
    monkeypatch.setenv("KCWILENSFIT_CACHE", str(tmp_path))
    import numpy as np
    from kcwilensfit import constants as K
    from kcwilensfit.templates import resolve_library
    for name in list(K.KC96_LIBRARIES) + list(K.KC96_TEMPLATES)[:2]:
        lib = resolve_library(name)
        d = np.load(lib.path)
        assert lib.is_empirical and d["templates"].shape[1] == len(d["ages"]) and d["lam"].size > 100
