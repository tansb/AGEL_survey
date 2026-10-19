import numpy as np
from astropy.wcs import WCS

from kcwilensfit.apertures import Aperture, aperture_in_arm, aperture_mask, parse_ds9_regions


def make_wcs(crval=(150.0, 2.0), crpix=(40.5, 40.5), scale=0.3, rot_deg=0.0):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = list(crval)
    w.wcs.crpix = list(crpix)
    s = scale / 3600.0
    c, sn = np.cos(np.deg2rad(rot_deg)), np.sin(np.deg2rad(rot_deg))
    w.wcs.cd = [[-s * c, s * sn], [-s * sn, -s * c * -1]]
    w.wcs.cd = np.array([[-s * c, s * sn], [s * sn, s * c]])
    return w


def test_circle_mask_area():
    m = aperture_mask((60, 60), Aperture(id="a", shape="circle", x=30.2, y=29.7, a=5.0))
    assert abs(m.sum() - np.pi * 25) < 0.5
    assert m.max() <= 1.0


def test_box_and_ellipse_area():
    b = aperture_mask((60, 60), Aperture(id="b", shape="box", x=30, y=30, a=8, b=4, angle=30))
    e = aperture_mask((60, 60), Aperture(id="e", shape="ellipse", x=30, y=30, a=6, b=3, angle=45))
    assert abs(b.sum() - 32) < 0.8
    assert abs(e.sum() - np.pi * 18) < 0.6


def test_fk5_region_is_mapped_through_each_cubes_own_wcs(tmp_path):
    """The same sky region must land at DIFFERENT pixels in two cubes whose WCSs differ."""
    ref = make_wcs(crval=(150.0, 2.0))
    # red cube: pointing differs by 0.9" in RA and the grid is rotated by 2 degrees
    red = make_wcs(crval=(150.0 + 0.9 / 3600 / np.cos(np.deg2rad(2.0)), 2.0), rot_deg=2.0)
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    sc = ref.pixel_to_world(45.0, 38.0)
    reg = tmp_path / "w.reg"
    reg.write_text("# Region file format: DS9 version 4.1\nfk5\n"
                   f"circle({sc.ra.deg:.7f},{sc.dec.deg:.7f},1.5\") # text={{lensgal}}\n"
                   f"ellipse({sc.ra.deg:.7f},{sc.dec.deg:.7f},2.4\",1.2\",30)\n")
    aps, warns = parse_ds9_regions(str(reg), ref, 0.3)
    assert len(aps) == 2 and aps[0].role == "deflector" and aps[0].id == "lensgal"
    assert abs(aps[0].x - 45.0) < 1e-3 and abs(aps[0].y - 38.0) < 1e-3
    assert abs(aps[0].a - 5.0) < 1e-6                      # 1.5" / 0.3"/pix
    a_blue = aperture_in_arm(aps[0], ref, 0.3, ref)
    a_red = aperture_in_arm(aps[0], red, 0.3, ref)
    assert abs(a_blue.x - 45.0) < 1e-3
    assert abs(a_red.x - a_blue.x) > 2.0                  # 0.9" ~ 3 pixels
    # round trip: the red pixel position corresponds to the same sky position
    s2 = red.pixel_to_world(a_red.x, a_red.y)
    assert s2.separation(sc).arcsec < 1e-3
    # the ellipse's major axis must point in the same SKY direction in both cubes
    e_ref = aperture_in_arm(aps[1], ref, 0.3, ref)
    e_red = aperture_in_arm(aps[1], red, 0.3, ref)

    def sky_dir(w, ang):
        d = w.wcs.cd @ np.array([np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))])
        return np.degrees(np.arctan2(d[1], d[0]))

    diff = (sky_dir(red, e_red.angle) - sky_dir(ref, e_ref.angle) + 180) % 360 - 180
    assert abs(diff) < 0.05


def test_image_region_is_one_indexed(tmp_path):
    reg = tmp_path / "i.reg"
    reg.write_text("physical\ncircle(51.0,52.0,3)\n")
    aps, _ = parse_ds9_regions(str(reg), None, 0.3)
    assert aps[0].x == 50.0 and aps[0].y == 51.0 and not aps[0].has_sky
