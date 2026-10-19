"""Wavelength-frame helpers: air <-> vacuum and velocity <-> redshift conversions.

Why this matters for a redshift pipeline
----------------------------------------
Air and vacuum wavelengths differ by ~ (n-1) ~ 2.8e-4, i.e. ~85 km/s. If a
spectrum is compared against rest-frame templates or line lists in the *other*
frame, that offset lands one-for-one in the fitted systemic velocity.  The
KCWI blue cubes distributed by the Keck DRP are on an air wavelength axis
(``CTYPE3 = AWAV``) while the KCRM red cubes are vacuum (``CTYPE3 = WAVE``);
this package therefore tracks the frame of every wavelength array explicitly
and converts with the (Ciddor 1996 / Morton 2000) dispersion formula used by
pPXF itself.
"""
from __future__ import annotations

import numpy as np

C_KMS = 299792.458  # speed of light [km/s]


def _refractive_index_term(sigma2):
    """(n - 1) for standard air, Morton (2000) form of Ciddor (1996); sigma = 1e4/lambda[A]."""
    return 0.0000834254 + 0.02406147 / (130.0 - sigma2) + 0.00015998 / (38.9 - sigma2)


def vac_to_air(lam_vac):
    """Vacuum -> air wavelength [Angstrom]."""
    lam_vac = np.asarray(lam_vac, dtype=float)
    s2 = (1.0e4 / lam_vac) ** 2
    return lam_vac / (1.0 + _refractive_index_term(s2))


def air_to_vac(lam_air):
    """Air -> vacuum wavelength [Angstrom] (Morton 2000 / Ciddor 1996)."""
    lam_air = np.asarray(lam_air, dtype=float)
    s2 = (1.0e4 / lam_air) ** 2
    return lam_air * (1.0 + _refractive_index_term(s2))


def velocity_to_z(v_kms, z_ref=0.0):
    """Relativistic-consistent log-lambda velocity -> redshift about z_ref.

    pPXF reports velocities as c*ln(lambda ratio), so the redshift is
    (1+z) = (1+z_ref) * exp(v/c)."""
    return (1.0 + z_ref) * np.exp(np.asarray(v_kms, dtype=float) / C_KMS) - 1.0


def z_to_velocity(z, z_ref=0.0):
    return C_KMS * np.log((1.0 + np.asarray(z, dtype=float)) / (1.0 + z_ref))


def dz_to_kms(dz, z):
    """Velocity difference [km/s] corresponding to a redshift difference at redshift z."""
    return C_KMS * np.asarray(dz, dtype=float) / (1.0 + np.asarray(z, dtype=float))


def kms_to_dz(dv_kms, z):
    return np.asarray(dv_kms, dtype=float) * (1.0 + np.asarray(z, dtype=float)) / C_KMS


def apply_velocity_correction(wave, v_corr_kms):
    """Shift an observed wavelength axis by a velocity correction.

    Convention: ``v_corr`` is the velocity to ADD to the observed radial velocity to
    obtain the heliocentric/barycentric one (the sign convention of the KCWI DRP
    ``VCORR`` keyword), so wavelengths are scaled by (1 + v_corr/c)."""
    return np.asarray(wave, dtype=float) * (1.0 + v_corr_kms / C_KMS)
