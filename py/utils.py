import numpy as np


def get_omgH(a):
    return a / (2.0 * (1.0 + np.sqrt(1.0 - a**2)))


def get_rc(r, a):
    return np.sqrt(r**2 + a**2 + 2 * a**2 / r)


def get_vr(r, prstar, pph, a):
    rc = get_rc(r, a)
    Ap = (1.0 - 2.0 / rc) * (1.0 + 2.0 / rc) / (1.0 + 2.0 / r)
    Horb = np.sqrt(Ap * (1.0 + pph**2 / rc**2 + prstar**2))
    return Ap * rc / r * prstar / Horb


def adv_time_in(a, r0, Shh=10):
    """
    Calculate initial advanced time at horizon for orbit with initial radius r0
    Teukode shifts HH time so it starts at 0; subtract this from HH time to compare with (BL) dynamics
    From Sec. 2.2-3 of https://arxiv.org/pdf/1406.5983
    """
    # Initial HH radial coord.
    x0 = r0 / (1.0 + r0 / Shh)
    # Inner and outer horizons
    rplus, rminus = 1.0 + np.sqrt(1.0 - a**2), 1.0 - np.sqrt(1.0 - a**2)
    # HH horizon (outer)
    xplus = rplus / (1.0 + rplus / Shh)
    # Initial tortoise
    rs0 = r0 + (rplus * np.log(r0 - rplus) - rminus * np.log(r0 - rminus)) / np.sqrt(
        1.0 - a**2
    )
    # Initial HH coord. time
    tau0 = -r0 + rs0 - x0 / (1 - x0 / Shh) + x0 + 4.0 * np.log(1.0 - x0 / Shh)
    # v - tau, at horizon
    v_delta = (
        xplus * (Shh + xplus) / (Shh - xplus)
        - 4.0 * np.log(1.0 - xplus / Shh)
        - 2.0 * np.log(2.0)
    )

    return v_delta + tau0


def pph_circ(a, r, retro=False):
    """
    Angular momentum of a circular orbit of radius r in Kerr geometry with spin a.
    retro = True for retrograde orbits.
    From Chandrasekhar, Chapter 61.b
    """
    u = 1.0 / r
    sign = -1 if retro else 1
    Q = 1.0 - 3.0 * u + sign * 2.0 * a * u ** (1.5)
    return sign * np.sqrt(r / Q) * (1.0 + a * a * u * u - sign * 2.0 * a * u ** (1.5))


def OmegaTS(r, omega, rdot, nu, a, model="base", emm=None):
    """
    Tidal angular frequency, according to some model
    """

    m1 = 0.5 * (1.0 + np.sqrt(1.0 - 4 * nu))
    m2 = 0.5 * (1.0 - np.sqrt(1.0 - 4 * nu))

    if model == "base":
        return omega
    elif model == "ns":
        if emm is None or 2 in emm:
            return omega * (
                1.0 - nu / r + (1.0 - m1 - nu) * rdot**2 - 16.0 / 3.0 * m1 * rdot / r
            )
        else:
            return omega
    elif model == "qc":
        return (
            omega
            + omega
            * (-1.0 + m1**2 + m2**2 * r * rdot**2 + m2**2 * r**3 * omega**2)
            / r
            / 2.0
            + m2 / r**3 * (m1 * a)
        )


def OmegaTE(r, omega, rdot, nu, a, model="base", emm=None):
    """
    Tidal angular frequency, according to some model
    """

    m1 = 0.5 * (1.0 + np.sqrt(1.0 - 4 * nu))
    m2 = 0.5 * (1.0 - np.sqrt(1.0 - 4 * nu))

    m_coeff = 3.0 if emm is None or (2 in emm and 0 in emm) else 9.0 / 4
    if emm == [1]:
        m_coeff = 0
    if model == "base":
        return (
            OmegaTS(r, omega, rdot, nu, a, model, emm)
            + m_coeff * rdot**2 / r**2 / omega
        )
    elif model == "ns":
        if emm == [1]:
            return (
                OmegaTS(r, omega, rdot, nu, a, model, emm)
                + 16.0 * rdot**2 / r**2 / omega
            )
        elif emm == [2]:
            romg_coeff = 7.0 / 3.0
        else:
            romg_coeff = 25.0 / 3.0
        return OmegaTS(r, omega, rdot, nu, a, model, emm) + m_coeff / r**2 / omega * (
            rdot**2
            * (
                1.0
                + (2.0 - 11.0 / 3 * m1 - 2.0 * nu) / r
                + r**2 * omega**2 * (romg_coeff + m1 + nu)
            )
            - 16.0 / 3 * m1 * (1.0 / r**2 - r * omega**2) * rdot
        )
    elif model == "qc":
        return OmegaTS(r, omega, rdot, nu, a, model, emm) + 3.0 * rdot**2 / r**2 / omega


def OmegaTS_p(r, pph, pr, nu, a, model="base", emm=None, newfact=None):
    """
    Tidal angular frequency, according to some model, as function of coordinates and momenta
    """

    m1 = 0.5 * (1.0 + np.sqrt(1.0 - 4 * nu))
    m2 = 0.5 * (1.0 - np.sqrt(1.0 - 4 * nu))

    if newfact is not None:
        K = pph**2 / pph_circ(a, r, retro=False) ** 2

    if "base" in model:
        if newfact is not None:
            return np.sqrt(K / r**3)
        else:
            return pph / r**2
    elif "ns" in model:
        if newfact is not None:
            if emm is None or 2 in emm:
                return np.sqrt(K / r**3) * (
                    1.0
                    + 0.5 / r * (1.0 - K * (1.0 + nu))
                    + 0.5 * pr**2 * (1.0 - 2.0 * m1 - 3.0 * nu)
                    - 16.0 / 3.0 * m1 * pr / r
                )
            else:
                return np.sqrt(K / r**3)
        else:
            if emm is None or 2 in emm:
                return (
                    pph
                    / r**2
                    * (
                        1.0
                        - 0.5 * (1.0 + nu) * pph**2 / r**2
                        - 1.0 / r
                        + 0.5 * (1.0 - 2 * m1 - 3.0 * nu) * pr**2
                        - 16.0 * m1 * pr / r
                    )
                )
            else:
                return pph / r**2


def OmegaTE_p(r, pph, pr, nu, a, model="base", emm=None, newfact=None):
    """
    Tidal angular frequency, according to some model, as function of coordinates and momenta
    """

    m1 = 0.5 * (1.0 + np.sqrt(1.0 - 4 * nu))
    m2 = 0.5 * (1.0 - np.sqrt(1.0 - 4 * nu))

    if newfact is not None:
        K = pph**2 / pph_circ(a, r, retro=False) ** 2

    m_coeff = 3.0 if emm is None or (2 in emm and 0 in emm) else 9.0 / 4
    if emm == [1]:
        m_coeff = 0
    if "base" in model:
        if newfact is not None:
            return OmegaTS_p(
                r, pph, pr, nu, a, model, emm, newfact=newfact
            ) + m_coeff * pr**2 / np.sqrt(K * r)
        else:
            return (
                OmegaTS_p(r, pph, pr, nu, a, model, emm, newfact=newfact)
                + m_coeff * pr**2 / pph
            )
    elif "ns" in model:
        if newfact is not None:
            if emm == [1]:
                return OmegaTS_p(
                    r, pph, pr, nu, a, model, emm, newfact=newfact
                ) + 16.0 * pr**2 / np.sqrt(K * r)
            elif emm == [2]:
                romg_coeff = 11.0 / 6.0
            else:
                romg_coeff = 47.0 / 6.0
            return OmegaTS_p(
                r, pph, pr, nu, a, model, emm, newfact=newfact
            ) + m_coeff / np.sqrt(K * r) * (
                pr**2
                * (
                    1.0
                    - 0.5 * (1.0 + nu) * pr**2
                    + (romg_coeff + 6.0 * m1 + 3.0 * nu) * K / r
                    + (-4.5 - 11.0 / 3.0 * m1 - nu) / r
                )
                + 16.0 / 3.0 * m1 * pr / r**2 * (K - 1.0)
            )
        else:
            if emm == [1]:
                return (
                    OmegaTS_p(r, pph, pr, nu, a, model, emm, newfact=newfact)
                    + 16.0 * pr**2 / pph
                )
            elif emm == [2]:
                romg_coeff = 11.0 / 6.0
            else:
                romg_coeff = 47.0 / 6.0
            return OmegaTS_p(
                r, pph, pr, nu, a, model, emm, newfact=newfact
            ) + m_coeff / pph * (
                pr**2
                * (
                    1.0
                    - 0.5 * (1.0 + nu) * pr**2
                    + (romg_coeff + 6.0 * m1 + 3.0 * nu) * pph**2 / r**2
                    + (-3.0 - 11.0 / 3.0 * m1 - nu) / r
                )
                + 16.0 / 3.0 * m1 * pr / r**3 * (pph**2 - r)
            )
