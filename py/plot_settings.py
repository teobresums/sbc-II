"""
Plot settings for consistent styling across plots.
Import it in your plotting scripts.

Example:
from plot_settings import *
"""

import matplotlib.pyplot as plt
import seaborn as sns

latex_fonts = {
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
}
plt.rcParams.update(latex_fonts)

# names for keys used in data files
keys_to_labels = {
    "a": r"$\hat{{a}}$",
    "logE": r"$\log_{10} E_0$",
    "logPphi": r"$\log_{10} p_{\varphi 0}$",
    "logdE": r"$\log_{10} |\dot{M}|_{\rm peak}$",
    "logdJ": r"$\log_{10} |\dot{S}|_{\rm peak}$",
    "dE": r"$ \dot{M}_{\rm peak}$",
    "dJ": r"$ \dot{S}_{\rm peak}$",
    "rmin": r"$r_{\rm min}$",
    "e0": r"$e_0$",
    "p0": r"$p_0$",
    "E": r"$E_0$",
    "Pphi": r"$p_{\varphi 0}$",
}

# color settings and utilities
colors = [
    "#332288",  # Deep Indigo
    "#88CCEE",  # Sky Blue
    "#44AA99",  # Teal
    "#D64C7F",  # Sand
    "#DDCC77",  # Muted Gold
]

# colors = [
#     "#C84B31",  # Burnt Sienna
#     "#E8A838",  # Amber
#     "#2B6777",  # Deep Teal
#     "#52AB98",  # Sage Green
#     "#8E5EA2",  # Muted Purple
# ]

# mapping orbit types to colors, markers and "kind" values
colors_orbit_types = {
    "circular": (colors[3], "o", 0),
    "eccentric": (colors[2], "s", 1),
    "hyperbolic": (colors[0], "^", 2),
}


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


# combinations for PN orders and their names
pn_orders_combinations = [
    (
        "LO",
        "LO",
        r"\dot{E}_\mathrm{PN}^\mathrm{LO}",
        r"\dot{J}_\mathrm{PN}^\mathrm{LO}",
    ),
    (
        "NLO",
        "NLO",
        r"\dot{E}_\mathrm{PN}^\mathrm{NLO}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NLO}",
    ),
    (
        "NNLO",
        "NNLO",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO}",
    ),
    (
        "fujita",
        "NNLO NC + Fujita C",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,Fuj}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,Fuj}",
    ),
    (
        "NNLO_fact",
        "NNLO fact",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact}",
    ),
    (
        "NNLO_fact_omg",
        r"NNLO fact, $\omega$",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega}",
    ),
    (
        "NNLO_fact_rdot",
        r"NNLO fact, $\dot{r}$",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r}}}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r}}}",
    ),
    (
        "NNLO_fact_CF",
        "NNLO fact, CF",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF}",
    ),
    (
        "NNLO_CF",
        "NNLO, CF",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,CF}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,CF}",
    ),
    (
        "NNLO_fact_CF_rdot",
        r"NNLO fact, CF, $\dot{r}$",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,\dot{r}}}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,\dot{r}}}",
    ),
    (
        "NNLO_fact_CF_omg",
        r"NNLO fact, CF, $\omega$",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,\omega}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,\omega}",
    ),
    (
        "NNLO_CF_CP12",
        r"NNLO, CF, P12",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,CF,\,CP12}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,CF,\,CP12}",
    ),
    (
        "NNLO_CF_CP21",
        r"NNLO, CF, P21",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,CF,\,CP21}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,CF,\,CP21}",
    ),
    (
        "NNLO_CF_CP03",
        r"NNLO, CF, P03",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,CF,\,CP03}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,CF,\,CP03}",
    ),
    (
        "NNLO_fact_CF_CP12",
        r"NNLO fact, CF, P12",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,CP12}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,CP12}",
    ),
    (
        "NNLO_fact_CF_CP21",
        r"NNLO fact, CF, P21",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,CP21}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,CP21}",
    ),
    (
        "NNLO_fact_CF_CP03",
        r"NNLO fact, CF, P03",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,CP03}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,CF,\,CP03}",
    ),
    (
        "NNLO_fact_rdot_CF_CP12",
        r"NNLO fact, $\dot{r}$, CF, P12",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r},\,CF,\,CP12}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r},\,CF,\,CP12}",
    ),
    (
        "NNLO_fact_rdot_CF_CP21",
        r"NNLO fact, $\dot{r}$, CF, P21",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r},\,CF,\,CP21}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r},\,CF,\,CP21}",
    ),
    (
        "NNLO_fact_rdot_CF_CP03",
        r"NNLO fact, $\dot{r}$, CF, P03",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r},\,CF,\,CP03}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\dot{r},\,CF,\,CP03}",
    ),
    (
        "NNLO_fact_omg_CF_CP12",
        r"NNLO fact, $\omega$, CF, P12",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega,\,CF,\,CP12}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega,\,CF,\,CP12}",
    ),
    (
        "NNLO_fact_omg_CF_CP21",
        r"NNLO fact, $\omega$, CF, P21",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega,\,CF,\,CP21}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega,\,CF,\,CP21}",
    ),
    (
        "NNLO_fact_omg_CF_CP03",
        r"NNLO fact, $\omega$, CF, P03",
        r"\dot{E}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega,\,CF,\,CP03}",
        r"\dot{J}_\mathrm{PN}^\mathrm{NNLO,\,fact,\,\omega,\,CF,\,CP03}",
    ),
]
pn_orders_names = {
    order: (e_name, j_name) for order, _, e_name, j_name in pn_orders_combinations
}
pn_orders_labels = {order: label for order, label, _, _ in pn_orders_combinations}
