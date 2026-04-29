"""
Load and plot the teukode output data, including particle trajectory and horizon fluxes.
Rossella Gamba and Danilo Chiaramello, 2026
"""

import numpy as np
import os, glob, scipy as sp
import matplotlib.pyplot as plt
import utils as ut


class PlotterTKO(object):
    """
    Class to plot the teukode output data
    """

    def __init__(
        self,
        a,
        path_dyn,  # Path to dynamics files
        path_tk,  # Path to output folders
        mu=1.0e-3,
    ) -> None:

        self.a = a
        self.mu = mu
        self.path_dyn = path_dyn
        self.path_tk = path_tk
        self.hflx = None
        self.dyn = None

        self.__read_sims__()

        pass

    def __read_sims__(self):
        """
        Find available simulations (resolution, cfl, m values) in tk folder
        """
        self.sims = {}
        globstr = os.path.join(self.path_tk, f"teuk_a{self.a:.4f}_r0*_*x*_m*/")
        sims = glob.glob(globstr)
        if len(sims) == 0:
            raise FileNotFoundError("PlotterTKO.__read_sims__(): no simulations found.")
        nx_max = 0
        for sim in sims:
            sub = sim[len(self.path_tk) :]
            nx, ny = (eval(sub.split("x")[kk].split("_")[kk - 1]) for kk in range(2))
            cfl = 2.0 if "cfl" not in sub else eval(sub.split("cfl")[1].split("/")[0])
            m = eval(sub.split("_m")[1].split("_")[0].split("/")[0])
            if (nx, ny, cfl) not in self.sims.keys():
                self.sims[(nx, ny, cfl)] = []
            self.sims[(nx, ny, cfl)].append(m)
            if nx > nx_max:
                nx_max = nx
        self.sims[(nx, ny, cfl)].sort()
        ny_max = max([sim[1] for sim in self.sims.keys() if sim[0] == nx_max])
        self.max_res = (nx_max, ny_max)

    def load_dyn(self):
        """
        Load the particle dynamics
        """
        self.dyn = {}
        for var in ["r", "ph", "pr", "pph", "H", "prstar", "Omg"]:
            dynf = os.path.join(self.path_dyn, f"{var}.txt")
            try:
                t, self.dyn[var] = np.loadtxt(dynf, unpack=True)
            except ValueError:
                print(f"Error loading {var}.txt")
            if "t" not in self.dyn.keys():
                self.dyn["t"] = t

    def plot_trajectory(self, kind="r", color="k", lb=None, ax=None):
        """
        Plot the trajectory
        """

        if "traj" not in kind:
            dynf = os.path.join(self.path_dyn, f"{kind}.txt")
            x, y = np.loadtxt(dynf, unpack=True)
        else:
            dyn_r = os.path.join(self.path_dyn, "r.txt")
            dyn_p = os.path.join(self.path_dyn, "ph.txt")
            _, r = np.loadtxt(dyn_r, unpack=True)
            _, p = np.loadtxt(dyn_p, unpack=True)
            x = r * np.cos(p)
            y = r * np.sin(p)

        if ax == None:
            _, ax = plt.subplots()

        ax.plot(x, y, color=color, label=lb)
        return ax, [x, y]

    def integrate_hflx(self, m=None, nx=None, ny=None, cfl=2.0):
        """
        Integrate dJ/dt and dm/dt.
        Resolution taken as highest available if not specified.
        Summing over all specified values, or all available if none in input.
        Skipping negative ms if the opposite is available.
        """
        if nx is None:
            nx, ny = self.max_res
            print(
                f"\x1b[1mPlotterTKO.integrate_hflx(): \x1b[0massuming highest resolution, (nx, ny) = ({nx}, {ny})"
            )
        if m is None:
            emm = self.sims[(nx, ny, cfl)]
        elif isinstance(m, list):
            emm = m
        else:
            emm = [m]
        cfls = "" if cfl == 2.0 else f"_cfl{cfl:.2g}"

        res = {}
        for mv in emm:
            if mv < 0 and -mv in emm:
                continue
            data_dir = os.path.join(
                self.path_tk,
                f"teuk_a{self.a:.4f}_r0{self.dyn['r'][0]:.3f}_{nx}x{ny}_proc4x2_m{mv}{cfls}/out1d",
            )

            # Multiply modes other than m = 0 by 2 to account for negative m values
            by2 = 1.0 if mv == 0 else 2.0

            for kind in ["E", "J"]:
                fname = f"d{kind}dt_hrz_td_lsum_m{mv}_Poisson.dat"
                if not os.path.exists(os.path.join(data_dir, fname)):
                    raise FileNotFoundError(
                        f"\x1b[31;1mPlotterTKO.integrate_hflx(): \x1b[22mNo {kind} flux file found for m = {mv}, (nx, ny, cfl) = ({nx}, {ny}, {cfl}).\x1b[0m"
                    )
                    continue
                t, dy = np.loadtxt(os.path.join(data_dir, fname), unpack=True)
                if "t" not in res.keys():
                    res["t"] = t
                else:
                    dy = np.interp(res["t"], t, dy)
                    t = res["t"]
                y = sp.integrate.cumulative_trapezoid(dy, t, initial=0)
                if mv not in res.keys():
                    res[mv] = {}
                res[mv][kind] = {"dy": dy, "y": y}
                if "sum" not in res.keys():
                    res["sum"] = {}
                if kind not in res["sum"].keys():
                    res["sum"][kind] = {
                        key: res[mv][kind][key] * by2 for key in res[mv][kind].keys()
                    }
                else:
                    res["sum"][kind]["dy"] += res[mv][kind]["dy"] * by2
                    res["sum"][kind]["y"] += res[mv][kind]["y"] * by2

        chi = self.a + res["sum"]["J"]["y"] / (1 + res["sum"]["E"]["y"]) ** 2
        res["chi"] = chi

        self.hflx = res
        return res

    def plot_hflx(
        self,
        kind="E",
        integrate=False,
        m=None,
        nx=None,
        ny=None,
        cfl=2.0,
        ax=None,
        color="k",
        lb=None,
        dx=0.0,
    ):
        """
        Plot the horizon flux data.
        Time axis is advanced time, coinciding with particle BL time, optionally further shifted by dx.
        Max available resolution used if not otherwise specified.
        Plotting sum of all m values if not otherwise specified.

        Returns axis of plot, time, plotted flux.
        """

        if self.dyn == None:
            self.load_dyn()
        if self.hflx == None:
            self.integrate_hflx(m=m, nx=nx, ny=ny, cfl=cfl)
        if nx is None:
            nx, ny = self.max_res
        if m is None:
            emm = "sum"
            print(
                f"PlotterTKO.plot_hflx(): no m value given; summing over m = {self.sims[(nx, ny, cfl)]}"
            )
        else:
            emm = m

        x = self.hflx["t"]
        if integrate:
            y = self.hflx["sum"][kind]["y"]
        else:
            y = self.hflx["sum"][kind]["dy"]

        v0 = ut.adv_time_in(self.a, self.dyn["r"][0])
        delta_x = v0 + dx

        if ax == None:
            _, ax = plt.subplots()

        ax.plot(x + delta_x, y, color=color, label=lb)
        return ax, [x + delta_x, y]
