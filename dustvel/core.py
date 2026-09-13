# -*- coding: utf-8 -*-
"""DustVel hisoblash yadrosi.

Uch bosqichli sxema:
    n        -> n + 1/3 : Krank-Nikolson   (3.25), (3.37)
    n + 1/3  -> n + 2/3 : BDF2             (3.29), (3.39)
    n + 2/3  -> n + 1   : BDF2             (3.31), (3.41)
Traektoriya (3.42)-(3.44) trapetsiya formulasi bilan integrallanadi.
"""
from dataclasses import dataclass
import numpy as np

G = 9.81


@dataclass
class Fraksiya:
    """Chang zarrachalari fraksiyasi."""
    nom: str
    d: float          # diametr, m
    rho_p: float      # zichlik, kg/m3

    @property
    def r(self) -> float:
        return 0.5 * self.d


@dataclass
class Muhit:
    """Atmosfera parametrlari."""
    mu_a: float = 1.81e-5
    rho_a: float = 1.225
    ux: float = 5.0
    uy: float = 0.0
    T_muhit: float = 293.15

    def qovushqoqlik(self) -> float:
        """Sazerlend formulasi bo'yicha haroratga tuzatma."""
        S, T0, mu0 = 110.4, 273.15, 1.716e-5
        T = self.T_muhit
        return mu0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)


@dataclass
class Natija:
    t: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    tau: float
    vt: float
    fraksiya: Fraksiya = None
    yer_t: float = None
    yer_x: float = None


class Hisoblagich:
    """Uch bosqichli CN + BDF2 + BDF2 algoritmi."""

    def __init__(self, muhit: Muhit, harorat_tuzatma: bool = True):
        self.muhit = muhit
        self.harorat_tuzatma = harorat_tuzatma

    def beta(self, fr: Fraksiya) -> float:
        mu = self.muhit.qovushqoqlik() if self.harorat_tuzatma else self.muhit.mu_a
        return 9.0 * mu / (2.0 * fr.rho_p * fr.r ** 2)

    def tau(self, fr: Fraksiya) -> float:
        return 1.0 / self.beta(fr)

    def chokish_tezligi(self, fr: Fraksiya) -> float:
        return G * self.tau(fr)

    def reynolds(self, fr: Fraksiya) -> float:
        mu = self.muhit.qovushqoqlik() if self.harorat_tuzatma else self.muhit.mu_a
        return self.muhit.rho_a * self.chokish_tezligi(fr) * fr.d / mu

    @staticmethod
    def _gorizontal(v0, b, dt, u):
        bd = b * dt
        v13 = ((6.0 - bd) * v0 + 2.0 * bd * u) / (6.0 + bd)
        v23 = (12.0 * v13 - 3.0 * v0 + 2.0 * bd * u) / (9.0 + 2.0 * bd)
        v1 = (12.0 * v23 - 3.0 * v13 + 2.0 * bd * u) / (9.0 + 2.0 * bd)
        return v13, v23, v1

    @staticmethod
    def _vertikal(v0, b, dt):
        bd = b * dt
        v13 = ((6.0 - bd) * v0 - 2.0 * G * dt) / (6.0 + bd)
        v23 = (12.0 * v13 - 3.0 * v0 - 2.0 * dt * G) / (9.0 + 2.0 * bd)
        v1 = (12.0 * v23 - 3.0 * v13 - 2.0 * dt * G) / (9.0 + 2.0 * bd)
        return v13, v23, v1

    def hisobla(self, fr: Fraksiya, T: float = 10.0, dt: float = 1e-3,
                H0: float = 20.0, yerda_toxta: bool = True,
                boshlangich: str = "harakatsiz") -> Natija:
        """boshlangich: 'harakatsiz' -> v(0)=0 (3.4);  'impuls' -> v(0)=u (3.45)."""
        b = self.beta(fr)
        n = int(round(T / dt))
        dtau = dt / 3.0
        u_x, u_y = self.muhit.ux, self.muhit.uy

        t = np.zeros(n + 1)
        vx = np.zeros(n + 1); vy = np.zeros(n + 1); vz = np.zeros(n + 1)
        X = np.zeros(n + 1); Y = np.zeros(n + 1); Z = np.zeros(n + 1)

        if boshlangich == "impuls":
            vx[0], vy[0], vz[0] = u_x, u_y, 0.0
        else:
            vx[0], vy[0], vz[0] = 0.0, 0.0, 0.0
        Z[0] = H0

        for k in range(n):
            ax13, ax23, ax1 = self._gorizontal(vx[k], b, dt, u_x)
            ay13, ay23, ay1 = self._gorizontal(vy[k], b, dt, u_y)
            az13, az23, az1 = self._vertikal(vz[k], b, dt)

            X[k + 1] = X[k] + dtau * (0.5 * vx[k] + ax13 + ax23 + 0.5 * ax1)
            Y[k + 1] = Y[k] + dtau * (0.5 * vy[k] + ay13 + ay23 + 0.5 * ay1)
            Z[k + 1] = Z[k] + dtau * (0.5 * vz[k] + az13 + az23 + 0.5 * az1)

            vx[k + 1], vy[k + 1], vz[k + 1] = ax1, ay1, az1
            t[k + 1] = t[k] + dt

            if yerda_toxta and Z[k + 1] <= 0.0:
                w = Z[k] / (Z[k] - Z[k + 1])
                yer_t = t[k] + w * dt
                yer_x = X[k] + w * (X[k + 1] - X[k])
                Z[k + 1] = 0.0
                sl = slice(0, k + 2)
                return Natija(t[sl], vx[sl], vy[sl], vz[sl], X[sl], Y[sl], Z[sl],
                              self.tau(fr), self.chokish_tezligi(fr), fr, yer_t, yer_x)

        return Natija(t, vx, vy, vz, X, Y, Z,
                      self.tau(fr), self.chokish_tezligi(fr), fr, None, None)

    def soddalashtirilgan(self, fr: Fraksiya, T=10.0, dt=1e-3, H0=20.0):
        """v_p = u soddalashtirishi."""
        vt = self.chokish_tezligi(fr)
        n = int(round(T / dt))
        t = np.linspace(0, T, n + 1)
        X = self.muhit.ux * t
        Z = np.maximum(H0 - vt * t, 0.0)
        return t, X, Z


def aniq_yechim(h: Hisoblagich, fr: Fraksiya, t):
    """Doimiy shamolda analitik yechim (3.50), (3.51)."""
    tau = h.tau(fr)
    u = h.muhit.ux
    vx = u * (1.0 - np.exp(-t / tau)) + u * np.exp(-t / tau)  # v(0)=u
    x = u * t
    return vx, x


def tahlil_qil(h: Hisoblagich, fr: Fraksiya, T=10.0, dt=1e-3, H0=20.0):
    """Ishlab chiqilgan va soddalashtirilgan modellarni taqqoslaydi."""
    n = h.hisobla(fr, T=T, dt=dt, H0=H0, yerda_toxta=False)
    t2, X2, Z2 = h.soddalashtirilgan(fr, T=T, dt=dt, H0=H0)
    m = min(len(n.t), len(t2))
    dx = X2[:m] - n.x[:m]
    nisbiy = np.divide(np.abs(dx), np.maximum(n.x[:m], 1e-12)) * 100.0
    return {
        "tau_ms": h.tau(fr) * 1e3,
        "vt": h.chokish_tezligi(fr),
        "Re_p": h.reynolds(fr),
        "dx_asimptotik": h.muhit.ux * h.tau(fr),
        "dx_max": float(np.max(np.abs(dx))),
        "nisbiy_oxirgi": float(nisbiy[-1]),
        "t": n.t[:m],
        "nisbiy": nisbiy,
    }
