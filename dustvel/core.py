# -*- coding: utf-8 -*-
"""DustVel hisoblash yadrosi (2.0).

Zarracha harakati tenglamalari (dissertatsiya, III bob, (3.2)-(3.5), (3.13)):

    dv_xi/dt = lambda(v) (u_xi - v_xi) - g delta_{xi z},
    lambda(v) = f(Re)/tau,  tau = rho_p d^2 / (18 mu_a),
    Re = rho_a d |u - v| / mu_a,  f(Re) = 1 + 0.15 Re^0.687.

Makroqadam Delta_tau uchta kichik qadamga bo'linadi (3.16):

    n       -> n + 1/3 : simmetrik oshkormas sxema (Krank-Nikolson)            (3.20)
    n + 1/3 -> n + 2/3 : ikkinchi tartibli uch qatlamli oshkormas sxema (BDF2)  (3.23)
    n + 2/3 -> n + 1   : ikkinchi tartibli uch qatlamli oshkormas sxema (BDF2)  (3.25)

Har bir kichik qadamda lambda(v) ga nisbatan chiziqsiz tenglamalar
komponentalar bo'yicha ketma-ket iteratsiya (3.26) bilan yechiladi va
(3.27) shart bajarilganda to'xtatiladi. Iteratsiyaning qisqartirish
koeffitsiyenti q < 0.687 (3.29), shuning uchun u istalgan Delta_tau da
yaqinlashadi.

Muvozanat cho'kish tezligi (3.6) oddiy iteratsiya (3.30) bilan topiladi.
Traektoriya har bir kichik qadamda trapetsiya formulasi bilan integrallanadi.

qarshilik="stoks" tanlansa f(Re) = 1 bo'ladi va algoritm II bobdagi Stoks
yaqinlashuviga (2.13) o'tadi (taqqoslash uchun saqlangan).
"""
from dataclasses import dataclass, field
import math
import numpy as np

G = 9.81
KARMAN = 0.4


def f_Re(Re):
    """Reynolds soniga bog'liq qarshilik tuzatmasi (3.4)."""
    return 1.0 + 0.15 * np.power(np.maximum(Re, 0.0), 0.687)


@dataclass
class Fraksiya:
    """Chang zarrachalari fraksiyasi."""
    nom: str
    d: float              # o'rta geometrik diametr, m
    rho_p: float          # zichlik, kg/m3
    p: float = None       # emissiyadagi massa ulushi (2.25), ixtiyoriy

    @property
    def r(self) -> float:
        return 0.5 * self.d


@dataclass
class Muhit:
    """Atmosfera parametrlari.

    profil = "doimiy" : shamol tezligi (ux, uy, uz) balandlik bo'yicha o'zgarmas;
    profil = "log"    : (ux, uy) — 10 m dagi shamol, gorizontal tezlik (2.23)
                        logarifmik profil bo'yicha, h_s dan yuqorida o'zgarmas.
    """
    mu_a: float = 1.81e-5
    rho_a: float = 1.225
    ux: float = 5.0
    uy: float = 0.0
    uz: float = 0.0
    T_muhit: float = 293.15
    profil: str = "doimiy"
    z0: float = 1.0e-4
    h_s: float = 100.0

    def qovushqoqlik(self) -> float:
        """Sazerlend formulasi bo'yicha haroratga tuzatma."""
        S, T0, mu0 = 110.4, 273.15, 1.716e-5
        T = self.T_muhit
        return mu0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)

    @property
    def u_star(self) -> float:
        U10 = np.hypot(self.ux, self.uy)
        return KARMAN * U10 / np.log(10.0 / self.z0)

    def shamol(self, z):
        """z balandlikdagi shamol tezligi vektori (u, v, w)."""
        if self.profil != "log":
            return (self.ux, self.uy, self.uz)
        U10 = math.hypot(self.ux, self.uy)
        if U10 == 0.0:
            return (0.0, 0.0, self.uz)
        zz = min(max(z, 0.0), self.h_s)
        U = self.u_star / KARMAN * math.log((zz + self.z0) / self.z0)
        return (U * self.ux / U10, U * self.uy / U10, self.uz)


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
    vt: float                       # muvozanat cho'kish tezligi w_s
    fraksiya: Fraksiya = None
    yer_t: float = None
    yer_x: float = None
    Re: np.ndarray = None           # har bir makroqadamdagi Reynolds soni
    iteratsiyalar: list = field(default_factory=list)

    @property
    def iter_ortacha(self) -> float:
        return float(np.mean(self.iteratsiyalar)) if self.iteratsiyalar else 0.0

    @property
    def iter_maks(self) -> int:
        return int(max(self.iteratsiyalar)) if self.iteratsiyalar else 0


class Hisoblagich:
    """Uch bosqichli (simmetrik + 2 x uch qatlamli) oshkormas algoritm."""

    def __init__(self, muhit: Muhit, harorat_tuzatma: bool = True,
                 qarshilik: str = "re", eps: float = 1e-10, smax: int = 50):
        if qarshilik not in ("re", "stoks"):
            raise ValueError("qarshilik: 're' yoki 'stoks'")
        self.muhit = muhit
        self.harorat_tuzatma = harorat_tuzatma
        self.qarshilik = qarshilik
        self.eps = eps
        self.smax = smax

    # ------------------------------------------------------------------
    # fizik parametrlar
    # ------------------------------------------------------------------
    def mu(self) -> float:
        return self.muhit.qovushqoqlik() if self.harorat_tuzatma else self.muhit.mu_a

    def tau(self, fr: Fraksiya) -> float:
        """Relaksatsiya vaqti tau = rho_p d^2/(18 mu_a)  (3.4)."""
        return fr.rho_p * fr.d ** 2 / (18.0 * self.mu())

    def beta(self, fr: Fraksiya) -> float:
        """Stoks qarshiligi koeffitsiyenti 1/tau."""
        return 1.0 / self.tau(fr)

    def _kRe(self, fr: Fraksiya) -> float:
        return self.muhit.rho_a * fr.d / self.mu()

    def f(self, fr: Fraksiya, nisbiy_tezlik):
        if self.qarshilik == "stoks":
            return np.ones_like(np.asarray(nisbiy_tezlik, dtype=float))
        return f_Re(self._kRe(fr) * np.asarray(nisbiy_tezlik, dtype=float))

    def lam(self, fr: Fraksiya, u, v) -> float:
        """lambda(v) = f(Re)/tau  (3.13)."""
        return float(self.f(fr, np.linalg.norm(np.asarray(u) - np.asarray(v)))) / self.tau(fr)

    def stoks_tezligi(self, fr: Fraksiya) -> float:
        """Stoks qonuni bo'yicha cho'kish tezligi w_g = g tau  (2.13)."""
        return G * self.tau(fr)

    def muvozanat_tezligi(self, fr: Fraksiya, eps: float = 1e-12, smax: int = 200):
        """Muvozanat cho'kish tezligi w_s: w_s f(Re(w_s)) = g tau  (3.6), (3.30).

        Qaytaradi: (w_s, iteratsiyalar soni).
        """
        gt = self.stoks_tezligi(fr)
        if self.qarshilik == "stoks":
            return gt, 0
        k = self._kRe(fr)
        w = gt
        for m in range(1, smax + 1):
            wn = gt / f_Re(k * w)
            if abs(wn - w) <= eps * max(wn, 1e-30):
                return float(wn), m
            w = wn
        return float(w), smax

    def chokish_tezligi(self, fr: Fraksiya) -> float:
        return self.muvozanat_tezligi(fr)[0]

    def reynolds(self, fr: Fraksiya) -> float:
        """Muvozanat cho'kish holatidagi Reynolds soni."""
        return self._kRe(fr) * self.chokish_tezligi(fr)

    # ------------------------------------------------------------------
    # kichik qadamlar (skalyar arifmetika — tezlik uchun)
    # ------------------------------------------------------------------
    def _lam_s(self, kRe, tau, u, v):
        if self.qarshilik == "stoks":
            return 1.0 / tau
        r = math.sqrt((u[0] - v[0]) ** 2 + (u[1] - v[1]) ** 2 + (u[2] - v[2]) ** 2)
        return (1.0 + 0.15 * (kRe * r) ** 0.687) / tau

    def _kichik_qadam(self, fr, dt, x0, v0, vm, sxema):
        """Bitta kichik qadam va (3.26) iteratsiya.

        sxema = "kn"  : (3.20), v0 = v^n;
        sxema = "bdf2": (3.23)/(3.25), v0 = v^k, vm = v^{k-1}.
        """
        h = dt / 3.0
        tau = self.tau(fr); kRe = self._kRe(fr)
        sham = self.muhit.shamol
        dzc = (0.0, 0.0, 1.0)
        if sxema == "kn":
            u0 = sham(x0[2])
            lam0 = self._lam_s(kRe, tau, u0, v0)
            v = list(v0)
        else:
            v = [2.0 * v0[c] - vm[c] for c in range(3)]
        eps, smax = self.eps, self.smax
        for m in range(1, smax + 1):
            dmax = 0.0
            for c in range(3):
                z_new = x0[2] + 0.5 * h * (v0[2] + v[2])
                u1 = sham(z_new)
                lam = self._lam_s(kRe, tau, u1, v)
                if sxema == "kn":
                    vn = ((6.0 - dt * lam0) * v0[c] + dt * (lam * u1[c] + lam0 * u0[c])
                          - 2.0 * dt * G * dzc[c]) / (6.0 + dt * lam)
                else:
                    vn = (12.0 * v0[c] - 3.0 * vm[c] + 2.0 * dt * lam * u1[c]
                          - 2.0 * dt * G * dzc[c]) / (9.0 + 2.0 * dt * lam)
                dmax = max(dmax, abs(vn - v[c]))
                v[c] = vn
            if dmax <= eps * max(1.0, math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)):
                break
        x = [x0[c] + 0.5 * h * (v0[c] + v[c]) for c in range(3)]
        return x, v, m

    def _kn(self, fr, dt, x0, v0):
        """n -> n+1/3: simmetrik oshkormas sxema (Krank-Nikolson)  (3.20)."""
        return self._kichik_qadam(fr, dt, x0, v0, None, "kn")

    def _bdf2(self, fr, dt, x1, v1, v0):
        """Ikkinchi tartibli uch qatlamli oshkormas sxema (BDF2)  (3.23), (3.25)."""
        return self._kichik_qadam(fr, dt, x1, v1, v0, "bdf2")

    def makroqadam(self, fr, dt, x, v):
        """Bitta makroqadam: (x^n, v^n) -> (x^{n+1}, v^{n+1})."""
        x13, v13, m1 = self._kn(fr, dt, x, v)
        x23, v23, m2 = self._bdf2(fr, dt, x13, v13, v)
        x1, v1, m3 = self._bdf2(fr, dt, x23, v23, v13)
        return x1, v1, (m1, m2, m3)

    # ------------------------------------------------------------------
    def hisobla(self, fr: Fraksiya, T: float = 10.0, dt: float = 1e-3,
                H0: float = 20.0, yerda_toxta: bool = True,
                boshlangich: str = "harakatsiz") -> Natija:
        """Tezliklar va traektoriyani hisoblaydi.

        boshlangich: 'harakatsiz' -> v(0) = 0;
                     'impuls'     -> v_x(0) = u, v_y(0) = v, v_z(0) = 0  (3.5).
        """
        n = int(round(T / dt))
        t = np.zeros(n + 1)
        V = np.zeros((n + 1, 3)); X = np.zeros((n + 1, 3)); Re = np.zeros(n + 1)
        X[0] = (0.0, 0.0, H0)
        if boshlangich == "impuls":
            u = self.muhit.shamol(H0)
            V[0] = (u[0], u[1], 0.0)
        its = []
        kRe = self._kRe(fr)
        Re[0] = kRe * np.linalg.norm(self.muhit.shamol(H0) - V[0])
        ws = self.chokish_tezligi(fr)

        xk, vk = list(X[0]), list(V[0])
        for k in range(n):
            xk, vk, m = self.makroqadam(fr, dt, xk, vk)
            X[k + 1] = xk; V[k + 1] = vk
            its.extend(m)
            t[k + 1] = t[k] + dt
            Re[k + 1] = kRe * np.linalg.norm(self.muhit.shamol(X[k + 1, 2]) - V[k + 1])
            if yerda_toxta and X[k + 1, 2] <= 0.0:
                w = X[k, 2] / (X[k, 2] - X[k + 1, 2])
                yer_t = t[k] + w * dt
                yer_x = X[k, 0] + w * (X[k + 1, 0] - X[k, 0])
                X[k + 1, 2] = 0.0
                s = slice(0, k + 2)
                return Natija(t[s], V[s, 0], V[s, 1], V[s, 2], X[s, 0], X[s, 1], X[s, 2],
                              self.tau(fr), ws, fr, yer_t, yer_x, Re[s], its)

        return Natija(t, V[:, 0], V[:, 1], V[:, 2], X[:, 0], X[:, 1], X[:, 2],
                      self.tau(fr), ws, fr, None, None, Re, its)

    def soddalashtirilgan(self, fr: Fraksiya, T=10.0, dt=1e-3, H0=20.0):
        """v_p = u soddalashtirishi (II bobdagi chiziqli model), cho'kish w_s bilan."""
        ws = self.chokish_tezligi(fr)
        n = int(round(T / dt))
        t = np.linspace(0, T, n + 1)
        X = np.zeros_like(t); Z = np.zeros_like(t)
        X[0], Z[0] = 0.0, H0
        for k in range(n):
            u = self.muhit.shamol(Z[k])
            X[k + 1] = X[k] + u[0] * dt
            Z[k + 1] = max(Z[k] - ws * dt, 0.0)
        return t, X, Z


def aniq_yechim_stoks(h: Hisoblagich, fr: Fraksiya, t, boshlangich="harakatsiz"):
    """Doimiy shamol va Stoks qarshiligida gorizontal harakatning aniq yechimi."""
    tau = h.tau(fr)
    u = h.muhit.ux
    if boshlangich == "impuls":
        return np.full_like(t, u), u * t
    vx = u * (1.0 - np.exp(-t / tau))
    x = u * (t - tau * (1.0 - np.exp(-t / tau)))
    return vx, x


def rk4_etalon(h: Hisoblagich, fr: Fraksiya, T: float, v0=(0.0, 0.0, 0.0), qadam=None):
    """Doimiy shamolda (3.13) tizimning juda kichik qadamli RK4 yechimi (etalon)."""
    tau = h.tau(fr)
    qadam = qadam or tau / 200.0
    U = h.muhit.shamol(0.0)

    def F(v):
        r = h.lam(fr, U, v) * (U - v)
        r[2] -= G
        return r

    v = np.array(v0, dtype=float)
    n = int(np.ceil(T / qadam)); hh = T / n
    for _ in range(n):
        k1 = F(v); k2 = F(v + hh / 2 * k1); k3 = F(v + hh / 2 * k2); k4 = F(v + hh * k3)
        v = v + hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return v


def tahlil_qil(h: Hisoblagich, fr: Fraksiya, T=10.0, dt=1e-3, H0=20.0):
    """Ishlab chiqilgan va soddalashtirilgan (v_p = u) modellarni taqqoslaydi."""
    n = h.hisobla(fr, T=T, dt=dt, H0=H0, yerda_toxta=False)
    t2, X2, Z2 = h.soddalashtirilgan(fr, T=T, dt=dt, H0=H0)
    m = min(len(n.t), len(t2))
    dx = X2[:m] - n.x[:m]
    nisbiy = np.divide(np.abs(dx), np.maximum(n.x[:m], 1e-12)) * 100.0
    return {
        "tau_ms": h.tau(fr) * 1e3,
        "w_g": h.stoks_tezligi(fr),
        "w_s": h.chokish_tezligi(fr),
        "Re_p": h.reynolds(fr),
        "dx_asimptotik": h.muhit.ux * h.tau(fr),
        "dx_max": float(np.max(np.abs(dx))),
        "nisbiy_oxirgi": float(nisbiy[-1]),
        "t": n.t[:m],
        "nisbiy": nisbiy,
    }
