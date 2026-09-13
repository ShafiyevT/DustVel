# -*- coding: utf-8 -*-
"""DustVel yadrosining testlari."""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dustvel import Fraksiya, Muhit, Hisoblagich

MU, RHO, UX = 1.81e-5, 2500.0, 5.0
FR = Fraksiya("test 50 um", 50e-6, RHO)


def _h():
    return Hisoblagich(Muhit(mu_a=MU, ux=UX), harorat_tuzatma=False)


def test_relaksatsiya_vaqti():
    h = _h()
    tau_kutilgan = RHO * FR.d ** 2 / (18 * MU)
    assert abs(h.tau(FR) - tau_kutilgan) / tau_kutilgan < 1e-12


def test_chokish_tezligi():
    h = _h()
    assert abs(h.chokish_tezligi(FR) - 0.1882) < 1e-3


def test_yaqinlashish_tartibi():
    """x(T) bo'yicha xatolik O(dt^2) tartibda kamayishi kerak."""
    h = _h(); tau = h.tau(FR); T = 1.0
    x_aniq = UX * (T - tau * (1 - np.exp(-T / tau)))
    xat = []
    for dt in (1e-2, 1e-3, 1e-4):
        n = h.hisobla(FR, T=T, dt=dt, H0=1e9, yerda_toxta=False)
        xat.append(abs(n.x[-1] - x_aniq))
    p1 = np.log10(xat[0] / xat[1])
    p2 = np.log10(xat[1] / xat[2])
    assert 1.9 < p1 < 2.1, f"p1={p1}"
    assert 1.9 < p2 < 2.1, f"p2={p2}"


def test_turgunlik_katta_qadamda():
    """Katta vaqt qadamida ham yechim chegaralangan bo'lishi kerak."""
    h = _h()
    n = h.hisobla(FR, T=100.0, dt=1.0, H0=1e9, yerda_toxta=False)
    assert np.all(np.isfinite(n.vx))
    assert abs(n.vx[-1] - UX) < 1e-6
    assert abs(n.vz[-1] + h.chokish_tezligi(FR)) < 1e-3


def test_impuls_sharti():
    h = _h()
    n = h.hisobla(FR, T=0.01, dt=1e-4, H0=20.0, boshlangich="impuls")
    assert abs(n.vx[0] - UX) < 1e-12
    assert abs(n.vz[0]) < 1e-12


if __name__ == "__main__":
    for nm, f in sorted(globals().items()):
        if nm.startswith("test_"):
            f(); print(f"  OK  {nm}")
    print("Barcha testlar muvaffaqiyatli o'tdi.")
