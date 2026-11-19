# -*- coding: utf-8 -*-
"""DustVel 2.0 yadrosining testlari."""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dustvel import (Fraksiya, Muhit, Hisoblagich, STANDART_FRAKSIYALAR,
                     rk4_etalon, f_Re)

MU = 1.81e-5
FR5 = Fraksiya("50-100 um", 70.711e-6, 2600.0)


def _h(qarshilik="re", ux=12.0, uy=3.2, **kw):
    return Hisoblagich(Muhit(mu_a=MU, ux=ux, uy=uy, **kw), harorat_tuzatma=False,
                       qarshilik=qarshilik)


def test_relaksatsiya_vaqti():
    h = _h()
    tau = FR5.rho_p * FR5.d ** 2 / (18 * MU)
    assert abs(h.tau(FR5) - tau) / tau < 1e-12


def test_muvozanat_chokish_tezliklari():
    """(3.6) tenglamaning yechimi: III bob, § 3.1 dagi qiymatlar."""
    h = _h()
    kutilgan = [2.26e-5, 1.504e-3, 1.315e-2, 6.945e-2, 0.3253]
    for fr, w in zip(STANDART_FRAKSIYALAR, kutilgan):
        ws = h.chokish_tezligi(fr)
        assert abs(ws - w) / w < 3e-3, (fr.nom, ws, w)
        Re = h.reynolds(fr)
        assert abs(ws * f_Re(Re) - 9.81 * h.tau(fr)) < 1e-12


def test_stoks_rejimi():
    """qarshilik='stoks' da w_s = g tau va 50-100 um uchun farq ~17 %."""
    hs, hr = _h("stoks"), _h("re")
    assert abs(hs.chokish_tezligi(FR5) - 9.81 * hs.tau(FR5)) < 1e-15
    farq = (hs.chokish_tezligi(FR5) - hr.chokish_tezligi(FR5)) / hs.chokish_tezligi(FR5)
    assert 0.16 < farq < 0.18


def test_yaqinlashish_tartibi_chiziqsiz():
    """f(Re) bilan: tezlik xatoligi O(dt^2) tartibda kamayadi."""
    h = _h(); T = 0.12
    ref = rk4_etalon(h, FR5, T, qadam=2e-6)
    xat = []
    for dt in (0.03, 0.015, 0.0075, 0.00375):
        n = h.hisobla(FR5, T=T, dt=dt, H0=1e9, yerda_toxta=False)
        xat.append(np.max(np.abs(np.array([n.vx[-1], n.vy[-1], n.vz[-1]]) - ref)))
    p = [np.log2(xat[k] / xat[k + 1]) for k in range(len(xat) - 1)]
    assert 1.9 < p[-1] < 2.2, p
    assert 1.9 < p[-2] < 2.6, p


def test_yaqinlashish_tartibi_stoks():
    """Stoks rejimida x(T) analitik yechim bilan solishtiriladi."""
    h = _h("stoks", ux=5.0, uy=0.0); fr = Fraksiya("50 um", 50e-6, 2500.0)
    tau = h.tau(fr); T = 1.0
    x_aniq = 5.0 * (T - tau * (1 - np.exp(-T / tau)))
    xat = []
    for dt in (1e-2, 1e-3, 1e-4):
        n = h.hisobla(fr, T=T, dt=dt, H0=1e9, yerda_toxta=False)
        xat.append(abs(n.x[-1] - x_aniq))
    assert 1.9 < np.log10(xat[0] / xat[1]) < 2.1
    assert 1.9 < np.log10(xat[1] / xat[2]) < 2.1


def test_iteratsiya_yaqinlashishi():
    """Katta qadamda ham (3.26) iteratsiya yaqinlashadi (q < 0.687)."""
    h = _h()
    n = h.hisobla(FR5, T=30.0, dt=10.0, H0=1e9, yerda_toxta=False,
                  boshlangich="impuls")
    assert n.iter_maks < 40
    assert np.all(np.isfinite(n.vz))


def test_turgunlik_katta_qadamda():
    """L-turg'unlik: katta qadamda tezlik muvozanat qiymatiga intiladi."""
    h = _h()
    n = h.hisobla(FR5, T=100.0, dt=1.0, H0=1e9, yerda_toxta=False)
    assert abs(n.vx[-1] - 12.0) < 1e-6 and abs(n.vy[-1] - 3.2) < 1e-6
    assert abs(n.vz[-1] + h.chokish_tezligi(FR5)) < 1e-6


def test_impuls_sharti():
    h = _h()
    n = h.hisobla(FR5, T=0.01, dt=1e-3, H0=20.0, boshlangich="impuls")
    assert abs(n.vx[0] - 12.0) < 1e-12 and abs(n.vy[0] - 3.2) < 1e-12
    assert abs(n.vz[0]) < 1e-12


def test_log_profil():
    """Logarifmik profilda zarracha pastga tushgan sari sekinlashadi."""
    h = _h(ux=8.0, uy=0.0, profil="log")
    n = h.hisobla(FR5, T=200.0, dt=1e-2, H0=20.0)
    assert n.yer_t is not None and n.yer_x > 0
    assert n.vx[-1] < n.vx[len(n.vx) // 4]


if __name__ == "__main__":
    for nm, f in sorted(globals().items()):
        if nm.startswith("test_"):
            f(); print(f"  OK  {nm}")
    print("Barcha testlar muvaffaqiyatli o'tdi.")
