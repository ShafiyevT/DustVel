# -*- coding: utf-8 -*-
"""DustVel — buyruqlar qatori interfeysi."""
import argparse, json, sys, os
import numpy as np

from .core import Fraksiya, Muhit, Hisoblagich, tahlil_qil
from . import STANDART_FRAKSIYALAR, __version__

BAR = "=" * 74


def _sarlavha():
    print(BAR)
    print(f"  DustVel v{__version__} — chang zarrachalari harakat tezliklarini hisoblash")
    print("  Buxoro davlat universiteti")
    print(BAR)


def _fraksiyalar(args):
    if args.fraksiya_fayl:
        with open(args.fraksiya_fayl, encoding="utf-8") as f:
            data = json.load(f)
        return [Fraksiya(d["nom"], d["d"], d["rho_p"]) for d in data]
    if args.d:
        return [Fraksiya(f"d = {args.d*1e6:.1f} um", args.d, args.rho_p)]
    return STANDART_FRAKSIYALAR


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="dustvel",
        description="Chang zarrachalarining harakat tezliklari va traektoriyasini hisoblash")
    p.add_argument("-u", "--ux", type=float, default=5.0, help="shamol tezligi, m/s")
    p.add_argument("-H", "--h0", type=float, default=20.0, help="chiqarish balandligi, m")
    p.add_argument("-T", "--vaqt", type=float, default=200.0, help="hisoblash davri, s")
    p.add_argument("-s", "--dt", type=float, default=1e-3, help="vaqt qadami, s")
    p.add_argument("-d", type=float, default=None, help="zarracha diametri, m")
    p.add_argument("--rho-p", type=float, default=2500.0, help="zarracha zichligi, kg/m3")
    p.add_argument("--harorat", type=float, default=293.15, help="havo harorati, K")
    p.add_argument("--fraksiya-fayl", default=None, help="fraksiyalar JSON fayli")
    p.add_argument("-o", "--chiqish", default=None, help="natijalarni CSV ga saqlash")
    p.add_argument("--taqqoslash", action="store_true",
                   help="soddalashtirilgan model bilan taqqoslash")
    p.add_argument("-v", "--version", action="version", version=f"DustVel {__version__}")
    a = p.parse_args(argv)

    muhit = Muhit(ux=a.ux, T_muhit=a.harorat)
    h = Hisoblagich(muhit)
    frs = _fraksiyalar(a)

    _sarlavha()
    print(f"  Shamol tezligi      : {a.ux:.2f} m/s")
    print(f"  Chiqarish balandligi: {a.h0:.1f} m")
    print(f"  Havo harorati       : {a.harorat:.2f} K")
    print(f"  Qovushqoqlik        : {muhit.qovushqoqlik():.4e} Pa*s")
    print(f"  Vaqt qadami         : {a.dt:.1e} s   (makroqadam)")
    print(f"  Kichik qadam        : {a.dt/3:.3e} s")
    print(BAR)
    print(f"{'Fraksiya':<16}{'tau, ms':>10}{'v_t, m/s':>12}{'Re_p':>8}"
          f"{'t_yer, s':>11}{'x_yer, m':>11}")
    print("-" * 74)

    natijalar = []
    for fr in frs:
        n = h.hisobla(fr, T=a.vaqt, dt=a.dt, H0=a.h0)
        natijalar.append((fr, n))
        ty = f"{n.yer_t:.1f}" if n.yer_t else "> T"
        xy = f"{n.yer_x:.1f}" if n.yer_x else "-"
        print(f"{fr.nom:<16}{h.tau(fr)*1e3:>10.4f}{h.chokish_tezligi(fr):>12.4e}"
              f"{h.reynolds(fr):>8.3f}{ty:>11}{xy:>11}")
    print(BAR)

    if a.taqqoslash:
        print("\n  Soddalashtirilgan model (v_p = u) bilan taqqoslash")
        print("-" * 74)
        print(f"{'Fraksiya':<16}{'u*tau, m':>12}{'max |dx|, m':>14}{'nisbiy, %':>12}")
        print("-" * 74)
        for fr in frs:
            r = tahlil_qil(h, fr, T=min(a.vaqt, 10.0), dt=a.dt, H0=a.h0)
            print(f"{fr.nom:<16}{r['dx_asimptotik']:>12.4f}"
                  f"{r['dx_max']:>14.4f}{r['nisbiy_oxirgi']:>12.4f}")
        print(BAR)

    for fr in frs:
        if h.reynolds(fr) > 1.0:
            print(f"  [ogohlantirish] {fr.nom}: Re_p = {h.reynolds(fr):.2f} > 1, "
                  f"Stoks yaqinlashuvi chegarasiga yaqin")

    if a.chiqish:
        os.makedirs(os.path.dirname(a.chiqish) or ".", exist_ok=True)
        with open(a.chiqish, "w", encoding="utf-8") as f:
            f.write("fraksiya,t,vx,vy,vz,x,y,z\n")
            for fr, n in natijalar:
                step = max(1, len(n.t) // 2000)
                for i in range(0, len(n.t), step):
                    f.write(f"{fr.nom},{n.t[i]:.6f},{n.vx[i]:.6f},{n.vy[i]:.6f},"
                            f"{n.vz[i]:.6f},{n.x[i]:.4f},{n.y[i]:.4f},{n.z[i]:.4f}\n")
        print(f"  Natijalar saqlandi: {a.chiqish}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
