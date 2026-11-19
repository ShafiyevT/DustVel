# -*- coding: utf-8 -*-
"""DustVel — buyruqlar qatori interfeysi."""
import argparse, json, sys, os
import numpy as np

from .core import Fraksiya, Muhit, Hisoblagich, tahlil_qil
from . import STANDART_FRAKSIYALAR, __version__

BAR = "=" * 86


def _sarlavha():
    print(BAR)
    print(f"  DustVel v{__version__} — chang zarrachalari harakat tezliklarini hisoblash")
    print("  Buxoro davlat universiteti")
    print(BAR)


def _fraksiyalar(args):
    if args.fraksiya_fayl:
        with open(args.fraksiya_fayl, encoding="utf-8") as f:
            data = json.load(f)
        return [Fraksiya(d["nom"], d["d"], d["rho_p"], d.get("p")) for d in data]
    if args.d:
        return [Fraksiya(f"d = {args.d*1e6:.1f} um", args.d, args.rho_p)]
    return STANDART_FRAKSIYALAR


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="dustvel",
        description="Chang zarrachalarining harakat tezliklari va traektoriyasini hisoblash")
    p.add_argument("-u", "--ux", type=float, default=5.0, help="shamol tezligi (x), m/s")
    p.add_argument("--uy", type=float, default=0.0, help="shamol tezligi (y), m/s")
    p.add_argument("--profil", choices=["doimiy", "log"], default="doimiy",
                   help="shamol profili: doimiy yoki logarifmik (2.23)")
    p.add_argument("--z0", type=float, default=1e-4, help="g'adir-budurlik uzunligi, m")
    p.add_argument("--hs", type=float, default=100.0, help="yer usti qatlami balandligi, m")
    p.add_argument("-H", "--h0", type=float, default=20.0, help="chiqarish balandligi, m")
    p.add_argument("-T", "--vaqt", type=float, default=300.0, help="hisoblash davri, s")
    p.add_argument("-s", "--dt", type=float, default=1e-2, help="makroqadam, s")
    p.add_argument("-d", type=float, default=None, help="zarracha diametri, m")
    p.add_argument("--rho-p", type=float, default=2600.0, help="zarracha zichligi, kg/m3")
    p.add_argument("--harorat", type=float, default=293.15, help="havo harorati, K")
    p.add_argument("--qarshilik", choices=["re", "stoks"], default="re",
                   help="qarshilik qonuni: re — f(Re)=1+0.15Re^0.687 (standart), stoks — f=1")
    p.add_argument("--eps", type=float, default=1e-10, help="iteratsiya aniqligi (3.27)")
    p.add_argument("--boshlangich", choices=["harakatsiz", "impuls"], default="harakatsiz",
                   help="boshlang'ich shart: v(0)=0 yoki (3.5) v_x(0)=u, v_z(0)=0")
    p.add_argument("--fraksiya-fayl", default=None, help="fraksiyalar JSON fayli")
    p.add_argument("-o", "--chiqish", default=None, help="natijalarni CSV ga saqlash")
    p.add_argument("--taqqoslash", action="store_true",
                   help="soddalashtirilgan (v_p = u) model bilan taqqoslash")
    p.add_argument("-v", "--version", action="version", version=f"DustVel {__version__}")
    a = p.parse_args(argv)

    muhit = Muhit(ux=a.ux, uy=a.uy, T_muhit=a.harorat, profil=a.profil, z0=a.z0, h_s=a.hs)
    h = Hisoblagich(muhit, qarshilik=a.qarshilik, eps=a.eps)
    frs = _fraksiyalar(a)

    _sarlavha()
    print(f"  Shamol tezligi      : ({a.ux:.2f}, {a.uy:.2f}) m/s,  profil: {a.profil}")
    print(f"  Chiqarish balandligi: {a.h0:.1f} m")
    print(f"  Havo harorati       : {a.harorat:.2f} K")
    print(f"  Qovushqoqlik        : {muhit.qovushqoqlik():.4e} Pa*s")
    print(f"  Qarshilik           : "
          f"{'f(Re) = 1 + 0.15 Re^0.687' if a.qarshilik == 're' else 'Stoks (f = 1)'}")
    print(f"  Makroqadam          : {a.dt:.1e} s   (kichik qadam {a.dt/3:.3e} s)")
    print(BAR)
    print(f"{'Fraksiya':<14}{'tau, ms':>10}{'w_g, m/s':>11}{'w_s, m/s':>11}{'farq,%':>8}"
          f"{'Re_s':>7}{'iter':>6}{'t_yer, s':>10}{'x_yer, m':>10}")
    print("-" * 86)

    natijalar = []
    for fr in frs:
        n = h.hisobla(fr, T=a.vaqt, dt=a.dt, H0=a.h0, boshlangich=a.boshlangich)
        natijalar.append((fr, n))
        wg = h.stoks_tezligi(fr); ws = n.vt
        ty = f"{n.yer_t:.1f}" if n.yer_t else "> T"
        xy = f"{n.yer_x:.1f}" if n.yer_x else "-"
        print(f"{fr.nom:<14}{h.tau(fr)*1e3:>10.4f}{wg:>11.3e}{ws:>11.3e}"
              f"{(wg-ws)/wg*100:>8.2f}{h.reynolds(fr):>7.2f}{n.iter_ortacha:>6.2f}"
              f"{ty:>10}{xy:>10}")
    print(BAR)
    print("  w_g — Stoks qonuni bo'yicha (2.13), w_s — muvozanat cho'kish tezligi (3.6);")
    print("  iter — kichik qadamdagi o'rtacha iteratsiyalar soni (3.26)")

    if a.taqqoslash:
        print("\n  Soddalashtirilgan model (v_p = u) bilan taqqoslash")
        print("-" * 86)
        print(f"{'Fraksiya':<14}{'u*tau, m':>12}{'max |dx|, m':>14}{'nisbiy, %':>12}")
        print("-" * 86)
        for fr in frs:
            r = tahlil_qil(h, fr, T=min(a.vaqt, 10.0), dt=a.dt, H0=a.h0)
            print(f"{fr.nom:<14}{r['dx_asimptotik']:>12.4f}"
                  f"{r['dx_max']:>14.4f}{r['nisbiy_oxirgi']:>12.4f}")
        print(BAR)

    for fr in frs:
        Re = h.reynolds(fr)
        if a.qarshilik == "stoks" and Re > 0.3:
            print(f"  [ogohlantirish] {fr.nom}: Re = {Re:.2f}; Stoks qonuni cho'kish tezligini "
                  f"{(h.stoks_tezligi(fr)/Hisoblagich(muhit).chokish_tezligi(fr)-1)*100:.1f} % "
                  f"oshirib baholaydi — --qarshilik re tavsiya etiladi")
        if Re > 800.0:
            print(f"  [ogohlantirish] {fr.nom}: Re = {Re:.0f} > 800, f(Re) formulasi "
                  f"qo'llanish sohasidan tashqarida")

    if a.chiqish:
        os.makedirs(os.path.dirname(a.chiqish) or ".", exist_ok=True)
        with open(a.chiqish, "w", encoding="utf-8") as f:
            f.write("fraksiya,t,vx,vy,vz,x,y,z,Re\n")
            for fr, n in natijalar:
                step = max(1, len(n.t) // 2000)
                for i in range(0, len(n.t), step):
                    f.write(f"{fr.nom},{n.t[i]:.6f},{n.vx[i]:.6f},{n.vy[i]:.6f},"
                            f"{n.vz[i]:.6f},{n.x[i]:.4f},{n.y[i]:.4f},{n.z[i]:.4f},"
                            f"{n.Re[i]:.4f}\n")
        print(f"  Natijalar saqlandi: {a.chiqish}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
