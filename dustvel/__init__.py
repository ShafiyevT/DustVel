# -*- coding: utf-8 -*-
"""
DustVel — chang zarrachalarining atmosferadagi harakat tezliklarini
hisoblovchi dasturiy vosita.

Hisoblash yadrosi (dissertatsiya, III bob, § 3.2):
    n -> n+1/3       : simmetrik oshkormas sxema (Krank-Nikolson)
    n+1/3 -> n+2/3   : ikkinchi tartibli uch qatlamli oshkormas sxema (BDF2)
    n+2/3 -> n+1     : ikkinchi tartibli uch qatlamli oshkormas sxema (BDF2)
Qarshilik Reynolds soniga bog'liq: f(Re) = 1 + 0.15 Re^0.687; chiziqsiz
tenglamalar komponentalar bo'yicha iteratsiya bilan yechiladi.
Vaqt bo'yicha aniqlik tartibi O(dt^2), sxema A- va L-turg'un.
"""
__version__ = "2.0.0"

from .core import (Fraksiya, Muhit, Hisoblagich, Natija, tahlil_qil,
                   f_Re, rk4_etalon, aniq_yechim_stoks)

# Dissertatsiyada qabul qilingan besh fraksiya (II bob, 2.2-jadval):
# o'rta geometrik diametr, zichlik va emissiyadagi massa ulushi p_i (2.25)
STANDART_FRAKSIYALAR = [
    Fraksiya("0.2-2.5 um",  0.707e-6, 1500.0, 0.049),
    Fraksiya("2.5-10 um",   5.000e-6, 2000.0, 0.327),
    Fraksiya("10-20 um",   14.142e-6, 2200.0, 0.088),
    Fraksiya("20-50 um",   31.623e-6, 2400.0, 0.327),
    Fraksiya("50-100 um",  70.711e-6, 2600.0, 0.208),
]

__all__ = ["Fraksiya", "Muhit", "Hisoblagich", "Natija", "tahlil_qil",
           "f_Re", "rk4_etalon", "aniq_yechim_stoks",
           "STANDART_FRAKSIYALAR", "__version__"]
