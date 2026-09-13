# -*- coding: utf-8 -*-
"""
DustVel — chang zarrachalarining atmosferadagi harakat tezliklarini
hisoblovchi dasturiy vosita.

Hisoblash yadrosi uch bosqichli sxemaga asoslangan:
    n -> n+1/3       : Krank-Nikolson
    n+1/3 -> n+2/3   : BDF2
    n+2/3 -> n+1     : BDF2
Vaqt bo'yicha aniqlik tartibi O(dt^2), sxema shartsiz turg'un.
"""
__version__ = "1.0.0"

from .core import Fraksiya, Muhit, Hisoblagich, Natija, tahlil_qil

STANDART_FRAKSIYALAR = [
    Fraksiya("PM2.5", 2.5e-6, 1500.0),
    Fraksiya("PM10", 10.0e-6, 2000.0),
    Fraksiya("Yirik (50 um)", 50.0e-6, 2500.0),
]

__all__ = ["Fraksiya", "Muhit", "Hisoblagich", "Natija",
           "tahlil_qil", "STANDART_FRAKSIYALAR", "__version__"]
