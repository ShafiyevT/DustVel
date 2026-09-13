# -*- coding: utf-8 -*-
"""python -m dustvel [--gui]"""
import sys

if "--gui" in sys.argv:
    from .gui import ishga_tushir
    ishga_tushir()
else:
    from .cli import main
    sys.exit(main())
