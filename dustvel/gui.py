# -*- coding: utf-8 -*-
"""DustVel — grafik foydalanuvchi interfeysi (matplotlib widget asosida)."""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons

from .core import Fraksiya, Muhit, Hisoblagich
from . import STANDART_FRAKSIYALAR, __version__

RANG = {"PM2.5": "#1f6fb2", "PM10": "#2c7a4b", "Yirik (50 um)": "#c0392b"}


class Oyna:
    """Asosiy ish oynasi."""

    def __init__(self):
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
        self.fig = plt.figure(figsize=(13.6, 7.6))
        self.fig.canvas.manager.set_window_title(
            f"DustVel {__version__} — chang zarrachalari harakati")
        self.fig.patch.set_facecolor("#f4f6f8")

        # --- sarlavha paneli ---
        self.fig.text(0.012, 0.962, "DustVel", fontsize=17, fontweight="bold",
                      color="#14507d")
        self.fig.text(0.012, 0.935,
                      "Chang zarrachalarining atmosferadagi harakat tezliklarini hisoblash",
                      fontsize=9.5, color="#333")
        self.fig.text(0.012, 0.915, "Buxoro davlat universiteti", fontsize=8.5,
                      color="#666", style="italic")
        self.fig.patches.append(plt.Rectangle((0, 0.905), 1, 0.003,
                                              transform=self.fig.transFigure,
                                              color="#14507d", zorder=0.5))

        # --- boshqaruv paneli ---
        self.fig.patches.append(plt.Rectangle((0.008, 0.022), 0.245, 0.868,
                                              transform=self.fig.transFigure,
                                              facecolor="white", edgecolor="#c8d2da",
                                              zorder=0))
        self.fig.text(0.020, 0.862, "KIRISH PARAMETRLARI", fontsize=9.5,
                      fontweight="bold", color="#14507d")

        self.s_u = Slider(self._ax([0.055, 0.812, 0.175, 0.020]),
                          "shamol\n$u_x$, m/s", 0.5, 20.0, valinit=5.0, valfmt="%.1f",
                          color="#1f6fb2")
        self.s_h = Slider(self._ax([0.055, 0.762, 0.175, 0.020]),
                          "balandlik\n$H_0$, m", 1.0, 60.0, valinit=20.0, valfmt="%.0f",
                          color="#2c7a4b")
        self.s_T = Slider(self._ax([0.055, 0.712, 0.175, 0.020]),
                          "harorat\n$T$, °C", -10.0, 45.0, valinit=20.0, valfmt="%.0f",
                          color="#d77a3a")
        self.s_dt = Slider(self._ax([0.055, 0.662, 0.175, 0.020]),
                           "qadam\n$\\log_{10}\\Delta t$", -4.0, -1.0, valinit=-3.0,
                           valfmt="%.1f", color="#8e44ad")

        self.fig.text(0.020, 0.624, "FRAKSIYALAR", fontsize=9.5,
                      fontweight="bold", color="#14507d")
        ax_ch = self._ax([0.030, 0.512, 0.200, 0.100])
        ax_ch.set_facecolor("#fbfcfd")
        self.chk = CheckButtons(ax_ch, [f.nom for f in STANDART_FRAKSIYALAR],
                                [True, True, True])

        self.fig.text(0.020, 0.478, "BOSHLANGʻICH SHART", fontsize=9.5,
                      fontweight="bold", color="#14507d")
        ax_bs = self._ax([0.030, 0.406, 0.200, 0.060])
        ax_bs.set_facecolor("#fbfcfd")
        self.rb_bs = RadioButtons(ax_bs, ["harakatsiz  v(0)=0", "impuls  v(0)=u"], active=0)

        self.fig.text(0.020, 0.372, "KO'RINISH", fontsize=9.5,
                      fontweight="bold", color="#14507d")
        ax_rb = self._ax([0.030, 0.262, 0.200, 0.098])
        ax_rb.set_facecolor("#fbfcfd")
        self.rb = RadioButtons(ax_rb, ["Tezliklar", "Traektoriya", "Taqqoslash"], active=0)

        self.b_hisob = Button(self._ax([0.030, 0.202, 0.095, 0.038]),
                              "Hisoblash", color="#1f6fb2", hovercolor="#155a92")
        self.b_hisob.label.set_color("white"); self.b_hisob.label.set_fontweight("bold")
        self.b_saqla = Button(self._ax([0.135, 0.202, 0.095, 0.038]),
                              "CSV saqlash", color="#e8eef3", hovercolor="#d3dde5")

        self.fig.text(0.020, 0.166, "HISOBLANGAN KATTALIKLAR", fontsize=9.5,
                      fontweight="bold", color="#14507d")
        self.txt = self.fig.text(0.022, 0.038, "", fontsize=6.5, family="DejaVu Sans Mono",
                                 va="bottom", color="#222")

        # --- grafik maydonlari ---
        self.ax1 = self.fig.add_axes([0.315, 0.565, 0.312, 0.310])
        self.ax2 = self.fig.add_axes([0.672, 0.565, 0.312, 0.310])
        self.ax3 = self.fig.add_axes([0.315, 0.120, 0.669, 0.335])

        self.status = self.fig.text(0.315, 0.020, "", fontsize=8.5, color="#2c7a4b")

        for w in (self.s_u, self.s_h, self.s_T, self.s_dt):
            w.on_changed(lambda v: None)
        self.b_hisob.on_clicked(self.hisobla)
        self.rb.on_clicked(lambda l: self.hisobla(None))
        self.rb_bs.on_clicked(lambda l: self.hisobla(None))
        self.chk.on_clicked(lambda l: self.hisobla(None))

        self.natijalar = []
        self.hisobla(None)

    def _ax(self, rect, z=3):
        a = self.fig.add_axes(rect)
        a.set_zorder(z)
        return a

    # ------------------------------------------------------------------
    def _tanlangan(self):
        st = self.chk.get_status()
        return [f for f, s in zip(STANDART_FRAKSIYALAR, st) if s]

    def hisobla(self, event):
        ux = self.s_u.val
        H0 = self.s_h.val
        Tc = self.s_T.val + 273.15
        dt = 10 ** self.s_dt.val

        muhit = Muhit(ux=ux, T_muhit=Tc)
        h = Hisoblagich(muhit)
        frs = self._tanlangan()

        self.natijalar = []
        for fr in frs:
            bs = "impuls" if self.rb_bs.value_selected.startswith("impuls") else "harakatsiz"
            n = h.hisobla(fr, T=300.0, dt=dt, H0=H0, boshlangich=bs)
            self.natijalar.append((fr, n, h))

        self._chiz(h, frs)
        self._matn(h, frs)
        self.status.set_text(
            f"Hisoblash yakunlandi  |  qovushqoqlik {muhit.qovushqoqlik():.3e} Pa·s  |  "
            f"makroqadam {dt:.1e} s, kichik qadam {dt/3:.2e} s")
        self.fig.canvas.draw_idle()

    def _chiz(self, h, frs):
        rejim = self.rb.value_selected
        for ax in (self.ax1, self.ax2, self.ax3):
            ax.clear()

        # 1-panel: gorizontal tezlik relaksatsiyasi
        for fr, n, _ in self.natijalar:
            m = n.t <= 0.15
            self.ax1.plot(n.t[m] * 1e3, n.vx[m], lw=1.9,
                          color=RANG.get(fr.nom, "#555"), label=fr.nom)
        self.ax1.axhline(h.muhit.ux, ls=":", lw=1.1, color="#888")
        self.ax1.set_xlabel("$t$, ms"); self.ax1.set_ylabel("$v_x$, m/s")
        self.ax1.set_title("Gorizontal tezlik", fontsize=10)
        self.ax1.grid(alpha=0.3); self.ax1.legend(fontsize=7.5)

        # 2-panel: vertikal cho'kish tezligi
        for fr, n, _ in self.natijalar:
            m = n.t > 0
            self.ax2.loglog(n.t[m], np.abs(n.vz[m]), lw=1.9,
                            color=RANG.get(fr.nom, "#555"), label=fr.nom)
        self.ax2.set_xlabel("$t$, s"); self.ax2.set_ylabel("$|v_z|$, m/s")
        self.ax2.set_title("Vertikal cho'kish tezligi", fontsize=10)
        self.ax2.grid(alpha=0.3, which="both"); self.ax2.legend(fontsize=7.5)

        # 3-panel: rejimga qarab
        if rejim == "Tezliklar":
            for fr, n, _ in self.natijalar:
                m = n.t <= 1.0
                self.ax3.plot(n.t[m], n.vz[m], lw=2.0,
                              color=RANG.get(fr.nom, "#555"), label=fr.nom)
            self.ax3.set_xlabel("$t$, s"); self.ax3.set_ylabel("$v_z$, m/s")
            self.ax3.set_title("Vertikal tezlikning barqarorlashuvi", fontsize=10)
        elif rejim == "Traektoriya":
            for fr, n, _ in self.natijalar:
                self.ax3.plot(n.x, n.z, lw=2.0,
                              color=RANG.get(fr.nom, "#555"), label=fr.nom)
            self.ax3.axhline(0, color="#8a6a3c", lw=2.0)
            self.ax3.set_xlabel("$x$, m"); self.ax3.set_ylabel("$z$, m")
            self.ax3.set_title("Zarrachalar traektoriyasi", fontsize=10)
            self.ax3.set_ylim(bottom=-1)
        else:
            for fr, n, hh in self.natijalar:
                t2, X2, _ = hh.soddalashtirilgan(fr, T=n.t[-1], dt=n.t[1] - n.t[0],
                                                 H0=n.z[0])
                m = min(len(n.t), len(t2))
                dx = np.abs(X2[:m] - n.x[:m])
                rel = dx / np.maximum(n.x[:m], 1e-9) * 100
                sel = n.t[:m] > 1e-3
                self.ax3.loglog(n.t[:m][sel], rel[sel], lw=2.0,
                                color=RANG.get(fr.nom, "#555"), label=fr.nom)
            self.ax3.set_xlabel("$t$, s")
            self.ax3.set_ylabel("nisbiy farq, %")
            self.ax3.set_title("Soddalashtirilgan model ($v_p=u$) bilan farq", fontsize=10)
        self.ax3.grid(alpha=0.3, which="both"); self.ax3.legend(fontsize=8)

    def _matn(self, h, frs):
        s = []
        s.append(f"{'fraksiya':<13}{'tau,ms':>8}{'v_t,m/s':>10}{'t_yer,s':>9}{'x_yer,m':>9}")
        s.append("-" * 49)
        for fr, n, _ in self.natijalar:
            ty = f"{n.yer_t:.1f}" if n.yer_t else ">300"
            xy = f"{n.yer_x:.0f}" if n.yer_x else "-"
            s.append(f"{fr.nom:<13}{h.tau(fr)*1e3:>8.3f}"
                     f"{h.chokish_tezligi(fr):>10.2e}{ty:>9}{xy:>9}")
        re_list = ", ".join(f"{fr.nom.split()[0]}: {h.reynolds(fr):.3f}"
                            for fr, _, _ in self.natijalar)
        s.append("")
        s.append(f"Re_p -> {re_list}")
        self.txt.set_text("\n".join(s))


def ishga_tushir():
    o = Oyna()
    plt.show()
    return o


if __name__ == "__main__":
    ishga_tushir()
