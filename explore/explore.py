"""Prvotni explorace: rocni prehled, teplotni zavislost, denni profily."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import etl

OUT = Path(__file__).resolve().parent

# paleta (dataviz reference, light mode)
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "sans-serif",
})

df = etl.load()
daily = df.groupby("date_local").agg(
    baseload=("baseload", "mean"), temp=("temp", "mean"),
    sun=("sun", "mean"), wind=("wind", "mean"), daytype=("daytype", "max"),
)
daily.index = np.array(daily.index, dtype="datetime64[D]")

# --- 1. rocni prehled: spotreba a teplota (dva panely, spolecna casova osa) ---
fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 5.5), sharex=True,
                             gridspec_kw={"hspace": 0.12})
a1.plot(daily.index, daily.baseload, color=BLUE, lw=1.2)
a1.set_ylabel("spotřeba (denní průměr)")
a2.plot(daily.index, daily.temp, color=ORANGE, lw=1.2)
a2.set_ylabel("teplota [°C] (denní průměr)")
a2.axhline(0, color=BASELINE, lw=0.8)
fig.suptitle("Roční přehled: denní průměry spotřeby a teploty", x=0.12, ha="left")
fig.savefig(OUT / "01_prehled.png", dpi=130, bbox_inches="tight")

# --- 2. energeticka signatura: denni spotreba vs. denni teplota, po typech dne ---
fig, ax = plt.subplots(figsize=(8.5, 6))
labels = ["Po–Čt", "pátek/most", "sobota", "neděle/svátek"]
# linie snesou 4 sloty (adjacent pairlist), scatter jen 3 (all-pairs) -> slouceni prac. dnu
DTSET = [(0, BLUE, labels[0]), (1, ORANGE, labels[1]), (2, AQUA, labels[2]), (3, YELLOW, labels[3])]
for sel, color, lab in [(daily.daytype <= 1, BLUE, "pracovní (vč. Pá)"),
                        (daily.daytype == 2, ORANGE, labels[2]),
                        (daily.daytype == 3, AQUA, labels[3])]:
    d = daily[sel]
    ax.scatter(d.temp, d.baseload, s=14, color=color, alpha=0.65, lw=0, label=lab)
ax.set_xlabel("denní průměrná teplota [°C]")
ax.set_ylabel("denní průměrná spotřeba")
ax.set_title("Energetická signatura portfolia", loc="left")
ax.legend(frameon=False)
fig.savefig(OUT / "02_signatura.png", dpi=130, bbox_inches="tight")

# --- 3. teplotni zavislost po hodinach dne (pracovni dny) ---
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
for ax, hh in zip(axes, [3, 12, 19]):
    d = df[(df.daytype == 0) & (df.tod == hh)]
    ax.scatter(d.temp, d.baseload, s=8, color=BLUE, alpha=0.5, lw=0)
    ax.set_title(f"{hh}:00", loc="left", color=INK2)
    ax.set_xlabel("teplota [°C]")
axes[0].set_ylabel("spotřeba (15min)")
fig.suptitle("Teplotní závislost po hodinách dne — pracovní dny", x=0.09, ha="left")
fig.savefig(OUT / "03_po_hodinach.png", dpi=130, bbox_inches="tight")

# --- 4. denni profily: zima vs. leto, po typech dne ---
month = df.index.tz_convert(etl.TZ).month
seasons = {"zima (12–2)": month.isin([12, 1, 2]), "léto (6–8)": month.isin([6, 7, 8])}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
for ax, (season, mask) in zip(axes, seasons.items()):
    for dt, color, lab in DTSET:
        prof = df[mask & (df.daytype == dt)].groupby("tod")["baseload"].mean()
        ax.plot(prof.index, prof.values, color=color, lw=2, label=lab)
    ax.set_title(season, loc="left", color=INK2)
    ax.set_xlabel("čas dne [h]")
    ax.set_xticks([0, 6, 12, 18, 24])
axes[0].set_ylabel("průměrná spotřeba")
axes[0].legend(frameon=False)
fig.suptitle("Denní profily podle typu dne", x=0.09, ha="left")
fig.savefig(OUT / "04_profily.png", dpi=130, bbox_inches="tight")

# --- 5. vliv osvitu: signatura obarvena osvitem (topna sezona) ---
fig, ax = plt.subplots(figsize=(8.5, 6))
d = daily[daily.temp < 15]
sc = ax.scatter(d.temp, d.baseload, s=16, c=d.sun, cmap=plt.matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ),
                alpha=0.85, lw=0)
cb = fig.colorbar(sc, ax=ax, label="osvit (denní průměr) [W/m²]")
cb.outline.set_visible(False)
ax.set_xlabel("denní průměrná teplota [°C]")
ax.set_ylabel("denní průměrná spotřeba")
ax.set_title("Topná sezóna: role osvitu (dny < 15 °C)", loc="left")
fig.savefig(OUT / "05_osvit.png", dpi=130, bbox_inches="tight")

# tistene souhrny pro report
dead = daily[(daily.temp > 15) & (daily.temp < 20) & (daily.daytype == 0)]
cold = daily[(daily.temp < 0) & (daily.daytype == 0)]
hot = daily[(daily.temp > 25) & (daily.daytype == 0)]
print("prac. dny v mrtvem pasmu 15-20C:", len(dead), "| prumer spotreby:", round(dead.baseload.mean(), 1))
print("prac. dny pod 0C:", len(cold), "| prumer:", round(cold.baseload.mean(), 1))
print("prac. dny nad 25C:", len(hot), "| prumer:", round(hot.baseload.mean(), 1))
print("pomer vikend/prac (mrtve pasmo):",
      round(daily[(daily.temp.between(15, 20)) & (daily.daytype == 3)].baseload.mean()
            / max(dead.baseload.mean(), 1e-9), 3))
