"""Je cast spotreby rizena svetlem (slunce), ne hodinami? Analyza bokem modelu.

Hypoteza: spotrebu ridi dva casove procesy — **rozvrh** (prace, skoly, obchody;
drzi se hodin) a **denni svetlo** (verejne i vnitrni osvetleni; drzi se
slunce). Testuje se na spotrebe ocistene o pocasi: y - topeni - chlazeni + FVE
(joint parametry modelu), takze zbyva baze + sum.

1. **Zmena casu** (prirozeny experiment): pres noc se hodiny posunou o hodinu
   vuci slunci, pocasi, sezona i rozvrh zustanou. Slozka rizena svetlem se na
   hodinach posune, slozka rizena rozvrhem ne. Rozdil profilu 14 dni po a pred
   zmenou minus stejny rozdil o 14 dni driv (odecte sezonni drift).
2. **Mesicni profily** pracovnich dnu na ose hodin a na ose "cas od zapadu
   slunce" — pokud vecerni narust sleduje slunce, na slunecni ose se srovnaji.
3. **Velikost po letech**: soumrakova slozka z kazde zmeny casu zvlast
   (usporne zdroje by ji v case zmensovaly). Regrese "tma x aktivita(cas dne)"
   pres cely den se neosvedcila — v poledne tma nenastava, krivka tam neni
   identifikovana; pro clen modelu omezit bazi na okna soumraku a usvitu.

Vychod/zapad a vyska slunce: src/sun.py (NOAA aproximace pro stred CR).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import base
import cooling
import etl
import fit
import heating
import pv
from sun import darkness, elevation, sun_times

OUT = Path(__file__).resolve().parent

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "sans-serif",
})


# --- data ----------------------------------------------------------------------

df = etl.load()
bp = base.load(base.MODEL_DIR / "base_joint.npz")
hp = heating.load(heating.MODEL_DIR / "heating_joint.npz")
cp = fit.load("cooling_joint.npz", cooling.PARAM_NAMES + ["rmse", "r2"])
pp = fit.load("pv_joint.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"])
yc = (df["baseload"].to_numpy() - heating.predict(df, hp) - cooling.predict(df, cp)
      + pv.predict(df, pp))
mid = df.index + pd.Timedelta(minutes=7.5)          # stred intervalu
loc = df.index.tz_convert(etl.TZ)
d = pd.DataFrame({
    "yc": yc, "den": pd.to_datetime(df["date_local"]).to_numpy(), "tod": df["tod"].to_numpy(),
    "slot": (df["tod"].to_numpy() * 4).round().astype(int), "dt": df["daytype"].to_numpy(),
    "po_volnu": df["po_volnu"].to_numpy(), "rok": loc.year, "mes": loc.month,
    "elev": elevation(mid), "pz": df["prazdniny"].to_numpy(),
}, index=df.index)
d["tma"] = darkness(d["elev"].to_numpy())
# odchylka od poledni urovne dne (10-15 h je svetlo cely rok) — odstrani uroven
noon = d[(d.tod >= 10) & (d.tod < 15)].groupby("den")["yc"].mean()
d["dev"] = d["yc"] - d["den"].map(noon)
regular = (d.dt == 0) & ~d.po_volnu & ~d.pz      # Ut-Ct mimo prazdniny

# --- 1. zmena casu -------------------------------------------------------------

def dst_dates() -> list[tuple[pd.Timestamp, str]]:
    out = []
    for y in range(2022, 2027):
        for m, kind in ((3, "jaro"), (10, "podzim")):
            last = pd.Timestamp(y, m, 31)
            sun = last - pd.Timedelta(days=(last.dayofweek + 1) % 7)
            out.append((sun, kind))
    return out


def window_profile(a: pd.Timestamp, b: pd.Timestamp) -> pd.DataFrame:
    sel = regular & (d.den >= a) & (d.den <= b)
    return d[sel].groupby(["den", "slot"])["dev"].mean().unstack()


W = 14
events = []
for sw, kind in dst_dates():
    if sw - pd.Timedelta(days=2 * W) < d.den.min() or sw + pd.Timedelta(days=W) > d.den.max():
        continue
    pre = window_profile(sw - pd.Timedelta(days=2 * W), sw - pd.Timedelta(days=W + 1))
    bef = window_profile(sw - pd.Timedelta(days=W), sw - pd.Timedelta(days=1))
    aft = window_profile(sw + pd.Timedelta(days=1), sw + pd.Timedelta(days=W))
    did = (aft.mean() - bef.mean()) - (bef.mean() - pre.mean())
    plac = bef.mean() - pre.mean()
    events.append({"datum": sw, "kind": kind, "did": did, "placebo": plac,
                   "se": np.sqrt(aft.var() / len(aft) + 4 * bef.var() / len(bef) + pre.var() / len(pre))})

print("zmeny casu v datech:", [(e["datum"].date(), e["kind"]) for e in events])
rise, sset = sun_times(pd.date_range(d.den.min(), d.den.max(), freq="D"))

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
for ax, kind, col in ((axes[0], "jaro", AQUA), (axes[1], "podzim", ORANGE)):
    ev = [e for e in events if e["kind"] == kind]
    did = pd.concat([e["did"] for e in ev], axis=1).mean(axis=1)
    se = np.sqrt(pd.concat([e["se"] ** 2 for e in ev], axis=1).sum(axis=1)) / len(ev)
    plac = pd.concat([e["placebo"] for e in ev], axis=1).mean(axis=1)
    x = did.index / 4
    ax.fill_between(x, did - 2 * se, did + 2 * se, color=col, alpha=0.18, lw=0)
    ax.plot(x, did, color=col, lw=2, label="rozdíl rozdílů kolem změny času")
    ax.plot(x, plac - plac.mean(), color=MUTED, lw=1, label="placebo: stejné okno o 14 dní dřív")
    ax.axhline(0, color=BASELINE, lw=1)
    # vychod/zapad v mistnim case tesne pred a po zmene
    for e in ev:
        for off, ls in ((-3, ":"), (3, "--")):
            day = e["datum"] + pd.Timedelta(days=off)
            for tt in (rise[day], sset[day]):
                tl = tt.tz_convert(etl.TZ)
                ax.axvline(tl.hour + tl.minute / 60, color=INK2, lw=0.6, ls=ls, alpha=0.35)
    ax.set_title(f"{kind}: {len(ev)} změn ({', '.join(str(e['datum'].year) for e in ev)})", loc="left")
    ax.set_xlabel("místní čas [h]  (svislé čáry: východ/západ slunce před ⋯ a po - - změně)")
    ax.set_xticks(range(0, 25, 3))
axes[0].set_ylabel("změna profilu (po − před) − (před − dřív)")
axes[0].legend(frameon=False, fontsize=9, loc="upper left")
fig.suptitle("1  Změna času: posune se část spotřeby se sluncem?  (Út–Čt, očištěno o počasí, ±2 SE)", x=0.01, ha="left")
fig.tight_layout()
fig.savefig(OUT / "26_svetlo_zmena_casu.png", dpi=130, bbox_inches="tight")

# --- 2. mesicni profily: hodiny vs. cas od zapadu --------------------------------

# vcetne prazdnin (jinak chybi cervenec a srpen)
mp = d[(d.dt == 0) & ~d.po_volnu].groupby(["mes", "slot"])["dev"].mean().unstack()
sun_loc = sset.dt.tz_convert(etl.TZ)
sset_h = (sun_loc.dt.hour + sun_loc.dt.minute / 60).groupby(sun_loc.index.month).mean()
cmap = plt.get_cmap("twilight")
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for m in range(1, 13):
    c = cmap((m - 0.5) / 12)
    x = mp.columns / 4
    axes[0].plot(x, mp.loc[m], color=c, lw=1.4, label=f"{m}")
    axes[1].plot(x - sset_h[m], mp.loc[m], color=c, lw=1.4)
axes[0].set_xlim(12, 24); axes[0].set_xlabel("místní čas [h]")
axes[1].set_xlim(-6, 6); axes[1].set_xlabel("hodiny od západu slunce")
axes[1].axvline(0, color=INK2, lw=0.8)
axes[0].set_ylabel("odchylka od polední úrovně dne")
axes[0].legend(title="měsíc", ncol=2, frameon=False, fontsize=8)
axes[0].set_title("na ose hodin", loc="left"); axes[1].set_title("na ose slunce (zarovnáno na západ)", loc="left")
fig.suptitle("2  Odpolední a večerní profil Út–Čt po měsících (očištěno o počasí)", x=0.01, ha="left")
fig.tight_layout()
fig.savefig(OUT / "27_svetlo_mesice.png", dpi=130, bbox_inches="tight")

# rozptyl mesicnich profilu ve vecernim okne: hodiny vs slunce
def spread(shift: bool) -> float:
    grid = np.arange(-3, 3.01, 0.25) if shift else np.arange(15, 21.01, 0.25)
    rows = []
    for m in range(1, 13):
        x = mp.columns / 4 - (sset_h[m] if shift else 0)
        rows.append(np.interp(grid, x, mp.loc[m].to_numpy()))
    R = np.array(rows)
    return float(np.mean(R.std(axis=0)))


print(f"rozptyl mesicnich vecernich profilu: osa hodin {spread(False):.1f}, osa slunce {spread(True):.1f}")

print("\n1) rozdil rozdilu (prumer pres zmeny), vybrane hodiny:")
for kind in ("jaro", "podzim"):
    ev = [e for e in events if e["kind"] == kind]
    did = pd.concat([e["did"] for e in ev], axis=1).mean(axis=1)
    pl = pd.concat([e["placebo"] for e in ev], axis=1).mean(axis=1)
    hrs = [6, 7, 8, 16, 17, 18, 19, 20, 21]
    print(f"   {kind:6s} DiD: " + " ".join(f"{h}h {did[did.index // 4 == h].mean():+5.1f}" for h in hrs))
    print(f"   {kind:6s} plac:" + " ".join(f"{h}h {(pl - pl.mean())[pl.index // 4 == h].mean():+5.1f}" for h in hrs))

# velikost "soumrakove" slozky po jednotlivych zmenach casu: prumer DiD v okne
# mezi zapadem pred a po zmene (tam se tma/svetlo prohodi), znamenko tak, aby
# kladna hodnota = spotreba navic za tmy
print("\n4) soumrakova slozka po zmenach casu (prumer v hodine mezi zapady, kladne = navic za tmy):")
rows = []
for e in events:
    a = sset[e["datum"] - pd.Timedelta(days=3)].tz_convert(etl.TZ)
    b = sset[e["datum"] + pd.Timedelta(days=3)].tz_convert(etl.TZ)
    lo, hi = sorted([a.hour + a.minute / 60, b.hour + b.minute / 60])
    sl = (e["did"].index / 4 >= lo) & (e["did"].index / 4 < hi)
    sign = -1 if e["kind"] == "jaro" else 1     # na jare se tma posune pozdeji
    v = sign * e["did"][sl].mean()
    se = float(np.sqrt((e["se"][sl] ** 2).mean() / sl.sum()))
    rows.append((e["datum"].date(), e["kind"], v, se))
    print(f"   {e['datum'].date()} {e['kind']:6s} {lo:5.2f}-{hi:5.2f} h: {v:+6.1f} (SE {se:4.1f})")
