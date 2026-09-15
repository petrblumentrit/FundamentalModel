"""Validace mimo vzorek (krok 7): vynechany rok (interpolace) + dopredna (extrapolace).

Vysledky foldu se ukladaji do models/validace.pkl, aby se grafy daly
prekreslit bez opakovani fitu (`--plot-only`).
"""
import pickle
import sys
import time
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
import validate

OUT = Path(__file__).resolve().parent
CACHE = base.MODEL_DIR / "validace.pkl"

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

FOLDS = [("vynechany", 2023), ("vynechany", 2024), ("vynechany", 2025),
         ("dopredna", 2024), ("dopredna", 2025), ("dopredna", 2026)]

df = etl.load()
yr = validate.year_local(df)
mon = df.index.tz_convert(etl.TZ).month.to_numpy()
y = df["baseload"].to_numpy()

# reference: joint fit ze vsech dat (in-sample)
full = (base.load(base.MODEL_DIR / "base_joint.npz"),
        heating.load(heating.MODEL_DIR / "heating_joint.npz"),
        fit.load("cooling_joint.npz", cooling.PARAM_NAMES + ["rmse", "r2"]),
        fit.load("pv_joint.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"]))
pred_full = validate.predict(df, full)


def monthly_rmse(e: np.ndarray, m: np.ndarray) -> pd.Series:
    return pd.Series(e).groupby(m).apply(lambda x: np.sqrt((x**2).mean()))


def fmt(d, k, w):
    return f"{d[k]:>{w}.1f}" if d else " " * w


res = {}
if "--plot-only" in sys.argv and CACHE.exists():
    res = pickle.load(open(CACHE, "rb"))
else:
    for kind, year in FOLDS:
        train = (yr != year) if kind == "vynechany" else (yr < year)
        test = yr == year
        t = time.time()
        params = validate.fit_masked(df, train)
        pred = validate.predict(df, params)
        sc = validate.score(y[test], pred[test])
        res[(kind, year)] = {
            "params": params, "pred_test": pred[test],
            "score": sc, "n_train_days": int(train.sum() / 96),
            "monthly": monthly_rmse(y[test] - pred[test], mon[test]),
            "summary": validate.summary_params(params, df),
        }
        print(f"{kind:>10} {year}: trenink {train.sum()/96:.0f} dnu, "
              f"RMSE {sc['rmse']:.1f}, bias {sc['bias']:+.1f}, "
              f"bez biasu {sc['rmse_bez_biasu']:.1f}, R2 {sc['r2']:.3f}  "
              f"({time.time() - t:.0f} s)", flush=True)
        with open(CACHE, "wb") as fh:
            pickle.dump(res, fh)

# --- tabulky ---
years = sorted({yv for _, yv in FOLDS})
print(f"\n{'rok':<6}{'in-sample':>11}{'vynechany':>11}{'bias':>7}"
      f"{'dopredna':>10}{'bias':>7}{'bez biasu':>11}")
for yv in years:
    ins = validate.score(y[yr == yv], pred_full[yr == yv])["rmse"]
    lo = res.get(("vynechany", yv), {}).get("score")
    fw = res.get(("dopredna", yv), {}).get("score")
    print(f"{yv:<6}{ins:>11.1f}{fmt(lo, 'rmse', 11)}{fmt(lo, 'bias', 7)}"
          f"{fmt(fw, 'rmse', 10)}{fmt(fw, 'bias', 7)}{fmt(fw, 'rmse_bez_biasu', 11)}")

print("\nRMSE po mesicich (in-sample / vynechany rok / dopredna):")
for yv in years:
    ins = monthly_rmse(y[yr == yv] - pred_full[yr == yv], mon[yr == yv])
    lo = res.get(("vynechany", yv), {}).get("monthly")
    fw = res.get(("dopredna", yv), {}).get("monthly")
    row = [f"{yv}:"]
    for mm in ins.index:
        parts = [f"{ins[mm]:.0f}",
                 f"{lo[mm]:.0f}" if lo is not None else "-",
                 f"{fw[mm]:.0f}" if fw is not None else "-"]
        row.append(f"{mm:>2}:{'/'.join(parts)}")
    print("  " + " ".join(row))

keys = ["top_t_b", "top_s", "top_w", "top_tau", "top_alpha", "top_a",
        "chl_t_bc", "chl_s_c", "chl_w_c", "chl_tau_c", "chl_alpha_c",
        "fve_gamma", "fve_spicka"]
sum_full = validate.summary_params(full, df)
cols = [("vse", sum_full)] + [(f"{k[:3]}{yv}", r["summary"]) for (k, yv), r in res.items()]
print(f"\n{'parametr':<12}" + "".join(f"{c:>10}" for c, _ in cols))
for k in keys:
    print(f"{k:<12}" + "".join(f"{s[k]:>10.3f}" for _, s in cols))

# --- 22. mesicni RMSE: in-sample vs. mimo vzorek ---
fig, axes = plt.subplots(1, len(years), figsize=(3.4 * len(years), 3.8), sharey=True)
for ax, yv in zip(axes, years):
    ins = monthly_rmse(y[yr == yv] - pred_full[yr == yv], mon[yr == yv])
    x = ins.index.to_numpy()
    ax.bar(x - 0.28, ins.values, 0.26, color=BASELINE, label="in-sample")
    lo = res.get(("vynechany", yv))
    if lo is not None:
        ax.bar(x, lo["monthly"].reindex(x).values, 0.26, color=BLUE, label="vynechaný rok")
    fw = res.get(("dopredna", yv))
    if fw is not None:
        ax.bar(x + 0.28, fw["monthly"].reindex(x).values, 0.26, color=ORANGE, label="dopředná")
    ax.set_title(str(yv), loc="left")
    ax.set_xticks(x)
    ax.set_xlabel("měsíc")
axes[0].set_ylabel("RMSE")
handles, lbls = axes[1].get_legend_handles_labels()
fig.legend(handles, lbls, frameon=False, ncol=3, loc="upper right", bbox_to_anchor=(0.99, 1.04))
fig.suptitle("Chyba mimo vzorek po měsících", x=0.01, y=1.03, ha="left")
fig.savefig(OUT / "22_validace_mesice.png", dpi=130, bbox_inches="tight")

# --- 23. stabilita parametru pres foldy ---
labels = {"top_t_b": "T_b [°C]", "top_s": "s [°C]", "top_w": "w", "top_tau": "τ [h]",
          "top_alpha": "α COP", "top_a": "a osvit", "chl_t_bc": "T_bc [°C]",
          "chl_s_c": "s_c [°C]", "chl_tau_c": "τ_c [h]", "chl_w_c": "w_c",
          "chl_alpha_c": "α EER", "fve_spicka": "FVE špička"}
fig, axes = plt.subplots(3, 4, figsize=(12.5, 8))
fold_names = [f"{'V' if k == 'vynechany' else 'D'}{yv}" for k, yv in res]
colors = [BLUE if k == "vynechany" else ORANGE for k, _ in res]
for ax, key in zip(axes.flat, labels):
    vals = [r["summary"][key] for r in res.values()]
    ax.axhline(sum_full[key], color=INK2, lw=1.2, ls="--")
    ax.scatter(range(len(vals)), vals, c=colors, s=40, zorder=3)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(fold_names, rotation=45, fontsize=8)
    ax.set_title(labels[key], loc="left", fontsize=10)
fig.suptitle("Parametry tvaru přes foldy (čárkovaně = fit ze všech dat; "
             "modře vynechaný rok, oranžově dopředná)", x=0.01, ha="left")
fig.tight_layout()
fig.savefig(OUT / "23_validace_parametry.png", dpi=130, bbox_inches="tight")

# --- 24. prubeh testovaneho roku: mereno vs. predpoved ---
show = [yv for yv in (2025, 2026) if ("dopredna", yv) in res]
fig, axes = plt.subplots(len(show), 1, figsize=(11.5, 3.6 * len(show)), squeeze=False)
for ax, yv in zip(axes[:, 0], show):
    test = yr == yv
    d = pd.DataFrame({"y": y[test], "fw": res[("dopredna", yv)]["pred_test"]},
                     index=df.index[test]).resample("D").mean()
    ax.plot(d.index, d.y.rolling(7, center=True).mean(), color=INK, lw=1.4, label="měřeno")
    ax.plot(d.index, d.fw.rolling(7, center=True).mean(), color=ORANGE, lw=1.6,
            label="dopředná predikce (úroveň a FVE drží konec tréninku)")
    if ("vynechany", yv) in res:
        d["lo"] = pd.Series(res[("vynechany", yv)]["pred_test"],
                            index=df.index[test]).resample("D").mean()
        ax.plot(d.index, d.lo.rolling(7, center=True).mean(), color=BLUE, lw=1.6,
                label="vynechaný rok (úroveň interpolovaná)")
    ax.set_ylabel("denní průměr (7denní)")
    ax.set_title(f"Testovaný rok {yv}", loc="left")
    ax.legend(frameon=False, ncol=3, fontsize=9)
fig.savefig(OUT / "24_validace_prubeh.png", dpi=130, bbox_inches="tight")
# --- 25. uroven baze a kapacita FVE: fit ze vsech dat vs. foldy ---
import pv  # noqa: E402
days = pd.date_range(df.index[0], df.index[-1], freq="D")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4))
a1.plot(days, base.level_curve(days, full[0]), color=INK, lw=2, label="vše")
a2.plot(days, pv.capacity_curve(days, full[3]) * 845, color=INK, lw=2, label="vše")
for (kind, yv), r in res.items():
    col = BLUE if kind == "vynechany" else ORANGE
    lab = f"{'vynechaný' if kind == 'vynechany' else 'dopředná'} {yv}"
    a1.plot(days, base.level_curve(days, r["params"][0]), color=col, lw=1.1, alpha=0.8, label=lab)
    a2.plot(days, pv.capacity_curve(days, r["params"][3]) * 845, color=col, lw=1.1, alpha=0.8, label=lab)
a1.set_ylabel("bazální úroveň")
a1.set_title("Úroveň báze po foldech", loc="left")
a2.set_ylabel("výkon FVE při 845 W/m²")
a2.set_title("Kapacita FVE po foldech", loc="left")
a2.legend(frameon=False, fontsize=8, ncol=2)
fig.savefig(OUT / "25_validace_uroven_fve.png", dpi=130, bbox_inches="tight")
print("\ngrafy 22-25 ulozeny")
