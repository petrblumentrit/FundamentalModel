"""Zatizeni CR (ENTSO-E / CEPS) vs portfolio: chovani a vyuziti predpovedi CEPS.

A) Fundamentalni model napasovany na zatizeni CR (skutecne zatizeni ENTSO-E,
   hodinove -> 15 min) se stejnym meteo (jedna stanice — pro celou CR jen
   priblizne). Srovnani s portfoliem v relativnich jednotkach: topna citlivost
   a jeji trend, uroven (krize 2022), FVE (decentralni vyroba snizuje merene
   zatizeni), osvetleni.
B) Pomuze denni predpoved CEPS zpresnit predikci portfolia?
   1. sdilene vykyvy: korelace nevysvetlenych odchylek CR a portfolia
   2. "prekvapeni CEPS" x4 = predpoved CEPS / model CR - 1 (model CR fitovany
      jen z dat pred backtestem — bez uniku z budoucnosti) jako dalsi vstup
      korekcni vrstvy; online odhad jako u ostatnich vstupu korekce.
   Pozor: backtest pouziva skutecne meteo, predpoved CEPS predpovezene — jeji
   prinos v provozu bude spis vyssi. Dostupnost predpovedi v 10:00 overit.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import backtest
import base
import correction
import cooling
import etl
import fit
import heating
import ote
import pv
import validate

ROOT = Path(__file__).resolve().parent.parent
# --bt=slozka: backtest k vyhodnoceni (vychozi simulace/); --meteo=predpoved:
# model CR pro "prekvapeni CEPS" pocita s predpovedi pocasi (jinak by do testu
# realistickeho backtestu unikla informace o skutecnem pocasi)
BT_DIR = ROOT / "simulace"
METEO_FC = False
for _a in sys.argv:
    if _a.startswith("--bt="):
        BT_DIR = ROOT / _a.split("=", 1)[1]
    if _a == "--meteo=predpoved":
        METEO_FC = True
rm = lambda x: float(np.sqrt(np.nanmean(np.square(x))))

df = etl.load()
ent = ote.entsoe_load()
# 15min osa: hodinovy prumer opakovany ve ctvrthodinach (2025+ jsou to
# prumery skutecnych ctvrthodin — pro srovnani tvaru a urovni staci)
cz = ent["actual"].reindex(df.index.floor("h")).to_numpy()
ok = ~np.isnan(cz)
dfc = df[ok].copy()
dfc["baseload"] = cz[ok]
print(f"CR: {dfc.index[0].date()} - {dfc.index[-1].date()}, prumer {dfc.baseload.mean():.0f} MW; "
      f"portfolio prumer {df.baseload.mean():.0f}")

# --- A) model na zatizeni CR (vsechna data) vs portfolio ----------------------
full_p = (base.load(base.MODEL_DIR / "base_joint.npz"), heating.load(heating.MODEL_DIR / "heating_joint.npz"),
          fit.load("cooling_joint.npz", cooling.PARAM_NAMES + ["rmse", "r2"]),
          fit.load("pv_joint.npz", ["gamma", "span_days", "knot_days", "rmse", "r2"]))
full_c = validate.fit_masked(dfc, np.ones(len(dfc), bool))


def describe(params, d: pd.DataFrame, label: str) -> dict:
    bp, hp, cp, pp = params
    mean = d.baseload.mean()
    loc = d.index.tz_convert(etl.TZ)
    H = heating.predict(d, hp); C = cooling.predict(d, cp); G = pv.predict(d, pp)
    win = np.isin(loc.month, [12, 1, 2]); summ = np.isin(loc.month, [6, 7, 8])
    k_mean = float(np.mean(heating.k_basis(d) @ hp["k_coef"]))
    tr = 100 * heating.trend_profile(d.index[-1], hp, np.array([2.0, 12.0, 18.0])) / k_mean
    days = pd.date_range(d.index[0], d.index[-1], freq="D")
    lev = pd.Series(base.level_curve(days, bp), index=days.tz_convert(etl.TZ)).groupby(lambda x: x.year).mean()
    cap = pv.capacity_curve(pd.DatetimeIndex([d.index[0], d.index[-1]]), pp) * 845
    light = base.light_curve(bp, np.array([17.0, 20.0]))
    resid = d.baseload.to_numpy() - validate.predict(d, params)
    out = {
        "RMSE % prumeru": 100 * rm(resid) / mean,
        "T_b [C]": hp["t_b"], "s [C]": hp["s"], "tau [h]": hp["tau"], "w": hp["w"],
        "topeni zima % zatizeni": 100 * H[win].mean() / d.baseload[win].mean(),
        "chlazeni leto % zatizeni": 100 * C[summ].mean() / d.baseload[summ].mean(),
        "trend topeni noc/poledne/18h %": "/".join(f"{x:+.0f}" for x in tr),
        "FVE spicka zacatek->konec % prumeru": f"{100 * cap[0] / mean:.1f} -> {100 * cap[1] / mean:.1f}",
        "osvetleni 17h/20h % prumeru": "/".join(f"{100 * x / mean:.1f}" for x in light),
        "uroven po letech (2022=100)": " ".join(f"{int(y)}:{100 * v / lev.iloc[0]:.1f}" for y, v in lev.items()),
    }
    return out


A = pd.DataFrame({"portfolio": describe(full_p, df, "portfolio"), "CR (ENTSO-E)": describe(full_c, dfc, "CR")})
pd.set_option("display.width", 220, "display.max_colwidth", 80)
print("\nA) model portfolia vs model zatizeni CR (relativne):")
print(A.to_string())

# --- B1) sdilene vykyvy (in-sample rezidua, bez svatku a Vanoc) ----------------
e_p = pd.Series(df.baseload.to_numpy() - validate.predict(df, full_p), index=df.index)
e_c = pd.Series(dfc.baseload.to_numpy() - validate.predict(dfc, full_c), index=dfc.index)
cal = df[["date_local", "daytype", "most"]].copy()
d_loc = pd.to_datetime(cal.date_local)
hol = etl._cz_holidays(2021, 2027)
special = d_loc.dt.date.isin(hol) | cal.most | ((d_loc.dt.month == 12) & (d_loc.dt.day >= 20)) | ((d_loc.dt.month == 1) & (d_loc.dt.day <= 6))
j = pd.DataFrame({"p": e_p / df.baseload.mean(), "c": e_c / dfc.baseload.mean()}).join(cal.date_local).dropna()
j = j[~special.reindex(j.index).to_numpy()]
dly = j.groupby("date_local")[["p", "c"]].mean()
print(f"\nB1) korelace nevysvetlenych odchylek portfolio vs CR (mimo svatky a Vanoce): "
      f"15min {j.p.corr(j.c):.2f}, denni prumery {dly.p.corr(dly.c):.2f}")

# --- B2) prekvapeni CEPS jako vstup korekce ------------------------------------
bt = pd.read_csv(BT_DIR / "predikce_D1.csv", sep=";", decimal=",")
bt.index = pd.to_datetime(bt.timestamp_utc, utc=True)
start = bt.vydano.iloc[0]
mask_c = (dfc.index < backtest.cutoff(pd.Timestamp(str(start)[:10])).tz_convert("UTC"))
pre_c = validate.fit_masked(dfc, mask_c)            # model CR jen z dat pred backtestem
dfc_pred = dfc
if METEO_FC:
    import meteo_forecast as mf
    until = backtest.cutoff(pd.Timestamp(str(start)[:10])).tz_convert("UTC")
    arch = mf.load_archive()
    fc = mf.to_15min(arch, dfc.index[(dfc.index >= until) & (dfc.index <= arch.index.max())], mf.calibrate(arch, df, until))
    dfc_pred = dfc.copy()
    for c in ("temp", "sun", "wind"):
        dfc_pred.loc[fc.index, c] = fc[c].to_numpy()
    print("   model CR pro prekvapeni CEPS: predpoved pocasi od", until)
cz_model = pd.Series(validate.predict(dfc_pred, pre_c), index=dfc.index)
da = ent["da"].reindex(bt.index.floor("h")).to_numpy()
x4 = da / cz_model.reindex(bt.index).to_numpy() - 1
print(f"\n   prekvapeni CEPS x4 v backtestu: prumer {np.nanmean(x4) * 100:+.2f} %, sm. odch. {np.nanstd(x4) * 100:.2f} %")

e_bt = (bt.skutecnost - bt.predikce).to_numpy()
c3, _ = correction.online(bt)
ek = e_bt - c3
X = pd.DataFrame({"ek": ek, "x4p": x4 * bt.predikce.to_numpy()}, index=bt.index)
X["den"] = pd.to_datetime(pd.to_datetime(bt.timestamp_mistni, utc=True).dt.tz_convert(etl.TZ).dt.date).to_numpy()
ok4 = ~np.isnan(X.x4p)
print(f"   korelace x4 * predikce s chybou po korekci: 15min {X.ek[ok4].corr(X.x4p[ok4]):.2f}, "
      f"denni {X[ok4].groupby('den')[['ek', 'x4p']].mean().corr().iloc[0, 1]:.2f}")

# online odhad koeficientu u x4 (jen z minulosti, jako ostatni vstupy korekce)
days = sorted(X.den.unique())
add = np.zeros(len(X))
coefs = []
for i, dd in enumerate(days):
    known = ok4 & (X.den < dd - pd.Timedelta(days=1)).to_numpy()
    tgt = ok4 & (X.den == dd).to_numpy()
    if known.sum() < 21 * 96 or not tgt.any():
        continue
    a = X.x4p[known].to_numpy(); b = X.ek[known].to_numpy()
    k = float(a @ b / (a @ a + 1e-3 * (a @ a)))
    add[tgt] = k * X.x4p[tgt].to_numpy()
    coefs.append(k)
e5 = ek - add
t = pd.to_datetime(bt.timestamp_mistni, utc=True).dt.tz_convert(etl.TZ)
xm = (((t.dt.month == 12) & (t.dt.day >= 21)) | ((t.dt.month == 1) & (t.dt.day <= 3))).to_numpy()
y = bt.skutecnost.to_numpy()
mape = lambda e: float(np.mean(np.abs(e) / y) * 100)
print(f"   koeficient u x4 (posledni): {coefs[-1]:.2f}")
print(f"   RMSE / MAPE: model {rm(e_bt):.2f} / {mape(e_bt):.2f} %, + korekce {rm(ek):.2f} / {mape(ek):.2f} %, "
      f"+ korekce + CEPS {rm(e5):.2f} / {mape(e5):.2f} %")
print(f"   mimo Vanoce: + korekce {rm(ek[~xm]):.2f}, + korekce + CEPS {rm(e5[~xm]):.2f}")
