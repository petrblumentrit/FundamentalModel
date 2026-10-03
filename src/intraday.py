"""Intraday predikce — klouzavy backtest s vydanim po kazdem 15min mereni.

Fundamentalni model nema autoregresni clen: predikce na cas t zavisi jen na
parametrech, kalendari a pocasi. Nove mereni se proto do predikce na nejblizsi
hodiny dostava dvema cestami:

1. **prepocet s cerstvym pocasim** — do konce dat s namerene pocasi (stav
   setrvacnostnich filtru), od s dal predpoved. Archiv predpovedi ma jedinou
   predpoved na den (z D 10:00 na D+1), takze zpresneni predpovedi s kratsim
   predstihem simulovat nejde — to ukaze az provoz;
2. **korekce zavisla na horizontu** — reziduova vrstva nad modelem, aditivni:

       korekce(s, h) = a_h*x1 + b_h*xden + c_h*vcera(h) + d_h*tyden(h)

   z rezidui r(u) = skutecnost - model pro u < s (model s parametry platnymi
   v s a namerenym pocasim; vse zname v okamziku vydani):

   - x1    = prumer r za posledni hodinu,
   - xden  = prumer r za poslednich 24 h,
   - vcera = r ve stejnem case dne jako cil, posledni znamy den (prumer hodiny
             koncici timto casem),
   - tyden = totez, prumer pres poslednich TOD_DAYS znamych dni.

   S predpovedi pocasi navic stejna ctverice z w(u) = model s namerenym
   pocasim - model s predpovedi pocasi (x1_w, xden_w, vcera_w, tyden_w): chyba
   predpovedi pocasi prepoctena modelem na spotrebu. Je setrvacna (predpoved
   o 2 C teplejsi nez stanice zustane teplejsi i dalsi hodiny) a v r neni —
   r je pocitane z namereneho pocasi. Bez techto clenu korekce na 15 min
   snizila RMSE jen 16,7 -> 12,0, protoze chybu pocasi v cili nevidela.

   Vaha posledni hodiny s horizontem klesa, vaha tvaru podle casu dne roste —
   koeficienty se proto odhaduji pro kazdy horizont zvlast, **online**: jen
   z vydani, jejichz cil byl v okamziku s uz zmereny.

Konvence: s = konec dostupnych dat (prvni nezmereny interval zacina v s);
horizont h = 0 je interval [s, s+15 min), tj. konec ciloveho intervalu je
(h+1)*15 min za koncem dat. Skutecny cas vydani = s + zpozdeni dat; pri
zpozdeni mereni o L se vysledky ctou na horizontu posunutem o L.

Parametry: linearni prefit jednou denne z dat do D 09:00 (backtest.refit),
plati pro s v [D 09:00, D+1 09:00); tvar tydne. Cile konci koncem dne D+1
(dal predpoved pocasi v archivu v okamziku s neni).
"""
import numpy as np
import pandas as pd

import backtest
import config
import correction
import etl
import meteo_forecast
import sun

HORIZON = 156      # [intervaly] 39 h: od D 09:00 do konce D+1
DAY = 96           # intervalu za 24 h (posun o den v UTC; kolem zmeny casu o hodinu vedle)
_C = config.model()["korekce"]["intraday"]   # hodnoty v config/model.yaml
LAST = _C["posledni_intervaly"]   # [intervaly] okno "posledni hodina"
TOD_DAYS = _C["tvar_dni"]         # [dny] okno tvaru chyby podle casu dne
FEATURES = ("x1", "xden", "vcera", "tyden", "x1_w", "xden_w", "vcera_w", "tyden_w")
W_DAYS = TOD_DAYS + 2   # [dny] historie, na ktere se pocita model s predpovedi pocasi
WARMUP = _C["rozjezd_dni"] * DAY  # [intervaly] doba od zacatku simulace, nez se korekce zacne pouzivat
MIN_ROWS = 21      # minimalni pocet znamych vydani na danem horizontu
RIDGE = _C["ridge"]       # relativne ke stope X'X
METEO_COLS = ("temp", "sun", "wind", "sero")


def _forecast_columns(recent: pd.DataFrame, df: pd.DataFrame, meteo: pd.DataFrame,
                      start: pd.Timestamp) -> dict:
    """Predpoved pocasi od `start` na ose `recent` (jako backtest.forecast_day):
    teplota a vitr bez klouzaveho biasu, sero jako stredni hodnota pres
    rozdeleni jasnosti. Pred `start` a kde predpoved chybi, zustava mereni."""
    fut = recent.index >= start
    fc = meteo.reindex(recent.index[fut])
    if backtest.DEBIAS:
        adj = meteo_forecast.debias(meteo, df, start, backtest.DEBIAS).reindex(fc.index)
        fc[list(backtest.DEBIAS)] = adj[list(backtest.DEBIAS)]
    out = {c: recent[c].to_numpy(float).copy() for c in METEO_COLS}
    for c in ("temp", "sun", "wind"):
        out[c][fut] = fc[c].fillna(recent.loc[fut, c]).to_numpy()
    if "sero_den" in fc:
        g = (1.0 - recent.loc[fut, "tma"]) * fc["sero_den"].to_numpy()
        out["sero"][fut] = g.fillna(recent.loc[fut, "sero"]).to_numpy()
    else:
        out["sero"][fut] = sun.gloom(recent.loc[fut, "tma"], out["sun"][fut])
    return out


def inputs(v: np.ndarray, i: int, j: np.ndarray) -> dict:
    """Vstupy korekce z rady v (rezidua nebo chyba pocasi) zname pred pozici i
    (konec dat) pro cile na pozicich j >= i."""
    # stejny cas dne v poslednich znamych dnech; hodina koncici casem cile
    back = j[:, None] - DAY * (((j - i) // DAY + 1)[:, None] + np.arange(TOD_DAYS))
    tod = np.mean([v[back - q] for q in range(LAST)], axis=0)
    return {"x1": v[i - LAST:i].mean(), "xden": v[i - DAY:i].mean(),
            "vcera": tod[:, 0], "tyden": tod.mean(axis=1)}


def apply(coefs: pd.DataFrame, x: dict) -> np.ndarray:
    """Korekce z ulozenych koeficientu (radky = horizonty cilu, sloupce = cleny;
    explore/intraday.py -> koeficienty.csv) a vstupu x (inputs)."""
    return sum(coefs[f].to_numpy() * x[f] for f in coefs.columns)


def day_issues(df: pd.DataFrame, day: pd.Timestamp, shape: tuple,
               window: int | None = None, meteo: pd.DataFrame | None = None) -> dict:
    """Vsechna vydani s v [D 09:00, D+1 09:00): predikce modelu na horizonty
    0..HORIZON-1 a vstupy korekce. Pole (pocet vydani x HORIZON), NaN za koncem
    dne D+1 a za koncem dat."""
    params = backtest.refit(df, backtest.train_mask(df, day, window), shape, fix_shape=True)
    last = day + pd.Timedelta(days=1)
    dates = pd.to_datetime(df["date_local"])
    recent = df[((dates >= last - pd.Timedelta(days=backtest.PREDICT_HISTORY)) & (dates <= last)).to_numpy()].copy()
    start = backtest.cutoff(day).tz_convert("UTC")
    end = backtest.cutoff(last).tz_convert("UTC")
    n = len(recent)
    y = recent["baseload"].to_numpy(float)
    issues = np.flatnonzero((recent.index >= start) & (recent.index < end))
    measured = {c: recent[c].to_numpy(float).copy() for c in METEO_COLS}
    forecast = _forecast_columns(recent, df, meteo, start) if meteo is not None else None
    pos = np.arange(n)
    pw = None
    if meteo is not None:
        # model s predpovedi pocasi i v minulych dnech (kazdy den ma v archivu
        # jednu predpoved); bias teploty a vetru podle stavu pred W_DAYS dny
        past = _forecast_columns(recent, df, meteo, start - pd.Timedelta(days=W_DAYS))
        for c in METEO_COLS:
            recent[c] = past[c]
        pw = backtest.components(recent, params)["predikce"].to_numpy()

    shape_ = (len(issues), HORIZON)
    two, one = ("pred", "skut", "vcera", "tyden"), ("x1", "xden")
    if pw is not None:
        two, one = two + ("vcera_w", "tyden_w"), one + ("x1_w", "xden_w")
    out = {k: np.full(shape_, np.nan, np.float32) for k in two}
    out.update({k: np.full(len(issues), np.nan, np.float32) for k in one})
    out["s"] = recent.index[issues].tz_convert(None).to_numpy("datetime64[s]")
    p = None
    for m, i in enumerate(issues):
        if forecast is not None:
            for c in METEO_COLS:
                recent[c] = np.where(pos < i, measured[c], forecast[c])
        if forecast is not None or p is None:    # s namerenym pocasim predikce na s nezavisi
            p = backtest.components(recent, params)["predikce"].to_numpy()
        r = y - p
        j = np.arange(i, min(i + HORIZON, n))
        h = j - i
        out["pred"][m, h] = p[j]
        out["skut"][m, h] = y[j]
        for v, suf in ((r, ""), (None if pw is None else p - pw, "_w")):
            if v is None:
                continue
            x = inputs(v, i, j)
            for k in ("x1", "xden"):
                out[k + suf][m] = x[k]
            for k in ("vcera", "tyden"):
                out[k + suf][m, h] = x[k]
    return out


def stack(parts: list[dict]) -> dict:
    """Spojeni dnu do jedne casove osy vydani (serazeno podle s)."""
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    order = np.argsort(out["s"])
    out = {k: v[order] for k, v in out.items()}
    assert (np.diff(out["s"]) == np.timedelta64(15, "m")).all(), "vydani nejsou po 15 min"
    return out


def available(res: dict, feats: tuple[str, ...] = FEATURES) -> tuple[str, ...]:
    """Cleny korekce pritomne ve vysledcich (cleny _w jen s predpovedi pocasi)."""
    return tuple(f for f in feats if f in res)


def online(res: dict, feats: tuple[str, ...] = FEATURES) -> tuple[np.ndarray, np.ndarray]:
    """Korekce (vydani x horizont) a koeficienty (vydani x horizont x clen).

    Pro vydani i a horizont h se koeficienty odhaduji LS (ridge ke nule) z
    vydani i' <= i - h - 1 — jejich cil na horizontu h byl v s uz zmereny.
    Prvnich WARMUP intervalu je korekce nulova."""
    e = res["skut"].astype(float) - res["pred"].astype(float)
    n, nh = e.shape
    feats = available(res, feats)
    k = len(feats)
    corr = np.zeros((n, nh))
    coefs = np.full((n, nh, k), np.nan)
    eye = np.eye(k)
    for h in range(nh):
        X = np.column_stack([res[f][:, h] if res[f].ndim == 2 else res[f] for f in feats]).astype(float)
        ok = ~np.isnan(X).any(axis=1) & ~np.isnan(e[:, h])
        Xz, ez = np.where(ok[:, None], X, 0.0), np.where(ok, e[:, h], 0.0)
        G = np.cumsum(Xz[:, :, None] * Xz[:, None, :], axis=0)
        b = np.cumsum(Xz * ez[:, None], axis=0)
        cnt = np.cumsum(ok)
        use = np.arange(n) - h - 1            # posledni vydani se znamym cilem
        # dlouhe horizonty ma jen cast vydani (cile konci s dnem D+1), proto
        # podminka na uplynuly cas a zvlast na pocet znamych vydani
        act = (use >= WARMUP) & ok
        act[act] = cnt[use[act]] >= MIN_ROWS
        if not act.any():
            continue
        Ga, ba = G[use[act]], b[use[act]]
        Ga = Ga + (RIDGE * np.trace(Ga, axis1=1, axis2=2) / k)[:, None, None] * eye
        c = np.linalg.solve(Ga, ba[:, :, None])[:, :, 0]
        coefs[act, h] = c
        corr[act, h] = (X[act] * c).sum(axis=1)
    return corr, coefs


def by_horizon(res: dict, corr: np.ndarray, rows: np.ndarray | None = None) -> pd.DataFrame:
    """RMSE / MAE / bias / MAPE po horizontech, model a model + korekce."""
    y, f = res["skut"].astype(float), res["pred"].astype(float)
    if rows is not None:
        y, f, corr = y[rows], f[rows], corr[rows]
    e, ek = y - f, y - f - corr
    out = pd.DataFrame({
        "horizont_h": (np.arange(y.shape[1]) + 1) / 4,
        "n": (~np.isnan(e)).sum(axis=0),
        "rmse_model": np.sqrt(np.nanmean(e**2, axis=0)),
        "rmse_kor": np.sqrt(np.nanmean(ek**2, axis=0)),
        "mae_model": np.nanmean(np.abs(e), axis=0),
        "mae_kor": np.nanmean(np.abs(ek), axis=0),
        "bias_model": np.nanmean(e, axis=0),
        "bias_kor": np.nanmean(ek, axis=0),
        "mape_model": np.nanmean(np.abs(e) / y, axis=0) * 100,
        "mape_kor": np.nanmean(np.abs(ek) / y, axis=0) * 100,
    })
    return out.set_index("horizont_h")


def local_time(res: dict) -> pd.DatetimeIndex:
    """Konec dat s v mistnim case."""
    return pd.DatetimeIndex(res["s"]).tz_localize("UTC").tz_convert(etl.TZ)


def irregular_targets(res: dict) -> np.ndarray:
    """Maska (vydani x horizont): cil pada na svatek, most nebo Vanoce."""
    t = pd.DatetimeIndex(res["s"]).tz_localize("UTC")
    days = pd.date_range(t[0].tz_convert(etl.TZ).normalize().tz_localize(None) - pd.Timedelta(days=1),
                         periods=len(t) // DAY + 5, freq="D")
    irr = pd.Series(correction.irregular_days(days), index=days)
    out = np.zeros(res["pred"].shape, bool)
    for h in range(out.shape[1]):
        d = (t + pd.Timedelta(minutes=15 * h)).tz_convert(etl.TZ).normalize().tz_localize(None)
        out[:, h] = irr.reindex(d).fillna(False).to_numpy()
    return out
