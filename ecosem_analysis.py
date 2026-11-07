import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.diagnostic import het_breuschpagan
from statsmodels.stats.outliers_influence import OLSInfluence, variance_inflation_factor
from statsmodels.stats.stattools import jarque_bera

FIRST, LAST, COMPARE_END = 2014, 2023, 2022
YEARS = list(range(FIRST, LAST + 1))
BREAK_YEARS = list(range(2016, 2023))
CRITICAL_T = -1.65
EXPECTED = {
    "packaging_waste": ["country", "iso2", "year", "kg_per_capita", "obs_flag"],
    "gdp_pps": ["country", "iso2", "year", "gdp_pps_eu27_100"],
    "wtp": ["country", "iso2", "wtp_2014_qa10", "wtp_2024_qb8"],
    "eu27_aggregate": ["year", "kg_per_capita"],
}

class DataError(Exception):
    pass

def numeric(frame, column, sheet):
    converted = pd.to_numeric(frame[column], errors="coerce")
    if int(converted.isna().sum()) != int(frame[column].isna().sum()):
        raise DataError(f"Non-numeric entries in {sheet}.{column}")
    return converted

def load(path):
    source = Path(path)
    if not source.is_file():
        raise DataError(f"Input file not found: {source}")
    try:
        sheets = pd.read_excel(source, sheet_name=list(EXPECTED), engine="openpyxl")
    except ValueError as exc:
        raise DataError(f"Required sheets missing: {exc}") from exc
    except Exception as exc:
        raise DataError(f"Could not read {source}: {exc}") from exc
    for name, columns in EXPECTED.items():
        absent = [c for c in columns if c not in sheets[name].columns]
        if absent:
            raise DataError(f"Sheet {name} lacks columns: {', '.join(absent)}")
    return sheets

def prepare(sheets):
    waste, gdp, wtp, agg = (sheets[k].copy() for k in EXPECTED)
    waste["year"] = numeric(waste, "year", "packaging_waste").astype(int)
    waste["kg_per_capita"] = numeric(waste, "kg_per_capita", "packaging_waste")
    gdp["year"] = numeric(gdp, "year", "gdp_pps").astype(int)
    gdp["gdp_pps_eu27_100"] = numeric(gdp, "gdp_pps_eu27_100", "gdp_pps")
    for column in ("wtp_2014_qa10", "wtp_2024_qb8"):
        wtp[column] = numeric(wtp, column, "wtp")
    agg["year"] = numeric(agg, "year", "eu27_aggregate").astype(int)
    agg["kg_per_capita"] = numeric(agg, "kg_per_capita", "eu27_aggregate")
    if waste.duplicated(["country", "year"]).any() or gdp.duplicated(["country", "year"]).any():
        raise DataError("Duplicate country-year rows")
    if wtp["country"].duplicated().any():
        raise DataError("Duplicate countries in wtp")
    countries = set(waste["country"])
    if len(countries) != 27 or countries != set(gdp["country"]) or countries != set(wtp["country"]):
        raise DataError("Sheets must hold the same 27 countries")
    if not set(YEARS) <= set(waste["year"]) or not set(YEARS) <= set(gdp["year"]):
        raise DataError(f"Years {FIRST} to {LAST} must be present")
    if (waste["kg_per_capita"].dropna() <= 0).any() or (gdp["gdp_pps_eu27_100"] <= 0).any():
        raise DataError("Waste and GDP values must be positive")
    if gdp["gdp_pps_eu27_100"].isna().any() or wtp[["wtp_2014_qa10", "wtp_2024_qb8"]].isna().any().any():
        raise DataError("GDP and willingness-to-pay values must be complete")
    if not wtp[["wtp_2014_qa10", "wtp_2024_qb8"]].apply(lambda s: s.between(0, 1)).all().all():
        raise DataError("Willingness-to-pay shares must lie between 0 and 1")
    wide = waste.pivot(index="country", columns="year", values="kg_per_capita")[YEARS]
    flags = {(r.country, r.year) for r in waste.itertuples() if isinstance(r.obs_flag, str) and r.obs_flag.strip()}
    gdp_wide = gdp.pivot(index="country", columns="year", values="gdp_pps_eu27_100")[YEARS]
    return wide, gdp_wide, wtp.set_index("country")[["wtp_2014_qa10", "wtp_2024_qb8"]], flags, agg.set_index("year")["kg_per_capita"]

def build_panel(wide, gdp, flags):
    rows = []
    for country in wide.index:
        for t in YEARS[1:]:
            earlier, later = wide.at[country, t - 1], wide.at[country, t]
            if pd.notna(earlier) and pd.notna(later):
                rows.append({"country": country, "t": t, "growth": np.log(later) - np.log(earlier), "level": np.log(earlier),
                             "log_gdp": np.log(gdp.at[country, t - 1]), "flagged": (country, t) in flags or (country, t - 1) in flags})
    panel = pd.DataFrame(rows)
    if panel.empty:
        raise DataError("No usable country-year changes")
    return panel

def build_endpoint(wide, gdp, wtp):
    rows = []
    for country in wide.index:
        series = wide.loc[country].dropna()
        if FIRST not in series.index:
            raise DataError(f"{country} has no {FIRST} value")
        end = int(series.index.max())
        if end <= FIRST:
            raise DataError(f"{country} has no later year")
        rows.append({"country": country, "end_year": end, "growth": (np.log(series[end]) - np.log(series[FIRST])) / (end - FIRST),
                     "level": np.log(series[FIRST]), "log_gdp": np.log(gdp.at[country, FIRST]),
                     "wtp": wtp.at[country, "wtp_2014_qa10"], "qb8": wtp.at[country, "wtp_2024_qb8"]})
    return pd.DataFrame(rows).set_index("country")

def ols(y, X, groups=None, cov=None):
    model = sm.OLS(y, sm.add_constant(X, has_constant="add"))
    if groups is not None:
        return model.fit(cov_type="cluster", cov_kwds={"groups": groups})
    return model.fit(cov_type=cov) if cov else model.fit()

def speed(beta):
    if 1 + beta <= 0:
        return np.nan, np.nan
    lam = -np.log(1 + beta)
    return lam, (np.log(2) / lam if lam > 0 else np.inf)

def year_dummies(panel):
    return pd.get_dummies(panel["t"], prefix="y", drop_first=True, dtype=float)

def country_dummies(panel):
    return pd.get_dummies(panel["country"], prefix="c", drop_first=True, dtype=float)

def beta_models(panel):
    g, groups = panel["growth"], panel["country"]
    years, countries = year_dummies(panel), country_dummies(panel)
    specs = {
        "eq1_absolute": panel[["level"]],
        "eq2_income_year_effects": pd.concat([panel[["level", "log_gdp"]], years], axis=1),
        "eq2_income_no_year_effects": panel[["level", "log_gdp"]],
        "two_way_fixed_effects": pd.concat([panel[["level", "log_gdp"]], countries, years], axis=1),
        "country_fixed_effects": pd.concat([panel[["level"]], countries], axis=1),
    }
    rows = []
    for name, X in specs.items():
        fit = ols(g, X, groups)
        lam, half = speed(fit.params["level"])
        low, high = fit.conf_int().loc["level"]
        row = {"model": name, "beta": fit.params["level"], "se": fit.bse["level"], "p": fit.pvalues["level"], "ci_low": low,
               "ci_high": high, "lambda": lam, "half_life_years": half, "r2": fit.rsquared, "n": int(fit.nobs)}
        if "log_gdp" in X:
            row.update({"gamma": fit.params["log_gdp"], "gamma_se": fit.bse["log_gdp"], "gamma_p": fit.pvalues["log_gdp"]})
        rows.append(row)
    return pd.DataFrame(rows)

def nested_f(restricted, unrestricted):
    q = unrestricted.df_model - restricted.df_model
    f_value = ((restricted.ssr - unrestricted.ssr) / q) / (unrestricted.ssr / unrestricted.df_resid)
    return f_value, stats.f.sf(f_value, q, unrestricted.df_resid), int(q), int(unrestricted.df_resid)

def cross_section(endpoint):
    cols = {"A_level": ["level"], "B_income": ["level", "log_gdp"], "C_income_wtp2014": ["level", "log_gdp", "wtp"], "D_income_qb8": ["level", "log_gdp", "qb8"]}
    fits = {name: ols(endpoint["growth"], endpoint[c]) for name, c in cols.items()}
    rows = []
    for name, fit in fits.items():
        for term in fit.params.index:
            rows.append({"model": name, "term": term, "coef": fit.params[term], "se": fit.bse[term], "p": fit.pvalues[term],
                         "r2": fit.rsquared, "adj_r2": fit.rsquared_adj, "n": int(fit.nobs)})
    tests = []
    for label, r, u in [("income", "A_level", "B_income"), ("wtp2014_given_income", "B_income", "C_income_wtp2014"),
                        ("income_and_wtp2014", "A_level", "C_income_wtp2014"), ("qb8_given_income", "B_income", "D_income_qb8")]:
        f_value, p_value, q, df = nested_f(fits[r], fits[u])
        tests.append({"test": label, "F": f_value, "p": p_value, "df1": q, "df2": df})
    ci = fits["C_income_wtp2014"].conf_int().loc["wtp"]
    extra = [{"measure": "wtp2014_coef_ci_low", "value": ci[0]}, {"measure": "wtp2014_coef_ci_high", "value": ci[1]},
             {"measure": "wtp2014_p_hc3", "value": ols(endpoint["growth"], endpoint[cols["C_income_wtp2014"]], cov="HC3").pvalues["wtp"]},
             {"measure": "corr_wtp2014_qb8", "value": stats.pearsonr(endpoint["wtp"], endpoint["qb8"])[0]},
             {"measure": "spearman_wtp2014_qb8", "value": stats.spearmanr(endpoint["wtp"], endpoint["qb8"])[0]},
             {"measure": "corr_wtp2014_income", "value": stats.pearsonr(endpoint["wtp"], endpoint["log_gdp"])[0]},
             {"measure": "corr_wtp2014_initial_level", "value": stats.pearsonr(endpoint["wtp"], endpoint["level"])[0]}]
    for dropped in ("Estonia", "Luxembourg"):
        if dropped in endpoint.index:
            sub = endpoint.drop(dropped)
            fit = ols(sub["growth"], sub[cols["C_income_wtp2014"]])
            extra += [{"measure": f"drop_{dropped}_beta1", "value": fit.params["level"]}, {"measure": f"drop_{dropped}_beta1_p", "value": fit.pvalues["level"]},
                      {"measure": f"drop_{dropped}_wtp_p", "value": fit.pvalues["wtp"]}]
    return pd.DataFrame(rows), pd.DataFrame(tests), pd.DataFrame(extra)

def trend(values):
    x = np.arange(len(values), dtype=float)
    fit = sm.OLS(np.asarray(values, dtype=float), sm.add_constant(x)).fit()
    return float(fit.params[1]), float(fit.bse[1]), float(fit.pvalues[1]), float(fit.rsquared)

def dispersion(wide):
    rows = []
    for year in YEARS:
        values = wide[year].dropna()
        rows.append({"year": year, "countries": len(values), "mean_kg": values.mean(), "sd_kg": values.std(ddof=0),
                     "cv": values.std(ddof=0) / values.mean(), "sd_ln": np.log(values).std(ddof=0)})
    yearly = pd.DataFrame(rows)
    balanced = wide.dropna()
    bal = pd.DataFrame({"year": YEARS, "countries": len(balanced), "mean_kg": balanced.mean().values,
                        "cv": (balanced.std(ddof=0) / balanced.mean()).values, "sd_ln": np.log(balanced).std(ddof=0).values})
    trends = []
    for label, frame in (("all_available", yearly), ("balanced", bal)):
        for column in ("cv", "sd_ln", "mean_kg"):
            slope, se, p_value, r2 = trend(frame[column])
            trends.append({"sample": label, "measure": column, "slope_per_year": slope, "se": se, "p": p_value, "r2": r2, "countries_last_year": int(frame["countries"].iloc[-1])})
    return yearly, bal, pd.DataFrame(trends)

def direction(wide, agg):
    both = wide[[FIRST, COMPARE_END]].dropna()
    first, last = both[FIRST], both[COMPARE_END]
    m0, m1 = first.mean(), last.mean()
    summary = pd.DataFrame([{"window": f"{FIRST}-{COMPARE_END}", "countries": len(both), "mean_start": m0, "mean_end": m1, "mean_change_pct": 100 * (m1 / m0 - 1),
                             "sd_ln_start": np.log(first).std(ddof=0), "sd_ln_end": np.log(last).std(ddof=0),
                             "cv_start": first.std(ddof=0) / m0, "cv_end": last.std(ddof=0) / m1}])
    gap_start, gap_end = (np.log(first) - np.log(m0)).abs(), (np.log(last) - np.log(m1)).abs()
    countries = pd.DataFrame({"start_kg": first, "end_kg": last, "change_pct": 100 * (last / first - 1),
                              "gap_narrowed": gap_end < gap_start, "own_waste_fell": last < first})
    cells = pd.crosstab(countries["gap_narrowed"], countries["own_waste_fell"]).reindex(index=[True, False], columns=[True, False], fill_value=0)
    rho, rho_p = stats.spearmanr(first, countries["change_pct"])
    aggregate = pd.DataFrame({"year": agg.index, "eu27_kg_per_capita": agg.values})
    return summary, countries.reset_index(), cells.reset_index(), pd.DataFrame([{"spearman_initial_level_vs_change": rho, "p": rho_p}]), aggregate

def log_t(block, start):
    h = block / block.mean()
    dispersion_h = ((h - 1) ** 2).mean()
    t = np.arange(1, block.shape[1] + 1)
    y = np.log(dispersion_h.iloc[0] / dispersion_h.values) - 2 * np.log(np.log(t + 1))
    x = np.log(t)
    fit = sm.OLS(y[start:], sm.add_constant(x[start:])).fit(cov_type="HAC", cov_kwds={"maxlags": 1})
    return float(fit.params[1]), float(fit.tvalues[1]), dispersion_h

def find_clubs(balanced, start):
    order = balanced.iloc[:, -3:].mean(axis=1).sort_values(ascending=False).index.tolist()
    block = balanced.loc[order]
    remaining, clubs = list(order), []
    while len(remaining) > 1:
        t_stats = {k: log_t(block.loc[remaining[:k]], start)[1] for k in range(2, len(remaining) + 1)}
        eligible = [k for k, v in t_stats.items() if v > CRITICAL_T]
        if not eligible:
            clubs.append([remaining[0]])
            remaining = remaining[1:]
            continue
        core = remaining[:max(eligible, key=t_stats.get)]
        for country in remaining[len(core):]:
            if log_t(block.loc[core + [country]], start)[1] > CRITICAL_T:
                core.append(country)
        clubs.append(core)
        remaining = [c for c in remaining if c not in core]
    clubs.extend([[c] for c in remaining])
    return clubs

def phillips_sul(wide):
    balanced = wide.dropna()
    if len(balanced) < 5:
        raise DataError("Too few balanced countries for the log-t test")
    tests, members = [], []
    for label, start in (("full_ten_years", 0), ("drop_first_two_years", 2)):
        b, t_value, h = log_t(balanced, start)
        tests.append({"window": label, "countries": len(balanced), "b": b, "t": t_value, "one_sided_p": float(stats.norm.cdf(t_value)),
                      "reject_convergence": t_value < CRITICAL_T, "h_first": h.iloc[0], "h_last": h.iloc[-1]})
        for number, club in enumerate(find_clubs(balanced, start), start=1):
            club_b, club_t = (log_t(balanced.loc[club], start)[:2] if len(club) > 1 else (np.nan, np.nan))
            members.append({"window": label, "club": number, "size": len(club), "club_b": club_b, "club_t": club_t, "members": ", ".join(club)})
    return pd.DataFrame(tests), pd.DataFrame(members)

def ssr(y, X):
    return float(sm.OLS(y, sm.add_constant(X)).fit().ssr)

def timing(panel):
    g, groups = panel["growth"], panel["country"]
    years, n, k = year_dummies(panel), len(panel), 2
    rows = []
    for year in BREAK_YEARS:
        early = panel["t"] < year
        if early.sum() <= k or (~early).sum() <= k:
            continue
        rss_pooled = ssr(g, panel[["level"]])
        rss_1, rss_2 = ssr(g[early], panel.loc[early, ["level"]]), ssr(g[~early], panel.loc[~early, ["level"]])
        f_value = ((rss_pooled - rss_1 - rss_2) / k) / ((rss_1 + rss_2) / (n - 2 * k))
        design = pd.concat([panel[["level"]].assign(interaction=panel["level"] * (~early).astype(float)), years], axis=1)
        fit = ols(g, design, groups)
        rows.append({"break_year": year, "chow_F": f_value, "chow_p": stats.f.sf(f_value, k, n - 2 * k), "theta": fit.params["interaction"],
                     "theta_p": fit.pvalues["interaction"], "beta_before": sm.OLS(g[early], sm.add_constant(panel.loc[early, "level"])).fit().params["level"],
                     "beta_after": sm.OLS(g[~early], sm.add_constant(panel.loc[~early, "level"])).fit().params["level"],
                     "n_before": int(early.sum()), "n_after": int((~early).sum()), "bonferroni_alpha": 0.05 / len(BREAK_YEARS),
                     "mean_growth_before": g[early].mean(), "mean_growth_after": g[~early].mean()})
    return pd.DataFrame(rows)

def diagnostics(panel, endpoint):
    g = panel["growth"]
    design = pd.concat([panel[["level", "log_gdp"]], year_dummies(panel)], axis=1)
    exog = sm.add_constant(design)
    fit = sm.OLS(g, exog).fit()
    bp_stat, bp_p = het_breuschpagan(fit.resid, exog)[:2]
    jb_stat, jb_p, skew, kurt = jarque_bera(fit.resid)
    core = sm.add_constant(panel[["level", "log_gdp"]])
    standardised = (panel[["level", "log_gdp"]] - panel[["level", "log_gdp"]].mean()) / panel[["level", "log_gdp"]].std()
    cooks = OLSInfluence(fit).cooks_distance[0]
    worst = int(np.argmax(cooks))
    leave_out = {c: sm.OLS(g[panel["country"] != c], sm.add_constant(panel.loc[panel["country"] != c, "level"])).fit().params["level"] for c in panel["country"].unique()}
    rows = [("breusch_pagan_stat", bp_stat), ("breusch_pagan_p", bp_p), ("jarque_bera_stat", jb_stat), ("jarque_bera_p", jb_p), ("skewness", skew), ("kurtosis", kurt),
            ("vif_level", variance_inflation_factor(core.values, 1)), ("vif_log_gdp", variance_inflation_factor(core.values, 2)),
            ("condition_number_standardised", np.linalg.cond(sm.add_constant(standardised).values)), ("max_cooks_distance", cooks[worst]),
            ("max_cooks_obs", f"{panel['country'].iloc[worst]} {panel['t'].iloc[worst]}"), ("cooks_threshold_4_over_n", 4 / len(panel)),
            ("obs_above_threshold", int((cooks > 4 / len(panel)).sum())), ("leave_one_country_out_beta_min", min(leave_out.values())),
            ("leave_one_country_out_beta_max", max(leave_out.values()))]
    ecols = ["level", "log_gdp", "wtp"]
    efit = sm.OLS(endpoint["growth"], sm.add_constant(endpoint[ecols])).fit()
    ex = sm.add_constant(endpoint[ecols])
    ecooks = OLSInfluence(efit).cooks_distance[0]
    estd = (endpoint[ecols] - endpoint[ecols].mean()) / endpoint[ecols].std()
    rows += [("endpoint_breusch_pagan_p", het_breuschpagan(efit.resid, ex)[1]), ("endpoint_jarque_bera_p", jarque_bera(efit.resid)[1]),
             ("endpoint_vif_level", variance_inflation_factor(ex.values, 1)), ("endpoint_vif_log_gdp", variance_inflation_factor(ex.values, 2)),
             ("endpoint_vif_wtp", variance_inflation_factor(ex.values, 3)), ("endpoint_condition_number_standardised", np.linalg.cond(sm.add_constant(estd).values)),
             ("endpoint_max_cooks_distance", ecooks.max()), ("endpoint_max_cooks_country", endpoint.index[int(np.argmax(ecooks))]),
             ("endpoint_cooks_threshold_4_over_n", 4 / len(endpoint))]
    return pd.DataFrame(rows, columns=["diagnostic", "value"])

def flagged_check(panel):
    kept = panel[~panel["flagged"]]
    years = year_dummies(kept)
    absolute = ols(kept["growth"], kept[["level"]], kept["country"])
    conditional = ols(kept["growth"], pd.concat([kept[["level", "log_gdp"]], years], axis=1), kept["country"])
    return pd.DataFrame([{"dropped_transitions": int(panel["flagged"].sum()), "n": len(kept), "beta_eq1": absolute.params["level"], "p_eq1": absolute.pvalues["level"],
                          "beta_eq2": conditional.params["level"], "p_eq2": conditional.pvalues["level"], "gamma_p": conditional.pvalues["log_gdp"]}])

def emit(name, frame, out_dir):
    print(f"\n{name}")
    print(frame.round(4).to_string(index=False))
    frame.to_csv(out_dir / f"{name}.csv", index=False)

def parse_args(argv):
    parser = argparse.ArgumentParser(description="EU packaging waste convergence analysis")
    parser.add_argument("--data", default="ecosem_data.xlsx")
    parser.add_argument("--out", default="results")
    return parser.parse_args(argv)

def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    out_dir = Path(args.out)
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        wide, gdp, wtp, flags, agg = prepare(load(args.data))
        panel = build_panel(wide, gdp, flags)
        endpoint = build_endpoint(wide, gdp, wtp)
    except DataError as exc:
        print(f"Data error: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"File system error: {exc}", file=sys.stderr)
        return 1
    print(f"Panel: {len(panel)} country-year changes, {panel['country'].nunique()} countries; cross-section: {len(endpoint)} countries")
    stages = [
        ("beta_convergence", lambda: [("beta_convergence", beta_models(panel))]),
        ("cross_section", lambda: list(zip(("cross_section_models", "cross_section_f_tests", "cross_section_extras"), cross_section(endpoint)))),
        ("sigma", lambda: list(zip(("dispersion_all_available", "dispersion_balanced", "dispersion_trends"), dispersion(wide)))),
        ("direction", lambda: list(zip(("direction_summary", "direction_countries", "direction_cells", "direction_spearman", "eu27_aggregate"), direction(wide, agg)))),
        ("phillips_sul", lambda: list(zip(("phillips_sul_tests", "phillips_sul_clubs"), phillips_sul(wide)))),
        ("timing", lambda: [("timing_tests", timing(panel))]),
        ("diagnostics", lambda: [("diagnostics", diagnostics(panel, endpoint))]),
        ("flagged_check", lambda: [("flagged_check", flagged_check(panel))]),
    ]
    failed = []
    for label, runner in stages:
        try:
            for name, frame in runner():
                emit(name, frame, out_dir)
        except (DataError, ValueError, KeyError, IndexError, np.linalg.LinAlgError, ZeroDivisionError) as exc:
            failed.append(label)
            print(f"Stage {label} failed: {type(exc).__name__}: {exc}", file=sys.stderr)
    if failed:
        print(f"Completed with failures in: {', '.join(failed)}", file=sys.stderr)
        return 2
    print(f"\nResults written to {out_dir.resolve()}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
