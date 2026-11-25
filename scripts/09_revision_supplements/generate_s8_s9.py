from pathlib import Path
import os
import math
import pandas as pd
import numpy as np

ROOT = Path(os.environ.get("ASAH_ANALYSIS_ROOT", Path.cwd())).resolve()
OUT = Path(os.environ.get("ASAH_SUPPLEMENT_OUTPUT", ROOT / "derived_local" / "revision_supplements")).resolve()
OUT.mkdir(parents=True, exist_ok=True)


def implied_corr(ri, rs, resid, t, s):
    vals = [ri, rs, resid]
    if any(pd.isna(v) for v in vals):
        return np.nan
    ri, rs, resid = map(float, vals)
    vt = ri + t * t * rs + resid
    vs = ri + s * s * rs + resid
    if vt <= 0 or vs <= 0:
        return np.nan
    return (ri + t * s * rs) / math.sqrt(vt * vs)


feature_rows = []
summary_rows = []
pairs = {
    "D01-D04": (-3.0, 0.0),
    "D01-D07": (-3.0, 3.0),
    "D04-D07": (0.0, 3.0),
}

for layer, filename in [
    ("miRNA", "STAGE1B_FULL_MODEL_RESULTS_miRNA.csv"),
    ("GENCODE", "STAGE1B_FULL_MODEL_RESULTS_GENCODE.csv"),
]:
    src = ROOT / "revision_stage1" / "stage1b_full_refit" / filename
    df = pd.read_csv(src, low_memory=False)
    for model in ["M3", "M4"]:
        for _, row in df.iterrows():
            item = {
                "layer": layer,
                "model": model,
                "feature": row["feature"],
                "model_status": row[f"{model}_model_status"],
                "random_slope_validated": row[f"{model}_random_slope_validated"],
                "numerical_evaluable": row[f"{model}_numerical_evaluable"],
                "strict_core_robust": row["strict_core_robust"],
                "RS_boundary_any": row["any_RS_boundary"],
                "random_intercept_variance": row[f"{model}_random_intercept_variance"],
                "random_slope_variance": row[f"{model}_random_slope_variance"],
                "residual_variance": row[f"{model}_residual_variance"],
            }
            for label, (t, s) in pairs.items():
                item[f"rho_{label}"] = implied_corr(
                    row[f"{model}_random_intercept_variance"],
                    row[f"{model}_random_slope_variance"],
                    row[f"{model}_residual_variance"],
                    t,
                    s,
                )
            feature_rows.append(item)

feature = pd.DataFrame(feature_rows)
feature.to_csv(OUT / "Supplementary_Table_S8_model_implied_within_subject_correlations.csv", index=False)

for (layer, model), g in feature.groupby(["layer", "model"], sort=False):
    strata = {
        "all_fixed_analysis_features": pd.Series(True, index=g.index),
        "RS_boundary_features": g["model_status"].eq("RS_BOUNDARY"),
        "strict_core_features": g["strict_core_robust"].astype(str).str.upper().eq("TRUE"),
    }
    for stratum, mask in strata.items():
        x = g.loc[mask]
        for pair in pairs:
            v = pd.to_numeric(x[f"rho_{pair}"], errors="coerce").dropna()
            summary_rows.append({
                "layer": layer,
                "model": model,
                "stratum": stratum,
                "label_pair": pair,
                "n_features_in_stratum": len(x),
                "n_evaluable_correlations": len(v),
                "median": v.median() if len(v) else np.nan,
                "Q1": v.quantile(0.25) if len(v) else np.nan,
                "Q3": v.quantile(0.75) if len(v) else np.nan,
                "minimum": v.min() if len(v) else np.nan,
                "maximum": v.max() if len(v) else np.nan,
            })

summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / "Supplementary_Table_S8_summary.csv", index=False)

# Duplicate-correlation estimates are consensus/pooled correlations, one per model family and layer.
legacy = pd.read_csv(ROOT / "revision_stage0" / "stage0_legacy_refit_correlations.csv")
legacy.to_csv(OUT / "Supplementary_Table_S8_duplicateCorrelation_consensus.csv", index=False)

# Convert 153 paired rows to 306 sample rows. The source metric used here is the
# Atlas transcriptome-mapped read count that underlies the recovered QC criterion.
pair_df = pd.read_csv(ROOT / "revision_stage0c" / "data" / "FULL_153_pairs.csv")
long_rows = []
for _, r in pair_df.iterrows():
    for fluid in ["CSF", "plasma"]:
        long_rows.append({
            "pair_id": r["pair_id"],
            "subject": r["subject"],
            "collection_label": r["collection_label"],
            "biofluid": "CSF" if fluid == "CSF" else "Plasma",
            "GSM": r[f"{fluid}_GSM"],
            "expression_id": r[f"{fluid}_expression_id"],
            "source_QC_status": r[f"{fluid}_QC"],
            "source_reported_transcriptome_mapped_reads": r[f"{fluid}_transcriptome_reads"],
            "source_reported_reference_genome_reads": r[f"{fluid}_reference_genome_reads"],
        })
long = pd.DataFrame(long_rows)
long.to_csv(OUT / "SOURCE_REPORTED_MAPPED_READS_306_SAMPLES.csv", index=False)


def q1(x):
    return x.quantile(0.25)


def q3(x):
    return x.quantile(0.75)


def summarize(g, scope, label):
    x = pd.to_numeric(g["source_reported_transcriptome_mapped_reads"], errors="coerce")
    return {
        "scope": scope,
        "biofluid": g["biofluid"].iloc[0],
        "collection_label": label,
        "n_samples": len(g),
        "n_with_mapped_reads": x.notna().sum(),
        "mapped_reads_median": x.median(),
        "mapped_reads_Q1": q1(x),
        "mapped_reads_Q3": q3(x),
        "mapped_reads_IQR": q3(x) - q1(x),
        "mapped_reads_minimum": x.min(),
        "mapped_reads_maximum": x.max(),
        "QC_PASS_n": g["source_QC_status"].astype(str).str.upper().eq("PASS").sum(),
        "QC_FAIL_n": g["source_QC_status"].astype(str).str.upper().eq("FAIL").sum(),
        "metric_definition": "Atlas source-reported transcriptome-mapped reads; not total reads, feature-count library size, or an RPM denominator",
    }


mapped_summary = []
for fluid, g in long.groupby("biofluid", sort=False):
    mapped_summary.append(summarize(g, "overall_by_biofluid", "ALL_D01-D07"))
    for label, gg in g.groupby("collection_label", sort=True):
        mapped_summary.append(summarize(gg, "biofluid_by_collection_label", label))

mapped_summary = pd.DataFrame(mapped_summary)
mapped_summary.to_csv(OUT / "Supplementary_Table_S9_source_reported_mapped_reads_and_QC.csv", index=False)

print("S8 feature rows:", len(feature))
print("S8 summary rows:", len(summary))
print("Mapped sample rows:", len(long))
print(mapped_summary.to_string(index=False))
