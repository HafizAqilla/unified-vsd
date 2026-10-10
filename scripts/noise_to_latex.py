"""Turn the measurement-noise summary CSVs into LaTeX rows for the manuscript.

Usage: python scripts/noise_to_latex.py results/batches/noise_<id>

Prints the headline numbers for Section 3.5 and the body rows of
Tables tab:noise_params and tab:noise_post, ready to paste over the
\\pending{} placeholders in docs/manuscript/revision_20260927/main.tex.
The manuscript reports values and errors only; no pass/fail counts are
printed. Quantity names match the pre- and post-closure tables.
"""
import csv
import statistics
import sys
from pathlib import Path

PARAM = {
    "group.R_sys_scale": (r"$R_{\mathrm{sys}}$ scale", "{:.4f}"),
    "R.SVEN": (r"$R_{\mathrm{SVEN}}$ (mmHg\,s/mL)", "{:.4f}"),
    "group.R_pul_scale": (r"$R_{\mathrm{pul}}$ scale", "{:.4f}"),
    "C.SAR": (r"$C_{\mathrm{SAR}}$ (mL/mmHg)", "{:.4f}"),
    "C.PAR": (r"$C_{\mathrm{PAR}}$ (mL/mmHg)", "{:.4f}"),
    "E.RV.EB": (r"$E_{\mathrm{B,RV}}$ (mmHg/mL)", "{:.4f}"),
    "vsd.Cd": (r"$C_d$ (--)", "{:.4f}"),
}
# Same names, symbols and order as Tables tab:pre_targets and tab:post.
POST_ORDER = [
    ("SAP_max", "Systolic pressure, RFA (mmHg)"),
    ("SAP_min", "Diastolic pressure, RFA (mmHg)"),
    ("SAP_mean", "Mean arterial pressure, RFA (mmHg)"),
    ("PAP_max", "PA systolic pressure (mmHg)"),
    ("PAP_min", "PA diastolic pressure (mmHg)"),
    ("PAP_mean", "Mean PA pressure (mmHg)"),
    ("RAP_mean", "Mean RA pressure (mmHg)"),
    ("CO_Lmin", "Systemic flow $Q_s$ (L/min)"),
    ("QpQs", "$Q_p/Q_s$ (--)"),
    ("SVR", "SVR (WU)"),
    ("LVEDV", "LVEDV (mL)"),
    ("LVESV", "LVESV (mL)"),
    ("RVEDV", "RVEDV (mL)"),
    ("RVESV", "RVESV (mL)"),
    ("LVEF", "LVEF (--)"),
    ("RVEF", "RVEF (--)"),
]
# Reference post-closure prediction from the unperturbed calibration (Table tab:post).
POST_REF = {
    "SAP_max": 90.08, "SAP_min": 51.95, "SAP_mean": 70.19, "PAP_max": 23.92, "PAP_min": 9.29,
    "PAP_mean": 14.86, "RAP_mean": 5.41, "LVEDV": 41.10, "LVESV": 11.62, "RVEDV": 43.35,
    "RVESV": 13.78, "LVEF": 0.717, "RVEF": 0.682, "CO_Lmin": 3.17, "QpQs": 1.00, "SVR": 20.43,
}
EOL = " \\\\"


def rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt(x, digits=3):
    x = float(x)
    return f"{x:.{digits}g}" if abs(x) < 1 else f"{x:.{max(digits - 1, 1)}f}"


def main(batch):
    s = Path(batch) / "summary"
    runs = rows(s / "noise_runs.csv")
    done = [r for r in runs if r["Status"] == "completed"]
    print(f"% Replicates: {len(runs)} total, {len(done)} recalibrated, "
          f"{len(runs) - len(done)} could not start calibration")
    rmse = [float(r["GovernedRMSE_vsPerturbed"]) for r in done]
    if rmse:
        print(f"% RMSE of relative error (nine fitted targets) vs perturbed record: median "
              f"{statistics.median(rmse):.4f} (range {min(rmse):.4f}-{max(rmse):.4f})")

    print("\n% ---- tab:noise_params body ----")
    for r in rows(s / "noise_summary_parameters.csv"):
        name, f = PARAM.get(r["Parameter"], (r["Parameter"], "{:.4g}"))
        print(f"{name} & {f.format(float(r['Reference']))} & "
              f"{f.format(float(r['Mean']))} $\\pm$ {f.format(float(r['SD']))} & "
              f"{float(r['CV_pct']):.1f} & [{f.format(float(r['P2_5']))}, {f.format(float(r['P97_5']))}]" + EOL)

    print("\n% ---- tab:noise_post body (Quantity | Measured | Reference prediction | Mean +- SD | 95% interval) ----")
    post = {r["Metric"]: r for r in rows(s / "noise_summary_post.csv")}
    for key, label in POST_ORDER:
        if key not in post:
            continue
        r = post[key]
        print(f"{label} & {fmt(r['Measured'])} & {fmt(POST_REF[key])} & "
              f"{fmt(r['Mean'])} $\\pm$ {fmt(r['SD'])} & [{fmt(r['P2_5'])}, {fmt(r['P97_5'])}]" + EOL)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
