"""Turn the measurement-noise summary CSVs into LaTeX rows for the manuscript.

Usage: python scripts/noise_to_latex.py results/batches/noise_<id>

Prints the headline numbers for Section 3.5 and the body rows of
Tables tab:noise_params and tab:noise_post, ready to paste over the
\\pending{} placeholders in docs/manuscript/revision_20260927/main.tex.
Every replicate (including failed ones) is counted in the headline.
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
POST = {
    "SAP_max": "Systolic pressure, RFA (mmHg)", "SAP_min": "Diastolic pressure, RFA (mmHg)",
    "SAP_mean": "Mean arterial pressure, RFA (mmHg)", "PAP_max": "PA systolic pressure (mmHg)",
    "PAP_min": "PA diastolic pressure (mmHg)", "PAP_mean": "Mean PA pressure (mmHg)",
    "RAP_mean": "Mean RA pressure (mmHg)", "LVEDV": "LVEDV (mL)", "LVESV": "LVESV (mL)",
    "RVEDV": "RVEDV (mL)", "RVESV": "RVESV (mL)", "LVEF": "LVEF (--)", "RVEF": "RVEF (--)",
    "CO_Lmin": "Cardiac output (L/min)", "QpQs": "$Q_p/Q_s$ (--)", "SVR": "SVR (WU)",
}
# Reference post-closure prediction from the unperturbed calibration (Table tab:post).
POST_REF = {
    "SAP_max": 90.08, "SAP_min": 51.95, "SAP_mean": 70.19, "PAP_max": 23.92, "PAP_min": 9.29,
    "PAP_mean": 14.86, "RAP_mean": 5.41, "LVEDV": 41.10, "LVESV": 11.62, "RVEDV": 43.35,
    "RVESV": 13.78, "LVEF": 0.717, "RVEF": 0.682, "CO_Lmin": 3.17, "QpQs": 1.00, "SVR": 20.43,
}


def rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt(x, digits=3):
    return f"{float(x):.{digits}g}" if abs(float(x)) < 1 else f"{float(x):.{max(digits - 1, 1)}f}"


def main(batch):
    s = Path(batch) / "summary"
    runs = rows(s / "noise_runs.csv")
    done = [r for r in runs if r["Status"] == "completed"]
    print(f"% Replicates: {len(runs)} total, {len(done)} completed, {len(runs) - len(done)} failed")
    labels = {}
    for r in done:
        labels[r["CalibrationStatus"]] = labels.get(r["CalibrationStatus"], 0) + 1
    print("% Status counts:", labels)
    rmse = [float(r["GovernedRMSE_vsPerturbed"]) for r in done]
    gate = [float(r["GovernedWithinGate"]) for r in done]
    post10 = [float(r["PostWithin10"]) for r in done]
    if rmse:
        print(f"% Governed RMSE vs perturbed: median {statistics.median(rmse):.4f} "
              f"(range {min(rmse):.4f}-{max(rmse):.4f})")
        print(f"% Governed within gate: median {statistics.median(gate):g} of 9 "
              f"(range {min(gate):g}-{max(gate):g})")
        print(f"% Post within 10%: median {statistics.median(post10):g} of 16 "
              f"(range {min(post10):g}-{max(post10):g})")

    print("\n% ---- tab:noise_params body ----")
    for r in rows(s / "noise_summary_parameters.csv"):
        name, f = PARAM.get(r["Parameter"], (r["Parameter"], "{:.4g}"))
        print(f"{name} & {f.format(float(r['Reference']))} & "
              f"{f.format(float(r['Mean']))} $\\pm$ {f.format(float(r['SD']))} & "
              f"{float(r['CV_pct']):.1f} & [{f.format(float(r['P2_5']))}, {f.format(float(r['P97_5']))}] \\\\")

    print("\n% ---- tab:noise_post body ----")
    for r in rows(s / "noise_summary_post.csv"):
        m = r["Metric"]
        inside = "yes" if r["MeasuredInside95Band"].strip().lower() in ("1", "true") else "no"
        print(f"{POST.get(m, m)} & {fmt(r['Measured'])} & {fmt(POST_REF.get(m, float('nan')))} & "
              f"{fmt(r['Mean'])} [{fmt(r['P2_5'])}, {fmt(r['P97_5'])}] & "
              f"{100 * float(r['FractionWithin10']):.0f} & {inside} \\\\")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
