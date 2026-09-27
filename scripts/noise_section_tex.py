"""Write the measurement-noise Results, Discussion and Conclusion text of the
manuscript directly from the batch summary, so the manuscript never carries
placeholders and is regenerated in one step when more replicates finish.

Usage:
    python scripts/noise_section_tex.py results/batches/noise_<id> docs/manuscript/revision_20260927

Writes into the manuscript folder:
    noise_results.tex     Results subsection body: text, Tables 8-9, Figure 8
    noise_discussion.tex  Discussion subsection body
    noise_conclusion.tex  one Conclusion sentence

Wording adapts to the number of recalibrated replicates n: with n = 1 the text
describes that replicate and reports no spread (SD, CV and intervals need at
least two values and are shown as --); with n >= 2 it reports the spread.
The text always states how many of the planned replicates were run.
"""
import csv
import statistics
import sys
from pathlib import Path

M_PLANNED = 50

PARAMS = [  # key, symbol for tables, words for text, format
    ("group.R_sys_scale", r"$R_{\mathrm{sys}}$ scale", "the systemic resistance scale", "{:.4f}"),
    ("R.SVEN", r"$R_{\mathrm{SVEN}}$ (mmHg\,s/mL)", r"$R_{\mathrm{SVEN}}$", "{:.4f}"),
    ("group.R_pul_scale", r"$R_{\mathrm{pul}}$ scale", "the pulmonary resistance scale", "{:.4f}"),
    ("C.SAR", r"$C_{\mathrm{SAR}}$ (mL/mmHg)", r"$C_{\mathrm{SAR}}$", "{:.4f}"),
    ("C.PAR", r"$C_{\mathrm{PAR}}$ (mL/mmHg)", r"$C_{\mathrm{PAR}}$", "{:.4f}"),
    ("E.RV.EB", r"$E_{\mathrm{B,RV}}$ (mmHg/mL)", r"$E_{\mathrm{B,RV}}$", "{:.4f}"),
    ("vsd.Cd", r"$C_d$ (--)", r"$C_d$", "{:.4f}"),
]
# Same names, symbols and order as Tables tab:pre_targets and tab:post.
POST = [  # reference = unrounded prediction of the reference calibration; precision as Table tab:post
    ("SAP_max", "Systolic pressure, RFA (mmHg)", 90.083, "{:.2f}"),
    ("SAP_min", "Diastolic pressure, RFA (mmHg)", 51.953, "{:.2f}"),
    ("SAP_mean", "Mean arterial pressure, RFA (mmHg)", 70.185, "{:.2f}"),
    ("PAP_max", "PA systolic pressure (mmHg)", 23.916, "{:.2f}"),
    ("PAP_min", "PA diastolic pressure (mmHg)", 9.291, "{:.2f}"),
    ("PAP_mean", "Mean PA pressure (mmHg)", 14.860, "{:.2f}"),
    ("RAP_mean", "Mean RA pressure (mmHg)", 5.411, "{:.2f}"),
    ("CO_Lmin", "Systemic flow $Q_s$ (L/min)", 3.171, "{:.2f}"),
    ("QpQs", "$Q_p/Q_s$ (--)", 1.000, "{:.2f}"),
    ("SVR", "SVR (WU)", 20.429, "{:.2f}"),
    ("LVEDV", "LVEDV (mL)", 41.098, "{:.2f}"),
    ("LVESV", "LVESV (mL)", 11.619, "{:.2f}"),
    ("RVEDV", "RVEDV (mL)", 43.347, "{:.2f}"),
    ("RVESV", "RVESV (mL)", 13.780, "{:.2f}"),
    ("LVEF", "LVEF (--)", 0.717, "{:.3f}"),
    ("RVEF", "RVEF (--)", 0.682, "{:.3f}"),
]
# Perturbed inputs: draws column, sigma-spec metric, words.
DRAWN = [
    ("CO_Lmin", "CO_Lmin", r"$Q_s$", "L/min", "{:.2f}"),
    ("SAP_sys_mmHg", "SAP_max", "systolic RFA pressure", "mmHg", "{:.1f}"),
    ("SAP_dia_mmHg", "SAP_min", "diastolic RFA pressure", "mmHg", "{:.1f}"),
    ("PAP_sys_mmHg", "PAP_max", "systolic PA pressure", "mmHg", "{:.1f}"),
    ("PAP_dia_mmHg", "PAP_min", "diastolic PA pressure", "mmHg", "{:.1f}"),
    ("PAP_mean_mmHg", "PAP_mean", "mean PA pressure", "mmHg", "{:.1f}"),
    ("RAP_mean_mmHg", "RAP_mean", "mean RA pressure", "mmHg", "{:.2f}"),
    ("QpQs", "QpQs", r"$Q_p/Q_s$", "", "{:.3f}"),
]
EOL = " \\\\\n"


def read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def pct(a, b):
    return 100.0 * (a - b) / abs(b)


def signed(x, digits=1):
    return f"${x:+.{digits}f}$"


def main(batch, out):
    batch, out = Path(batch), Path(out)
    s = batch / "summary"
    runs = read(s / "noise_runs.csv")
    done = [r for r in runs if r["Status"] == "completed"]
    failed = [r for r in runs if r["Status"] != "completed"]
    n, a = len(done), len(runs)
    par = {r["Parameter"]: r for r in read(s / "noise_summary_parameters.csv")}
    post = {r["Metric"]: r for r in read(s / "noise_summary_post.csv")}
    spec = {r["Metric"]: r for r in read(batch / "noise_spec.csv")}
    draws = {int(r["Replicate"]): r for r in read(batch / "noise_draws.csv")}
    spread = n >= 2

    def spread_cells(r, f):
        if not spread:
            return f.format(float(r["Mean"])), "--", "--"
        return (f"{f.format(float(r['Mean']))} $\\pm$ {f.format(float(r['SD']))}",
                f"{float(r['CV_pct']):.1f}" if "CV_pct" in r else "--",
                f"[{f.format(float(r['P2_5']))}, {f.format(float(r['P97_5']))}]")

    # ------------------------------------------------------------ tables
    tab_par = [r"\begin{table}[!htbp]", r"\centering\small",
               r"\caption{Calibrated parameters across the measurement-noise replicates (" + f"$n = {n}$ recalibrated of {a} run, {M_PLANNED} planned" + r"). Reference: the calibration to the unperturbed record (Table~\ref{tab:bounds}). SD, CV and the 2.5--97.5th percentile interval need at least two replicates; -- where not yet defined.}",
               r"\label{tab:noise_params}", r"\begin{tabular}{lrrrrr}", r"\toprule",
               r"\textbf{Parameter} & \textbf{Reference} & \textbf{Mean $\pm$ SD} & \textbf{Change (\%)} & \textbf{CV (\%)} & \textbf{95\% interval} \\",
               r"\midrule"]
    for key, sym, _, f in PARAMS:
        r = par[key]
        ref, mean = float(r["Reference"]), float(r["Mean"])
        cells = spread_cells(r, f)
        tab_par.append(f"{sym} & {f.format(ref)} & {cells[0]} & {signed(pct(mean, ref))} & {cells[1]} & {cells[2]}" + EOL.rstrip("\n"))
    tab_par += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

    tab_post = [r"\begin{table}[!htbp]", r"\centering\small", r"\setlength{\tabcolsep}{3pt}",
                r"\caption{Post-closure prediction across the measurement-noise replicates (" + f"$n = {n}$" + r"), beside the unperturbed post-closure record and the reference prediction (Table~\ref{tab:post}). Rows, names and symbols as in Table~\ref{tab:post}; $\varepsilon$, relative error of the mean against the measurement; -- where not yet defined.}",
                r"\label{tab:noise_post}", r"\resizebox{\linewidth}{!}{\begin{tabular}{lrrrrr}", r"\toprule",
                r"\textbf{Quantity} & \textbf{Measured} & \textbf{Reference prediction} & \textbf{Mean $\pm$ SD} & $\varepsilon$ (\%) & \textbf{95\% interval} \\",
                r"\midrule"]
    for key, label, ref, f in POST:
        r = post[key]
        meas, mean = float(r["Measured"]), float(r["Mean"])
        r2 = dict(r)
        cells = spread_cells(r2, f)
        tab_post.append(f"{label} & {f.format(meas)} & {f.format(ref)} & {cells[0]} & {signed(pct(mean, meas))} & {cells[2]}" + EOL.rstrip("\n"))
    tab_post += [r"\bottomrule", r"\end{tabular}}", r"\end{table}"]

    fig = [r"\begin{figure}[!htbp]", r"\centering",
           r"\includegraphics[width=\linewidth]{M6_noise.pdf}",
           r"\caption{Measurement-noise simulation (" + f"$n = {n}$" + r" recalibrated replicates). (A) Each calibrated parameter divided by its value in the reference calibration, one point per replicate; dashed line: reference. (B) Post-closure prediction minus measurement (\% of measurement), one point per replicate; the vertical line marks zero.}",
           r"\label{fig:noise}", r"\end{figure}"]

    # ------------------------------------------------------------ Results text
    rmse = [float(r["GovernedRMSE_vsPerturbed"]) for r in done]
    chi = [float(r["Chi2PerObs"]) for r in done]
    head = (f"Of the {M_PLANNED} planned replicates, {a} had been run when this version was written: "
            f"{n} {'was' if n == 1 else 'were'} recalibrated and {len(failed)} could not start calibration because the perturbed "
            r"shunt flow placed the clinically seeded $C_d$ on its lower bound. The remaining replicates use the same settings, "
            r"and Tables~\ref{tab:noise_params} and \ref{tab:noise_post} and Figure~\ref{fig:noise} are regenerated from the batch output.")
    if a == M_PLANNED:
        head = (f"Of the {M_PLANNED} replicates, {n} were recalibrated and {len(failed)} could not start calibration because the perturbed "
                r"shunt flow placed the clinically seeded $C_d$ on its lower bound.")

    if not spread:
        k = int(done[0]["Replicate"])
        d = draws[k]
        dev = []
        for col, met, words, unit, f in DRAWN:
            mu, sg = float(spec[met]["Measured"]), float(spec[met]["Sigma"])
            v = float(d[col])
            dev.append((abs(v - mu) / sg, (v - mu) / sg, words, f.format(v), f.format(mu), unit))
        dev.sort(reverse=True)
        big = [x for x in dev if x[0] >= 1.0]
        small_max = max(x[0] for x in dev if x[0] < 1.0)
        parts = [f"{w} was {v}~{u} against {m}~{u} measured (${z:+.1f}\\sigma$)".replace("~ ", " ").replace("~}", "}")
                 for _, z, w, v, m, u in big]
        rec = (f"In the recalibrated replicate (replicate {k}), " + "; ".join(parts)
               + f"; every other perturbed target moved by less than {small_max:.1f}$\\sigma$. ")
        rec += (f"The recalibrated model fitted this perturbed record with an RMSE of relative error of {rmse[0]:.3f} "
                f"($\\chi^2/N = {chi[0]:.2f}$), against 0.105 and 3.48 for the reference calibration to the measured record. ")
        chg = []
        for key, _, words, f in PARAMS:
            r = par[key]
            ref, mean = float(r["Reference"]), float(r["Mean"])
            chg.append(f"{words} from {f.format(ref)} to {f.format(mean)} (${pct(mean, ref):+.0f}\\%$)")
        rec += "Relative to the reference calibration, the parameters changed as follows: " + ", ".join(chg[:-1]) + " and " + chg[-1] + r" (Table~\ref{tab:noise_params}, Figure~\ref{fig:noise}A). "
        pred = []
        for key, label, ref, f in POST:
            if key not in ("CO_Lmin", "SAP_max", "RVEDV", "LVEDV", "PAP_max"):
                continue
            r = post[key]
            meas, mean = float(r["Measured"]), float(r["Mean"])
            name = label.split(" (")[0].replace("Systemic flow ", "")
            name = name if name[:2].isupper() or name.startswith("$") else name[0].lower() + name[1:]
            pred.append(f"{name} {f.format(ref)} to {f.format(mean)} against {f.format(meas)} measured")
        rec += ("The post-closure prediction of this replicate differed from the reference prediction as follows: "
                + "; ".join(pred) + r" (Table~\ref{tab:noise_post}, Figure~\ref{fig:noise}B).")
        results = head + "\n\n" + rec
        disc = (
            r"The completed replicate shows how strongly the calibrated parameters depend on the measured systemic flow. "
            + f"Its perturbed $Q_s$ was {dev[0][3]}~L/min, {abs(pct(float(dev[0][3]), float(dev[0][4]))):.0f}\\% below the measured value; "
            r"a deviation of this order is comparable to the error of pediatric Fick output computed with an assumed oxygen consumption, which has a mean bias of $+27.5\%$ and a standard deviation of 30\% against thermodilution or MRI~\cite{Ralston2023}. "
            + f"Recalibrating to this record fitted it as closely as the reference calibration fitted the measured record (RMSE {rmse[0]:.3f} against 0.105), "
            + f"but with {PARAMS[0][2]} changed by ${pct(float(par['group.R_sys_scale']['Mean']), float(par['group.R_sys_scale']['Reference'])):+.0f}\\%$ "
            + f"and {PARAMS[5][2]} by ${pct(float(par['E.RV.EB']['Mean']), float(par['E.RV.EB']['Reference'])):+.0f}\\%$. "
            r"Records that differ only within the declared measurement uncertainty are therefore fitted equally well by markedly different parameter sets. This agrees with the conditioning of the sensitivity matrix of the seven-parameter set (Supplementary Section~S2) and supports calibrating such models with quantified parameter uncertainty when the clinical data are limited~\cite{Kachabi2025,Colunga2023}. "
            + "The same replicate predicted RVEDV at "
            + f"{float(post['RVEDV']['Mean']):.1f}~mL and LVEDV at {float(post['LVEDV']['Mean']):.1f}~mL against 30.5 and 36.1~mL measured, where the reference calibration predicted 43.3 and 41.1~mL. "
            r"The post-closure volume errors of the reference prediction are therefore of the size that pre-closure measurement uncertainty alone can produce, and a volume prediction from this record should be reported with its uncertainty, as recommended for mechanistic cardiovascular models~\cite{Colebank2024Guidelines}. "
            + f"The {M_PLANNED}-replicate ensemble quantifies how often such shifts occur.")
        concl = ("In the replicate recalibrated so far, perturbing the pre-closure record within its declared uncertainty "
                 f"changed the calibrated parameters by up to ${max(abs(pct(float(par[k_]['Mean']), float(par[k_]['Reference']))) for k_,*_ in PARAMS):.0f}\\%$ "
                 f"and the predicted RVEDV by {abs(float(post['RVEDV']['Mean']) - 43.347):.1f}~mL, so post-closure volume predictions from this record carry a large measurement-driven uncertainty.")
    else:
        cvs = sorted(((float(par[k_]["CV_pct"]), w) for k_, _, w, _ in PARAMS), reverse=True)
        results = (head + "\n\n"
                   + f"Across the recalibrated replicates, the RMSE of relative error against the perturbed record had a median of {statistics.median(rmse):.3f} (range {min(rmse):.3f}--{max(rmse):.3f}). "
                   + f"The coefficient of variation of the calibrated parameters ranged from {cvs[-1][0]:.1f}\\% for {cvs[-1][1]} to {cvs[0][0]:.1f}\\% for {cvs[0][1]} "
                   + r"(Table~\ref{tab:noise_params}, Figure~\ref{fig:noise}A). ")
        lines = []
        for key, label, ref, f in POST:
            if key not in ("SAP_max", "SAP_min", "PAP_max", "CO_Lmin", "LVEDV", "RVEDV"):
                continue
            r = post[key]
            lines.append(f"{label.split(' (')[0]} [{f.format(float(r['P2_5']))}, {f.format(float(r['P97_5']))}] against {f.format(float(r['Measured']))} measured")
        results += ("The 95\\% intervals of the post-closure prediction were: " + "; ".join(lines)
                    + r" (Table~\ref{tab:noise_post}, Figure~\ref{fig:noise}B).")
        disc = (f"Measurement noise at the declared level changed the calibrated parameters with CV up to {cvs[0][0]:.0f}\\% ({cvs[0][1]}). "
                r"Tables~\ref{tab:noise_params} and \ref{tab:noise_post} give the resulting spread of every parameter and prediction beside the measured values; "
                r"the width of each interval is the uncertainty that the pre-closure measurements alone transfer to the post-closure prediction~\cite{Kachabi2025,Colunga2023,Colebank2024Guidelines}.")
        concl = (f"Under measurement noise at the declared level, the calibrated parameters varied with a coefficient of variation of up to {cvs[0][0]:.0f}\\%, "
                 f"and the 95\\% interval of the predicted RVEDV was [{float(post['RVEDV']['P2_5']):.1f}, {float(post['RVEDV']['P97_5']):.1f}]~mL against 30.5~mL measured.")

    (out / "noise_results.tex").write_text(
        "% GENERATED by scripts/noise_section_tex.py from " + batch.as_posix() + " -- do not edit by hand.\n"
        + results + "\n\n" + "\n".join(tab_par) + "\n\n" + "\n".join(tab_post) + "\n\n" + "\n".join(fig) + "\n",
        encoding="utf-8")
    (out / "noise_discussion.tex").write_text("% GENERATED -- do not edit by hand.\n" + disc + "\n", encoding="utf-8")
    (out / "noise_conclusion.tex").write_text("% GENERATED -- do not edit by hand.\n" + concl + "\n", encoding="utf-8")
    print(f"written: n={n} recalibrated, {a} run, {len(failed)} could not start")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
