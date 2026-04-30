function metrics = compute_clinical_indices(sim, params)
% COMPUTE_CLINICAL_INDICES
% -----------------------------------------------------------------------
% Derive haemodynamic and geometry-proxy metrics from a steady-state run.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  2.0
% -----------------------------------------------------------------------

t  = sim.t(:);
XV = sim.V;
sidx = params.idx;
mLs_to_Lmin = params.conv.mLs_to_Lmin;

[P, Q] = reconstruct_hemodynamic_signals(t, XV, params);

T_HB = 60 / params.HR;
T1 = t(end);
T0 = T1 - T_HB;
time_mask = (t >= T0) & (t <= T1);

tc = t(time_mask);
Pc = struct();
Qc = struct();
flds_P = fieldnames(P);
flds_Q = fieldnames(Q);
for i = 1:numel(flds_P), Pc.(flds_P{i}) = P.(flds_P{i})(time_mask); end
for i = 1:numel(flds_Q), Qc.(flds_Q{i}) = Q.(flds_Q{i})(time_mask); end

mean_t = @(y) trapz(tc, y) / max(tc(end) - tc(1), 1e-9);

metrics = struct();

metrics.RAP_mean = mean_t(Pc.RA);
metrics.RAP_min  = min(Pc.RA);
metrics.RAP_max  = max(Pc.RA);
metrics.LAP_mean = mean_t(Pc.LA);
metrics.LAP_min  = min(Pc.LA);
metrics.LAP_max  = max(Pc.LA);

metrics.PAP_min  = min(Pc.PAR);
metrics.PAP_max  = max(Pc.PAR);
metrics.PAP_mean = mean_t(Pc.PAR);
metrics.PVP_mean = mean_t(Pc.PVEN);

metrics.RVP_min  = min(Pc.RV);
metrics.RVP_max  = max(Pc.RV);
metrics.RVP_mean = mean_t(Pc.RV);
metrics.LVP_min  = min(Pc.LV);
metrics.LVP_max  = max(Pc.LV);
metrics.LVP_mean = mean_t(Pc.LV);

metrics.SAP_min  = min(Pc.SAR);
metrics.SAP_max  = max(Pc.SAR);
metrics.SAP_mean = mean_t(Pc.SAR);

Qsys_mLs  = mean_t(Qc.SVEN);
Qpul_mLs  = mean_t(Qc.PVEN);
Qsys_Lmin = Qsys_mLs * mLs_to_Lmin;
Qpul_Lmin = Qpul_mLs * mLs_to_Lmin;
metrics.Qs_mean_mLs = Qsys_mLs;
metrics.Qp_mean_mLs = Qpul_mLs;
metrics.Qs_Lmin = Qsys_Lmin;
metrics.Qp_Lmin = Qpul_Lmin;

metrics.SVR  = (metrics.SAP_mean - metrics.RAP_mean) / max(Qsys_Lmin, 1e-6);
metrics.PVR  = (metrics.PAP_mean - metrics.LAP_mean) / max(Qpul_Lmin, 1e-6);
metrics.QpQs = Qpul_Lmin / max(Qsys_Lmin, 1e-6);

V_LV_c = XV(time_mask, sidx.V_LV);
V_RV_c = XV(time_mask, sidx.V_RV);
V_LA_c = XV(time_mask, sidx.V_LA);

metrics.LVEDV = max(V_LV_c);
metrics.LVESV = min(V_LV_c);
metrics.RVEDV = max(V_RV_c);
metrics.RVESV = min(V_RV_c);

metrics.LVEF = (metrics.LVEDV - metrics.LVESV) / max(metrics.LVEDV, 1e-6);
metrics.RVEF = (metrics.RVEDV - metrics.RVESV) / max(metrics.RVEDV, 1e-6);
metrics.LVSV = metrics.LVEDV - metrics.LVESV;
metrics.RVSV = metrics.RVEDV - metrics.RVESV;
metrics.LVCO_Lmin = metrics.LVSV * params.HR / 1000;
metrics.RVCO_Lmin = metrics.RVSV * params.HR / 1000;

metrics.Q_shunt_mean_mLs = mean_t(Qc.VSD);
metrics.VSD_frac_pct = 100 * metrics.Q_shunt_mean_mLs / max(abs(metrics.Qp_mean_mLs), 1e-6);
% CO_Lmin is reported as effective systemic output (Qs). In unrepaired VSD,
% LV stroke output includes recirculated shunt volume and overstates the
% clinically reported systemic cardiac output.
metrics.CO_Lmin = Qsys_Lmin;

metrics.LVEDD_mm = teichholz_diameter_mm(metrics.LVEDV);
metrics.LVESD_mm = teichholz_diameter_mm(metrics.LVESV);
metrics.RV_diameter_proxy_mm = sphere_diameter_mm(metrics.RVEDV);
metrics.LA_size_proxy_mm = sphere_diameter_mm(max(V_LA_c));
metrics.LV_sphericity_index = metrics.LVESD_mm / max(metrics.LVEDD_mm, 1e-6);

end

function diameter_mm = teichholz_diameter_mm(volume_mL)
diameter_cm = max((6 * max(volume_mL, 0) / pi)^(1/3), 1e-6);
for iter = 1:6
    f = 7 * diameter_cm^3 / (2.4 + diameter_cm) - volume_mL;
    df = 7 * diameter_cm^2 * (7.2 + 2 * diameter_cm) / (2.4 + diameter_cm)^2;
    diameter_cm = max(diameter_cm - f / max(df, 1e-6), 1e-6);
end
diameter_mm = 10 * diameter_cm;
end

function diameter_mm = sphere_diameter_mm(volume_mL)
diameter_mm = 10 * (6 * max(volume_mL, 0) / pi)^(1/3);
end
