function [P, Q] = reconstruct_hemodynamic_signals(t, XV, params)
% RECONSTRUCT_HEMODYNAMIC_SIGNALS
% -----------------------------------------------------------------------
% Reconstructs pressure and flow signals from the 14-state model.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

sidx = params.idx;

V_RA   = XV(:, sidx.V_RA);
V_RV   = XV(:, sidx.V_RV);
V_LA   = XV(:, sidx.V_LA);
V_LV   = XV(:, sidx.V_LV);
V_SAR  = XV(:, sidx.V_SAR);
V_SVEN = XV(:, sidx.V_SVEN);
Q_SVEN = XV(:, sidx.Q_SVEN);
V_PAR  = XV(:, sidx.V_PAR);
V_PVEN = XV(:, sidx.V_PVEN);
Q_PVEN = XV(:, sidx.Q_PVEN);

[E_LV, E_RV, E_LA, E_RA] = elastance_model(t, params);

P.RA = max(E_RA .* (V_RA - params.V0.RA), -5);
P.RV = E_RV .* (V_RV - params.V0.RV);
P.LA = max(E_LA .* (V_LA - params.V0.LA), -5);
P.LV = E_LV .* (V_LV - params.V0.LV);
P.SAR = (V_SAR - params.V0.SAR) ./ params.C.SAR;
P.SVEN = max((V_SVEN - params.V0.SVEN) ./ params.C.SVEN, -5);
P.PAR = (V_PAR - params.V0.PAR) ./ params.C.PAR;
P.PVEN = max((V_PVEN - params.V0.PVEN) ./ params.C.PVEN, -5);

Q.TV = valve_model(P.RA, P.RV, params);
Q.PVv = valve_model(P.RV, P.PAR, params);
Q.MV = valve_model(P.LA, P.LV, params);
Q.AV = valve_model(P.LV, P.SAR, params);
Q.SVEN = Q_SVEN;
Q.PVEN = Q_PVEN;
Q.VSD = vsd_shunt_model(P.LV, P.RV, params);

end
