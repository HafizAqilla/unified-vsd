function comparison = compare_scaling_methods(params_ref, patient)
% COMPARE_SCALING_METHODS
% -----------------------------------------------------------------------
% Produces a compact Zhang-vs-Lundquist scaling comparison for one patient.
%
% INPUTS:
%   params_ref  - adult/reference parameter struct                      [-]
%   patient     - patient demographics struct                           [-]
%
% OUTPUTS:
%   comparison  - table of key scaling outcomes and scale factors        [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-14
% VERSION:  1.0
% -----------------------------------------------------------------------

params_zhang = apply_physiological_scaling(params_ref, patient, 'zhang');
params_lundquist = apply_physiological_scaling(params_ref, patient, 'lundquist');

rows = {
    'HR_bpm', params_zhang.HR, params_lundquist.HR, 'bpm'
    'BSA_patient_m2', params_zhang.scaling.BSA_patient, params_lundquist.scaling.BSA_patient, 'm2'
    'R_SAR_factor', params_zhang.R.SAR / params_ref.R.SAR, params_lundquist.R.SAR / params_ref.R.SAR, '-'
    'R_PAR_factor', params_zhang.R.PAR / params_ref.R.PAR, params_lundquist.R.PAR / params_ref.R.PAR, '-'
    'C_SAR_factor', params_zhang.C.SAR / params_ref.C.SAR, params_lundquist.C.SAR / params_ref.C.SAR, '-'
    'C_PAR_factor', params_zhang.C.PAR / params_ref.C.PAR, params_lundquist.C.PAR / params_ref.C.PAR, '-'
    'E_LV_EA_factor', params_zhang.E.LV.EA / params_ref.E.LV.EA, params_lundquist.E.LV.EA / params_ref.E.LV.EA, '-'
    'E_RV_EA_factor', params_zhang.E.RV.EA / params_ref.E.RV.EA, params_lundquist.E.RV.EA / params_ref.E.RV.EA, '-'
    'V0_LV_factor', params_zhang.V0.LV / params_ref.V0.LV, params_lundquist.V0.LV / params_ref.V0.LV, '-'
    'V0_RV_factor', params_zhang.V0.RV / params_ref.V0.RV, params_lundquist.V0.RV / params_ref.V0.RV, '-'
    'L_SAR_factor', NaN, params_lundquist.L.SAR / params_ref.L.SAR, '-'
    'Rvalve_open_factor', params_zhang.Rvalve.open / params_ref.Rvalve.open, ...
        params_lundquist.Rvalve.open / params_ref.Rvalve.open, '-'
    };

comparison = cell2table(rows, ...
    'VariableNames', {'Quantity','Zhang','Lundquist','Unit'});
end
