function uc = unit_conversion()
% UNIT_CONVERSION
% -----------------------------------------------------------------------
% Central repository for all unit conversion factors used in the
% cardiovascular simulation codebase.
%
% The project internal unit system (Guardrail §6.1):
%   Pressure      [mmHg]       (internal and clinical reporting)
%   Volume        [mL]         (internal and clinical reporting)
%   Flow (ODE)    [mL/s]       (internal);  [L/min] for clinical reporting
%   Time          [s]
%   Resistance    [mmHg·s/mL]  (internal);  Wood units for clinical reporting
%   Compliance    [mL/mmHg]
%   Elastance     [mmHg/mL]
%   Diameter      [mm]
%   Area          [mm²]
%   BSA           [m²]
%
% Unit conversions MUST only be applied via this function or via the
% pre-computed conversion constants in params.conv (default_parameters.m).
% NEVER scatter raw conversion numbers (e.g., 60/1000) in physics code.
%
% USAGE:
%   uc = unit_conversion();
%   Q_Lmin = Q_mLs * uc.mLs_to_Lmin;
%   R_WU   = R_mech * uc.R_to_WU;
%
% OUTPUTS:
%   uc  - struct of named, commented conversion factors
%
% REFERENCES:
%   [1] SI unit definitions.
%   [2] Wood (1947). Wood units for vascular resistance.
%   [3] Guardrail §6.1–6.5 — units policy.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% VERSION:  1.0
% -----------------------------------------------------------------------

uc = struct();

%% Flow conversions
uc.mLs_to_Lmin  = 60   / 1000;   % [L/min per mL/s]  = 0.06
uc.Lmin_to_mLs  = 1000 / 60;     % [mL/s per L/min]  ≈ 16.667

%% Resistance / Wood unit conversions
% 1 Wood unit [WU] = 1 mmHg / (L/min) = mmHg / (1000/60 mL/s) = 0.06 mmHg·s/mL
uc.WU_to_R      = 60   / 1000;   % [mmHg·s/mL  per WU]  = 0.06
uc.R_to_WU      = 1000 / 60;     % [WU  per mmHg·s/mL]  ≈ 16.667

%% Pressure conversions
uc.mmHg_to_Pa   = 133.322;       % [Pa per mmHg]  — exact value, 1 mmHg = 133.322 Pa
uc.Pa_to_mmHg   = 1 / 133.322;  % [mmHg per Pa]

%% Volume conversions
uc.mL_to_m3     = 1e-6;          % [m³ per mL]
uc.m3_to_mL     = 1e6;           % [mL per m³]

%% Area conversions
uc.mm2_to_m2    = 1e-6;          % [m² per mm²]
uc.m2_to_mm2    = 1e6;           % [mm² per m²]

%% Length conversions
uc.mm_to_m      = 1e-3;          % [m per mm]
uc.m_to_mm      = 1e3;           % [mm per m]

%% Ejection fraction — stored as fraction; multiply by 100 only for display
% EF_pct = EF_frac * uc.EF_frac_to_pct   (display only)
uc.EF_frac_to_pct   = 100;        % [% per fraction]
uc.EF_pct_to_frac   = 1 / 100;   % [fraction per %]

end
