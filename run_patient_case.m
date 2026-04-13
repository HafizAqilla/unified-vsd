% RUN_PATIENT_CASE
% Simulation and calibration setup for specific 3-year-old VSD patient.

clear; clc;

% 1. Load the empty template
clinical = patient_template();

% 2. Demographics
clinical.common.age_years = 3 + (2 / 12); % 3 years 2 months
clinical.common.weight_kg = 13.4;          % kg 
clinical.common.height_cm = 95;           % cm
clinical.common.sex       = 'F';          % Female

% 3. VSD defect data (Pre-surgery)
clinical.pre_surgery.VSD_diameter_mm   = 6.0;   % mm
clinical.pre_surgery.VSD_gradient_mmHg = 94;    % mmHg (Peak systolic VSD gradient)

% 4. Haemodynamic targets (Pre-surgery)
clinical.pre_surgery.PAP_sys_mmHg  = 20;      % mmHg (PAP Max)
clinical.pre_surgery.PAP_dia_mmHg  = 10;      % mmHg (PAP Min)
clinical.pre_surgery.PAP_mean_mmHg = 15;      % mmHg (PAP Mean)
clinical.pre_surgery.SAP_sys_mmHg  = 119;     % mmHg (SAP Max)
clinical.pre_surgery.SAP_dia_mmHg  = 83;      % mmHg (SAP Min)
clinical.pre_surgery.SAP_mean_mmHg = 95;      % mmHg (MAP)
clinical.pre_surgery.RAP_mean_mmHg = 5;       % mmHg
clinical.pre_surgery.QpQs          = 1.194;   % Pulmonary/Systemic Flow Ratio

% 5. Run the full pipeline
scenario = 'pre_surgery';
fprintf('Starting calibration for patient...\n');
main_run(scenario, clinical);
