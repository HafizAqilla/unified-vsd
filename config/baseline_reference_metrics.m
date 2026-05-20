function reference = baseline_reference_metrics()
% BASELINE_REFERENCE_METRICS
% -----------------------------------------------------------------------
% Returns the adult and pediatric scaling baseline metric snapshot used to
% reconcile the May 2026 Keisya/Hafiz baseline discussion.
%
% INPUTS:
%   none
%
% OUTPUTS:
%   reference - table of adult, Zhang, and Lundquist baseline metrics    [-]
%
% ASSUMPTIONS:
%   - Adult_ref is generated from default_parameters() without patient
%     scaling.
%   - Zhang_2019 and Lundquist_2025 are pediatric scaled-baseline outputs
%     for the spreadsheet reference patient in SHEET.xlsx.
%
% REFERENCES:
%   [1] SHEET.xlsx, sheet demo_metrics_20260518_144652.
%   [2] Screenshot shared in Codex thread, 2026-05-20.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-20
% VERSION:  1.0
% -----------------------------------------------------------------------

metric = { ...
    'Heart Rate'
    'Cardiac Output'
    'LV Stroke Volume'
    'RV Stroke Volume'
    'LVEF'
    'RVEF'
    'LVEDV'
    'LVESV'
    'RVEDV'
    'RVESV'
    'SVR'
    'PVR'
    'SBP'
    'DBP'
    'MAP'
    'PAP_sys'
    'PAP_dia'
    'PAP_mean'
    'LVESP'
    'LVEDP'
    'RVESP'
    'RVEDP'
    'LAP_mean'
    'RAP_mean'
    'PWP_mean'
    'Qp_Qs'};

unit = { ...
    'bpm'
    'L/min'
    'mL'
    'mL'
    '%'
    '%'
    'mL'
    'mL'
    'mL'
    'mL'
    'WU'
    'WU'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    '-'};

norm_adult = { ...
    '61-87'
    '4-8'
    '60-100'
    '--'
    '52-72'
    '50-66'
    '87-137'
    '31-51'
    '64.6-160.5'
    '18.5-81.2'
    '9.1-31.5'
    '0.3-2'
    '128-156'
    '80-96'
    '70-100'
    '15-30'
    '4-12'
    '8-20'
    '90-140'
    '3-12'
    '20-30'
    '<8'
    '2-12'
    '2-6'
    '--'
    '1.0 (normal)'};

adult_ref = [ ...
    75
    5.38797099834
    71.8396133112
    71.7762928572
    68.6593924986
    60.3230976656
    104.631880209
    32.7922668975
    118.986417533
    47.2101246758
    15.0481422596
    1.49984982493
    107.06725749
    65.9809285498
    86.0259678766
    21.6901637587
    12.1035780847
    16.5125607223
    107.66011476
    8.33844753264
    22.5049728405
    5.10919609536
    8.43420124093
    5.01748872314
    10.2295716992
    1.00052468422];

norm_pediatric = { ...
    '73-142'
    '2.2-4.85'
    '15-25'
    '--'
    '55-73'
    '45-78'
    '25-41'
    '10-16'
    '30.73-65.68'
    '4.2-26.18'
    '13.8-33'
    '<6'
    '90-103'
    '47-59'
    '61-74'
    '<32.9'
    '<14.95'
    '<21'
    '--'
    '--'
    '<35'
    '--'
    '2-10'
    '3-6'
    '--'
    '1.0 (normal)'};

zhang_2019 = [ ...
    121.549244752
    2.15426137243
    17.7233628792
    17.7423802146
    57.4382273069
    55.6439781394
    30.8563890464
    13.1330261672
    31.885535161
    14.1431549464
    32.2680643091
    4.62967905965
    99.7069281716
    51.0064161179
    73.9147946158
    22.1255569427
    10.4712235556
    15.8004288878
    100.282154174
    5.58135280528
    22.8009516651
    4.52924056124
    5.81823721765
    4.41723747805
    8.03557973988
    1.00110214842];

lundquist_2025 = [ ...
    105.376303938
    1.29231633427
    12.2638229466
    12.3424803453
    61.9840954403
    65.234900913
    19.7854350532
    7.52161210654
    18.9200568599
    6.57757651467
    42.0226026671
    4.21380482125
    68.6357539069
    48.1918351644
    58.2849727943
    12.5325942858
    7.70221315398
    10.0206951923
    69.0205603264
    4.47143713463
    13.0865754636
    3.71755336817
    4.56484073363
    3.64891851534
    5.77429141512
    0.995845533599];

reference = table(metric, unit, norm_adult, adult_ref, norm_pediatric, ...
    zhang_2019, lundquist_2025, ...
    'VariableNames', {'Metric','Unit','Norm_Adult','Adult_ref', ...
    'Norm_Pediatric','Zhang_2019','Lundquist_2025'});

reference.Properties.Description = ...
    'Baseline metrics snapshot from SHEET.xlsx demo_metrics_20260518_144652.';
reference.Properties.UserData.SourceFile = 'SHEET.xlsx';
reference.Properties.UserData.SourceSheet = 'demo_metrics_20260518_144652';
reference.Properties.UserData.SourceDate = '2026-05-20';
end
