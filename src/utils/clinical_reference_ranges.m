function ref = clinical_reference_ranges(population)
% CLINICAL_REFERENCE_RANGES
% -----------------------------------------------------------------------
% Returns a struct of healthy haemodynamic reference intervals for use
% in simulation comparison dashboards.  Each field is a two-element
% vector [lo, hi] representing the inclusive normal range for that
% metric, or a scalar when only a single representative value is known.
%
% INPUTS:
%   population - 'adult' or 'pediatric'   (case-insensitive string)
%
% OUTPUTS:
%   ref        - struct with fields (units match the summary table):
%
%     ref.HR        [bpm]        Heart rate
%     ref.CO        [L/min]      Cardiac output
%     ref.LVEDV     [mL]         LV end-diastolic volume
%     ref.LVESV     [mL]         LV end-systolic volume
%     ref.RVEDV     [mL]         RV end-diastolic volume
%     ref.RVESV     [mL]         RV end-systolic volume
%     ref.SV        [mL]         Stroke volume (LV)
%     ref.LVEF      [%]          LV ejection fraction
%     ref.RVEF      [%]          RV ejection fraction
%     ref.SBP       [mmHg]       Systemic arterial systolic
%     ref.DBP       [mmHg]       Systemic arterial diastolic
%     ref.MAP       [mmHg]       Mean arterial pressure
%     ref.RAP_mean  [mmHg]       Mean right atrial pressure
%     ref.LAP_mean  [mmHg]       Mean left atrial pressure
%     ref.PAP_mean  [mmHg]       Mean pulmonary artery pressure
%     ref.PAP_max   [mmHg]       Pulmonary artery systolic
%     ref.PAP_min   [mmHg]       Pulmonary artery diastolic
%     ref.SVR       [WU]         Systemic vascular resistance
%     ref.PVR       [WU]         Pulmonary vascular resistance
%     ref.RVESP     [mmHg]       RV end-systolic pressure
%     ref.RVEDP     [mmHg]       RV end-diastolic pressure
%     ref.LVESP     [mmHg]       LV end-systolic pressure
%     ref.LVEDP     [mmHg]       LV end-diastolic pressure
%
% USAGE:
%   ref_ad  = clinical_reference_ranges('adult');
%   ref_ped = clinical_reference_ranges('pediatric');
%
% NOTES:
%   - Ranges are [lo, hi]; use ref.HR(1) / ref.HR(2) for bounds.
%   - NaN indicates no established normative range for that population.
%   - These are healthy-normal references; VSD patients are expected
%     to deviate — use for context, not as optimisation targets.
%
% REFERENCES:
%   Adult values from tabulated clinical reference spreadsheet
%   (Columna B — "Adult ref"):
%     HR        74.2 ± 12.7 → ~[61, 87]
%     LVEDV     112.1 ± 24.8 → ~[87, 137]
%     LVESV     41  ± 10    → ~[31, 51]
%     RVEDV     64.56–160.49
%     RVESV     18.45–81.17
%     SV        60–100
%     LVEF      52–72 %
%     RVEF      58 ± 8 %  → ~[50, 66]
%     SBP       142 ± 14  → ~[128, 156]
%     DBP       88  ± 8   → ~[80, 96]
%     MAP       70–100
%     RAP_mean  2–6
%     LAP_mean  2–12
%     PAP_mean  8–20
%     PAP_max   15–30
%     PAP_min   4–12
%     CO        4–8
%     SVR       9.06–31.46
%     PVR       0.3–2.0
%     RVESP     20–30
%     RVEDP     0–8
%     LVESP     90–140
%     LVEDP     3–12
%
%   Pediatric values from same spreadsheet (Columna E — "Pediatric ref"):
%     HR        73–142
%     LVEDV     32.9 ± 7.5 → ~[25, 41]
%     LVESV     13.0 ± 3.1 → ~[10, 16]
%     RVEDV     17.09 ± 2.70 → ~[14, 20]
%     RVESV     8.71 ± 1.78 → ~[7, 10]
%     SV        19.9 ± 4.9 → ~[15, 25]
%     LVEF      64 ± 9 %   → ~[55, 73]
%     RVEF      63 ± 5 %   → ~[58, 68]
%     SBP       90–103
%     DBP       47–59
%     MAP       61–74
%     RAP_mean  3–6
%     LAP_mean  2–10
%     PAP_mean  0–21      (<21)
%     PAP_max   0–32.9    (<32.9)
%     PAP_min   0–14.95   (<14.95)
%     CO        ~[2.20, 4.85]
%     SVR       13.76–32.99
%     PVR       0–6        (<6)
%     RVESP     0–35       (<35)
%     RVEDP     NaN        (not specified)
%     LVESP     NaN        (not specified)
%     LVEDP     NaN        (not specified)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-13
% VERSION:  1.0
% -----------------------------------------------------------------------

switch lower(strtrim(population))

    % -----------------------------------------------------------------
    case 'adult'
    % -----------------------------------------------------------------
        ref.HR       = [61,   87];      % bpm   — 74.2 ± 12.7
        ref.CO       = [4.0,  8.0];     % L/min
        ref.LVEDV    = [87,   137];     % mL    — 112.1 ± 24.8
        ref.LVESV    = [31,   51];      % mL    — 41 ± 10
        ref.RVEDV    = [64.6, 160.5];   % mL
        ref.RVESV    = [18.5, 81.2];    % mL
        ref.SV       = [60,   100];     % mL
        ref.LVEF     = [52,   72];      % %
        ref.RVEF     = [50,   66];      % %    — 58 ± 8
        ref.SBP      = [128,  156];     % mmHg — 142 ± 14
        ref.DBP      = [80,   96];      % mmHg — 88 ± 8
        ref.MAP      = [70,   100];     % mmHg
        ref.RAP_mean = [2,    6];       % mmHg
        ref.LAP_mean = [2,    12];      % mmHg
        ref.PAP_mean = [8,    20];      % mmHg
        ref.PAP_max  = [15,   30];      % mmHg
        ref.PAP_min  = [4,    12];      % mmHg
        ref.SVR      = [9.1,  31.5];    % WU
        ref.PVR      = [0.3,  2.0];     % WU
        ref.RVESP    = [20,   30];      % mmHg
        ref.RVEDP    = [0,    8];       % mmHg
        ref.LVESP    = [90,   140];     % mmHg
        ref.LVEDP    = [3,    12];      % mmHg

    % -----------------------------------------------------------------
    case 'pediatric'
    % -----------------------------------------------------------------
        ref.HR       = [73,   142];     % bpm
        ref.CO       = [2.20, 4.85];    % L/min
        ref.LVEDV    = [25,   41];      % mL    — 32.9 ± 7.5
        ref.LVESV    = [10,   16];      % mL    — 13.0 ± 3.1
        ref.RVEDV    = [14,   20];      % mL    — 17.09 ± 2.70
        ref.RVESV    = [7,    10];      % mL    — 8.71 ± 1.78
        ref.SV       = [15,   25];      % mL    — 19.9 ± 4.9
        ref.LVEF     = [55,   73];      % %     — 64 ± 9
        ref.RVEF     = [58,   68];      % %     — 63 ± 5
        ref.SBP      = [90,   103];     % mmHg
        ref.DBP      = [47,   59];      % mmHg
        ref.MAP      = [61,   74];      % mmHg
        ref.RAP_mean = [3,    6];       % mmHg
        ref.LAP_mean = [2,    10];      % mmHg
        ref.PAP_mean = [0,    21];      % mmHg  (<21)
        ref.PAP_max  = [0,    32.9];    % mmHg  (<32.9)
        ref.PAP_min  = [0,    14.95];   % mmHg  (<14.95)
        ref.SVR      = [13.8, 33.0];    % WU
        ref.PVR      = [0,    6.0];     % WU    (<6)
        ref.RVESP    = [0,    35];      % mmHg  (<35)
        ref.RVEDP    = [NaN,  NaN];     % mmHg  — not specified
        ref.LVESP    = [NaN,  NaN];     % mmHg  — not specified
        ref.LVEDP    = [NaN,  NaN];     % mmHg  — not specified

    otherwise
        error('clinical_reference_ranges: unknown population ''%s''. Use ''adult'' or ''pediatric''.', population);
end

end  % clinical_reference_ranges
