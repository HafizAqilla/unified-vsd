function HR_bpm = clinical_heart_rate(clinical, scenario)
% CLINICAL_HEART_RATE
% -----------------------------------------------------------------------
% Heart rate recorded for one clinical scenario.
%
% The scenario's own measurement (clinical.<scenario>.HR) takes precedence,
% so a pre-closure and a post-closure state recorded at different rates each
% run at their own rate. clinical.common.HR is used only when the scenario
% carries no rate of its own (patients whose record has a single HR).
%
% INPUTS:
%   clinical - unified clinical struct                                [-]
%   scenario - 'pre_surgery' | 'post_surgery'                         [-]
%
% OUTPUTS:
%   HR_bpm   - heart rate, NaN when neither source is available    [bpm]
% -----------------------------------------------------------------------

HR_bpm = NaN;
if isfield(clinical, scenario) && isfield(clinical.(scenario), 'HR') && ...
        isfinite(clinical.(scenario).HR)
    HR_bpm = clinical.(scenario).HR;
elseif isfield(clinical, 'common') && isfield(clinical.common, 'HR') && ...
        isfinite(clinical.common.HR)
    HR_bpm = clinical.common.HR;
end
end
