function Y = gsa_pce_vsd_qoi(X, P)
% GSA_PCE_VSD_QOI
% -----------------------------------------------------------------------
% UQLab Quantity of Interest (QoI) wrapper for the Unified VSD model.
%
% This function is the equivalent of "QoI_CaseName.m" in the SoBioS
% toolbox pattern (Tosin, Côrtes, Cunha Jr, 2020). UQLab calls this
% function internally for each sample row during PCE training.
%
% INPUTS:
%   X   - [1 × d] row vector of parameter values for one sample
%         (UQLab passes one row at a time when isVectorized = false)
%   P   - struct of constant parameters (set in gsa_run_pce.m):
%           .p0     baseline parameter struct (params0)
%           .names  cell array of parameter names  {d × 1}
%           .mf     metric name to extract, e.g. 'QpQs'
%           .cfg    full GSA config struct (for reference)
%
% OUTPUT:
%   Y   - [1 × 1] scalar value of the requested metric for this sample
%
% REFERENCES:
%   SoBioS QoI_CaseName.m pattern — Tosin, Côrtes, Cunha Jr (2020)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-10
% VERSION:  1.0
% -----------------------------------------------------------------------

%% Apply sampled parameter values onto the baseline struct
params = P.p0;
for i = 1:numel(P.names)
    parts = strsplit(P.names{i}, '.');
    switch numel(parts)
        case 2
            params.(parts{1}).(parts{2})             = X(i);
        case 3
            params.(parts{1}).(parts{2}).(parts{3})  = X(i);
    end
end

%% Run cardiovascular model and extract requested metric
try
    sim     = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);

    if isfield(metrics, P.mf)
        Y = metrics.(P.mf);
    else
        Y = 0;   % metric not available for this param set
    end
catch
    Y = 0;       % failed simulation: return 0 (consistent with Sobol fallback)
end

end  % gsa_pce_vsd_qoi
