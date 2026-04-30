function Q_VSD = vsd_shunt_model(P_LV, P_RV, params)
% VSD_SHUNT_MODEL
% -----------------------------------------------------------------------
% Configurable ventricular septal defect shunt model.
%
% MODES:
%   linear_left_to_right_only   legacy smooth-diode linear resistance
%   linear_bidirectional        symmetric linear resistance
%   orifice_bidirectional       reduced-order orifice law with signed flow
%
% SIGN CONVENTION:
%   Positive Q_VSD = LV -> RV (left-to-right shunt)
%   Negative Q_VSD = RV -> LV (right-to-left shunt)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  2.0
% -----------------------------------------------------------------------

dP = P_LV - P_RV;
mode = lower(params.vsd.mode);

switch mode
    case 'linear_left_to_right_only'
        gate = 0.5 + 0.5 * tanh(dP / params.epsilon_vsd);
        Q_VSD = gate .* dP / max(params.R.vsd, 1e-6);

    case 'linear_bidirectional'
        Q_VSD = dP / max(params.R.vsd, 1e-6);

    case 'orifice_bidirectional'
        if params.vsd.area_mm2 <= 0
            Q_VSD = zeros(size(dP));
            return;
        end
        A_m2 = params.vsd.area_mm2 * 1e-6;
        dP_Pa = abs(dP) * params.conv.mmHg_to_Pa;
        q_mag_m3s = params.vsd.Cd * A_m2 .* sqrt(2 * dP_Pa / params.vsd.rho_blood);
        q_mag_mLs = q_mag_m3s * params.conv.m3_to_mL;
        signed_gate = dP ./ sqrt(dP.^2 + params.epsilon_vsd^2);
        Q_VSD = q_mag_mLs .* signed_gate;

    otherwise
        error('vsd_shunt_model:unknownMode', ...
            'Unsupported VSD mode: %s', params.vsd.mode);
end

end
