function [Ptout_dBm, Prout_dBm] = powercontrol(link_type, Pt_dBm, Pr_dBm, BS_sensitivity_dBm, MS_sensitivity_dBm, Ptmax_BS_dBm, Ptmin_BS_dBm, Ptmax_MS_dBm, Ptmin_MS_dBm, PCmargin_dB, delta)
% Apply Power Control Algorithm based on Delta
A_dB = Pt_dBm - Pr_dBm;
PrcDL_dBm = MS_sensitivity_dBm + PCmargin_dB;
PrcUL_dBm = BS_sensitivity_dBm + PCmargin_dB;
switch link_type
    case 1   % Downlink
        Ptout_dBm = PrcDL_dBm + A_dB*delta;
        Prout_dBm = PrcDL_dBm + A_dB*(delta-1);
        if(Ptout_dBm >= Ptmax_BS_dBm)
            Ptout_dBm = Ptmax_BS_dBm;
            Prout_dBm = Ptout_dBm - A_dB;
        elseif(Ptout_dBm <= Ptmin_BS_dBm)
            Ptout_dBm = Ptmin_BS_dBm;
            Prout_dBm = Ptout_dBm - A_dB;
        end
    case 2   % Uplink
        Ptout_dBm = PrcUL_dBm + A_dB*delta;
        Prout_dBm = PrcUL_dBm + A_dB*(delta-1);
        if(Ptout_dBm >= Ptmax_MS_dBm)
            Ptout_dBm = Ptmax_MS_dBm;
            Prout_dBm = Ptout_dBm - A_dB;
        elseif(Ptout_dBm <= Ptmin_MS_dBm)
            Ptout_dBm = Ptmin_MS_dBm;
            Prout_dBm = Ptout_dBm - A_dB;
        end
end
end