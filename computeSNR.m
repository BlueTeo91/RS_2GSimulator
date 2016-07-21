function SNR_dB = computeSNR(Pr_dBm,F,Rb)
% Compute SNR
% Compute System Equivalent Noise Temperature (°K)
Tsyst = 290+290*(F-1);

% Monolateral noise spectral density (Watts)
No = ((1.38e-23)*Tsyst)/2;

% Compute Signal-to-Noise Ratio (dB)
SNR_dB = Pr_dBm - (10*log10(2*No*Rb)+30);
end