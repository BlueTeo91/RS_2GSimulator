function Pr = propagation(Pt,fc,hBS,sigmadB,d)
% Original Okumura Hata Formula + Shadowing
% Pt is the transmitted power depending on traffic kind
L = 69.55 + 26.16*log10(fc) - 13.82*log10(hBS) + (44.9 - 6.55*log10(hBS))*log10(d);
shadowing = sigmadB*randn(1,4);     % Generation of 4 shadowing samples (for each BS)
Pr = Pt - (L + shadowing);          % Received Power dB
end