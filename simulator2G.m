%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%             2G SIMULATOR              %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all    % Close all open figures
clear        % Reset variables
clc          % Clear the command window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulator Parameters
tic

% Antennas Parameters
Ptmax_BS = 10;                                     % BS Max TX Power (Watts)
Ptmin_BS = 2;                                      % BS Min TX Power (Watts)
Ptmax_BS_dBm = (10*log10(Ptmax_BS))+30;            % BS Max TX Power (dBm)
Ptmin_BS_dBm = (10*log10(Ptmin_BS))+30;            % BS Min TX Power (dBm)

Ptmax_MS = 2;                                      % MS Max TX Power (Watts)
Ptmin_MS = 0.002;                                  % MS Min TX Power (Watts)
Ptmax_MS_dBm = (10*log10(Ptmax_MS))+30;            % MS Max TX Power (dBm)
Ptmin_MS_dBm = (10*log10(Ptmin_MS))+30;            % MS Min TX Power (dBm)

Prmin_BS_dBm = -104;                               % BS Sensitivity (dBm)
Prmin_MS_dBm = -102;                               % MS Sensitivity (dBm)

hBS = 30;                                          % BS height (meters)
fc = 900;                                          % Carrier Frequency (MHz)

% Propagation Parameters
sigmadB = 8;                                       % Shadowing St. Dev
LP0 = 0.99;                                        % Location Probability
Mf_dB = sigmadB*sqrt(2)*erfinv(2*LP0-1);           % Shadowing (Slow fading) Margin (dB)
Lmax = Ptmax_MS_dBm - (Prmin_BS_dBm + Mf_dB);      % Maximum Path Loss (dB)

% Original Okumura Model -> Cell Radius (m)
R = round((10^((Lmax-69.55-26.16*log10(fc)+13.82*log10(hBS))/(44.9-6.55*log10(hBS))))*1000);

% Network Parameters
p_DL = 0.45;                                       % Probability of Downlink State
p_UL = 0.45;                                       % Probability of Uplink State
p_IN = 0.1;                                        % Probability of Inactive State


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BS Deployment

filenameBS = 'BS_K=7.txt';                         % File with BS coordinates
[X_BS,Y_BS] = BSread(filenameBS);                  % Read from file
BSC = [X_BS,Y_BS];                                 % BS Coordinates
N_BS = length(BSC);                                % Number of BS (considering 2 interfering tiers)
if(strcmp(filenameBS,'BS_K=7.txt'))
    K = 7;                                         % Cluster Size
    N = 19;                                        % Normalizing factor (plot)
else
    
end
r = (1/sqrt(3))/N;                                 % Cell radius in the plot
Scale = R/r;                                       % Scaling Factor

plotBS(BSC(:,1),BSC(:,2),r,K);                     % Plot BS Deployment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mobile Station Deployment

N_MSe = 10000;                                     % Estimated Number of MS in the service area
Asquare = (0.85-0.15)^2;                           % Area of the square
a = r*sqrt(3)/2;                                   % Apothem of the hexagon
Aservice = (((6*r)*a)/2)*N_BS;                     % Area of the service area
N_MStot = round(N_MSe*(Asquare/Aservice));         % Estimated Number of MS in the square area

% MS coordinates generation
[X_MS,Y_MS] = uniformMS(0.15,0.85,0.15,0.85,N_MStot);

% Assign nearest BS to each MS
cellID = dsearchn([X_BS,Y_BS],delaunayn([X_BS,Y_BS]),[X_MS,Y_MS]);

% Compute distance between MS and nearest BS
distance = computeDistance(BSC(cellID,:),[X_MS, Y_MS]);

MSCtemp = [X_MS, Y_MS, distance*Scale];            % MS temporary coordinates and distance (m)
index = find(MSCtemp(:,3) < R);                    % Index of MSCtemp with distance less than the cell radius
MSC = MSCtemp(index,:);                            % Save in MSC only the MS inside the service area
N_MS = length(MSC);                                % Real Number of MS in the service area

plotMS(MSC(:,1),MSC(:,2));                         % Plot MS Deployment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mobile Station Deployment





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

