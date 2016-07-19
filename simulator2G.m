%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%             2G SIMULATOR              %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all    % Close all open figures
clear        % Reset variables
clc          % Clear the command window
tic          % Start stopwatch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulator Parameters

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

% Noise Figures
F_BS_dB = 5;                                       % BS noise figure (dB)
F_MS_dB = 10;                                      % MS noise figure (dB)

F_BS = 10^(F_BS_dB/10);                            % BS noise figure
F_MS = 10^(F_MS_dB/10);                            % MS noise figure

% Propagation Parameters
sigmadB = 8;                                       % Shadowing St. Dev
LP0 = 0.99;                                        % Location Probability
Mf_dB = sigmadB*sqrt(2)*erfinv(2*LP0-1);           % Shadowing (Slow fading) Margin (dB)
Lmax = Ptmax_MS_dBm - (Prmin_BS_dBm + Mf_dB);      % Maximum Path Loss (dB)

% Original Okumura Model -> Cell Radius (meters)
R = round((10^((Lmax-69.55-26.16*log10(fc)+13.82*log10(hBS))/(44.9-6.55*log10(hBS))))*1000);

% Network Parameters
N_MSe = 10000;                                     % Estimated Number of MS in the service area
Pcall_average = 0.5;                               % Average call probability
Pcall_StDev = 0.05;                                % Call probability standard deviation
p_DL = 0.45;                                       % Probability of Downlink State
p_UL = 0.45;                                       % Probability of Uplink State
p_IN = 0.1;                                        % Probability of Inactive State
Rb = 271e3;                                        % Bitrate (bit/s)

% Total number of Radio Resource Units available to the operator
N_RU = 700;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base Stations Deployment

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

%plotBS(BSC(:,1),BSC(:,2),r,K);                     % Plot BS Deployment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mobile Stations Deployment

% Area of the rectangle (square km)
Arect = (((0.8508-0.1519)*(0.8421-0.1579))*(Scale^2))/(1e6);
a = r*sqrt(3)/2;                                   % Apothem of the hexagon (plot)
Aservice = (((((6*r)*a)/2)*N_BS)*(Scale^2))/(1e6); % Area of the service area (square km)
N_MStot = round(N_MSe*(Arect/Aservice));           % Number of MS in the rectangular area

% MS coordinates generation
[X_MS,Y_MS] = uniformMS(0.1519,0.8508,0.1579,0.8421,N_MStot);

% Assign nearest BS to each MS and compute distance
[cellID,shortest_distance] = dsearchn([X_BS,Y_BS],delaunayn([X_BS,Y_BS]),[X_MS,Y_MS]);

% MS temporary matrix
MSCtemp = [X_MS, Y_MS, cellID, shortest_distance*Scale];
MSindex = find(MSCtemp(:,4) < R);                  % Index of MSCtemp with distance less than cell radius
MSC = MSCtemp(MSindex,1:2);                        % Save in MSC only the MS inside the service area
N_MS = length(MSC);                                % Real Number of MS in the service area

% Compute distance between each MS and each BS
distance = zeros(N_MS,N_BS);
for i = 1:N_MS
    for j = 1:N_BS
        distance(i,j) = (computeDistance(MSC(i,:),BSC(j,:)).*Scale);
    end
end

% Sort distance matrix in ascending order
[neighbouring_distance, neighbouring_index] = sort(distance,2);

% Add neighbouring cell ID and distances to MS matrix
MSC = [MSC, neighbouring_index(:,1:4), neighbouring_distance(:,1:4)];

% Count the number of MS in each cell
N_MS_eachcell = zeros(N_BS,1);
for i = 1:N_BS
    N_MS_eachcell(i,1) = sum(MSC(:,3) == i);
end

%plotMS(MSC(:,1),MSC(:,2));                         % Plot MS Deployment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Traffic Generation

% Probability of being in a call between [PcallMin, PcallMax]
Pcall = abs(Pcall_average+Pcall_StDev*randn(1,1));

% Generate random vector of 0s and 1s -> 1 = calling, 0 = not calling
% (with a percentage of Pcall users calling)
calling = (rand(N_MS,1) <= Pcall);
N_MScalling = sum(calling);                        % Number of calling MS
MSC = [MSC, calling, zeros(N_MS,1)];               % Add call state column to MS matrix

% Count the number of MS in each cell starting a call
N_MS_calling_eachcell = zeros(N_BS,1);
for i = 1:N_BS
    for k = 1:N_MS
        if(MSC(k,3) == i && MSC(k,11) == 1)
            N_MS_calling_eachcell(i,1) = N_MS_calling_eachcell(i,1) + 1;
        end
    end
end

% Assign traffic type based on probabilities defined above
% 1 = DOWNLINK
% 2 = UPLINK
% 3 = SILENT
traffic_kind = randsample(3,N_MScalling,true,[p_DL p_UL p_IN]);
j = 1:N_MScalling;
i = find(calling);
MSC(i,12) = traffic_kind(j);                       % Add traffic type column to MS matrix

% Count the number of MS in DOWNLINK state in each cell
N_MS_downlink_eachcell = zeros(N_BS,1);
for i = 1:N_BS
    for k = 1:N_MS
        if(MSC(k,3) == i && MSC(k,12) == 1)
            N_MS_downlink_eachcell(i,1) = N_MS_downlink_eachcell(i,1) + 1;
        end
    end
end

% Count the number of MS in UPLINK state in each cell
N_MS_uplink_eachcell = zeros(N_BS,1);
for i = 1:N_BS
    for k = 1:N_MS
        if(MSC(k,3) == i && MSC(k,12) == 2)
            N_MS_uplink_eachcell(i,1) = N_MS_uplink_eachcell(i,1) + 1;
        end
    end
end

% Count the number of MS in SILENT state in each cell
N_MS_inactive_eachcell = zeros(N_BS,1);
for i = 1:N_BS
    for k = 1:N_MS
        if(MSC(k,3) == i && MSC(k,12) == 3)
            N_MS_inactive_eachcell(i,1) = N_MS_inactive_eachcell(i,1) + 1;
        end
    end
end

% Add columns with number of MS in active state (DOWNLINK + UPLINK)
% BSC = [BSC, N_MS_downlink_eachcell + N_MS_uplink_eachcell];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Propagation

% Initialization (max values)
Pt_BS = Ptmax_BS_dBm;
Pt_MS = Ptmax_MS_dBm;

% Inizialize first column with traffic type
links(:,1) = MSC(:,12);
temp_power = zeros(N_MS,4);

% Compute received power for each link
for i=1:N_MS
    switch links(i,1)
        case 1      % Downlink
            temp_power(i,:) = propagation(Pt_BS,fc,hBS,sigmadB,MSC(i,7:10));
        case 2      % Uplink
            temp_power(i,:) = propagation(Pt_MS,fc,hBS,sigmadB,MSC(i,7:10));
        case 3      % Silent
            temp_power(i,:) = propagation(Pt_BS,fc,hBS,sigmadB,MSC(i,7:10));
        otherwise   % Inactive (Not Calling)
            temp_power(i,:) = propagation(Ptmax_BS_dBm,fc,hBS,sigmadB,MSC(i,7:10));
    end
end

% Sort power in descending order and correspondent BSid
temp_BSid = MSC(:,3:6);
[temp_power,power_index] = sort(temp_power,2,'descend');
for i = 1:N_MS
    temp_BSid(i,:) = temp_BSid(i,power_index(i,:));
end

% Add BSid and power ordered by received power to links matrix
links = [links, temp_BSid, temp_power];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SNR Computation

% Compute SNR for each link
SNR_dB = zeros(N_MS,4);
for i=1:N_MS
    switch links(i,1)
        case 1      % Downlink
            SNR_dB(i,:) = computeSNR(links(i,6:9),F_MS,Rb);
        case 2      % Uplink
            SNR_dB(i,:) = computeSNR(links(i,6:9),F_BS,Rb);
        case 3      % Silent
            SNR_dB(i,:) = computeSNR(links(i,6:9),F_MS,Rb);
        otherwise   % Inactive (Not Calling)
            SNR_dB(i,:) = computeSNR(links(i,6:9),F_MS,Rb);
    end
end

% Add SNR columns to links matrix
links = [links, SNR_dB];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUs Assignment

N_RU_cell = 700/K;                                 % Number of RRUs per cell

% Add number of available RUs to BSC matrix 
BSC = [BSC, N_RU_cell*ones(N_BS,1)];       




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIR Assignment





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc          % Stop stopwatch