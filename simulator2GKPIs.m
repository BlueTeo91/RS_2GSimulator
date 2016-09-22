%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%             2G SIMULATOR              %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all    % Close all open figures
clear        % Reset variables
clc          % Clear the command window
tic          % Start stopwatch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulator Parameters

% Number of Snapshots
snapshots = 1;

% Number of Snapshots between two new MS deployments
MS_update = 10;

% Antennas Parameters
Ptmax_BS = 10;                                     % BS Max TX Power (Watts)
Ptmin_BS = 0.01;                                   % BS Min TX Power (Watts)
Ptmax_BS_dBm = (10*log10(Ptmax_BS))+30;            % BS Max TX Power (dBm)
Ptmin_BS_dBm = (10*log10(Ptmin_BS))+30;            % BS Min TX Power (dBm)

Ptmax_MS = 2;                                      % MS Max TX Power (Watts)
Ptmin_MS = 0.00002;                                % MS Min TX Power (Watts)
Ptmax_MS_dBm = (10*log10(Ptmax_MS))+30;            % MS Max TX Power (dBm)
Ptmin_MS_dBm = (10*log10(Ptmin_MS))+30;            % MS Min TX Power (dBm)

Prmin_BS_dBm = -104;                               % BS Sensitivity (dBm)
Prmin_MS_dBm = -102;                               % MS Sensitivity (dBm)

hBS = 30;                                          % BS height (meters)
fc = 1800;                                         % Carrier Frequency (MHz)

% Noise Figures
F_BS_dB = 5;                                       % BS noise figure (dB)
F_MS_dB = 10;                                      % MS noise figure (dB)

F_BS = 10^(F_BS_dB/10);                            % BS noise figure
F_MS = 10^(F_MS_dB/10);                            % MS noise figure

% Propagation Parameters
sigmadB = 8;                                       % Shadowing St. Dev
Pcov = 0.95;                                       % Coverage Probability
Mf_dB = sigmadB*sqrt(2)*erfinv(2*Pcov-1);          % Shadowing (Slow fading) Margin (dB)

% Maximum Path Loss (dB)
Lmax = min(Ptmax_MS_dBm - (Prmin_BS_dBm + Mf_dB),Ptmax_BS_dBm - (Prmin_MS_dBm + Mf_dB));

% Original Okumura Model -> Cell Radius (meters)
R = round((10^((Lmax-69.55-26.16*log10(fc)+13.82*log10(hBS))/(44.9-6.55*log10(hBS))))*1000);

% Network Parameters
K = 7;                                             % Cluster Size (3 or 7)
N_MSe = 12750;                                     % Estimated Number of MS in the service area

Pcall_average = 1.0;                               % Average call probability
Pcall_StDev = 0.00;                                % Call probability standard deviation

p_DL = 0.5;                                        % Probability of Downlink State
p_UL = 0.5;                                        % Probability of Uplink State

Rb = 271e3;                                        % Bitrate (bit/s)

% Directed Retry
retry = 1;
if(retry==1)
    retryS = 'ON';
else
    retryS = 'OFF';
end

% Power Control Parameters
PCmargin_dB = 10;                                  % Power Control Margin (dB)
delta = 1;                                         % Delta [0,1]

% Total number of Radio Resource Units available to the operator
N_RU = 700;

% Outage Thresholds
SNR_Out_Thr_DL = computeSNR(Prmin_MS_dBm,F_MS,Rb) + 3;    % Outage SNR Downlink (dB)
SNR_Out_Thr_UL = computeSNR(Prmin_BS_dBm,F_BS,Rb) + 3;    % Outage SNR Uplink (dB)
SIR_Out_Thr = 10;                                         % Outage SIR (dB)

% Forced Termination Thresholds
SNR_FT_Thr_DL = computeSNR(Prmin_MS_dBm,F_MS,Rb) - 3;     % FT SNR Downlink (dB)
SNR_FT_Thr_UL = computeSNR(Prmin_BS_dBm,F_BS,Rb) - 3;     % FT SNR Uplink (dB)
SIR_FT_Thr = 5;                                           % FT SIR (dB)

% KPIs initialization
network_load_th_TOT = 0;
network_load_TOT = 0;
refCell_load_TOT = 0;
Avg_retries_TOT = 0;
retry_rate_TOT = 0;
blocking_rate_TOT = 0;
refCell_blocking_rate_TOT = 0;
outage_rate_TOT = 0;
forced_termination_rate_TOT = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base Stations Deployment

if(K==3)
    filenameBS = 'BS_K=3.txt';                     % File with BS coordinates (K=3)
    clustersize = '3';
elseif(K==7)
    filenameBS = 'BS_K=7.txt';                     % File with BS coordinates (K=7)
    clustersize = '7';
end

[X_BS,Y_BS] = BSread(filenameBS);                  % Read from file
BSC = [X_BS,Y_BS];                                 % BS Coordinates
N_BS = length(BSC(:,1));                           % Number of BS (considering 2 interfering tiers)
N = 19;                                            % Normalizing factor (plot)
r = (1/sqrt(3))/N;                                 % Cell radius in the plot
scale = R/r;                                       % Scaling Factor

%plotBS(BSC(:,1),BSC(:,2),r,K);                     % Plot BS Deployment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mobile Stations Deployment

for snap = 1:snapshots
    
    % Print to video
    fprintf('Snapshot %d/%d\n',snap,snapshots);
    
    if ((snap == 1) || (mod(snap,MS_update) == 0))
        
        % Print to video
        fprintf('MS Deployment\n');
        
        % Area of the rectangle (square km)
        if(K==3)
            Arect = (((0.7597-0.2887)*(0.7368-0.3158))*(scale^2))/(1e6);
        elseif(K==7)
            Arect = (((0.8508-0.1519)*(0.8421-0.1579))*(scale^2))/(1e6);
        end
        a = r*sqrt(3)/2;                                   % Apothem of the hexagon (plot)
        Aservice = (((((6*r)*a)/2)*N_BS)*(scale^2))/(1e6); % Area of the service area (square km)
        N_MStot = round(N_MSe*(Arect/Aservice));           % Number of MS in the rectangular area
        
        % MS coordinates generation
        if(K==3)
            [X_MS,Y_MS] = uniformMS(0.2887,0.7597,0.3158,0.7368,N_MStot);
        elseif(K==7)
            [X_MS,Y_MS] = uniformMS(0.1519,0.8508,0.1579,0.8421,N_MStot);
        end
        
        % Assign nearest BS to each MS and compute distance
        [cellID,shortest_distance] = dsearchn([X_BS,Y_BS],delaunayn([X_BS,Y_BS]),[X_MS,Y_MS]);
        
        % MS temporary matrix
        MSCtemp = [X_MS, Y_MS, cellID, shortest_distance*scale];
        MSindex = find(MSCtemp(:,4) < R);                  % Index of MSCtemp with distance less than cell radius
        MSC = MSCtemp(MSindex,1:2);                        % Save in MSC only the MS inside the service area
        N_MS = length(MSC(:,1));                           % Real Number of MS in the service area
        
        % Compute distance between each MS and each BS
        distance = zeros(N_MS,N_BS);
        for i = 1:N_MS
            for j = 1:N_BS
                distance(i,j) = (computeDistance(MSC(i,:),BSC(j,:)).*scale);
            end
        end
        
        % Sort distance matrix in ascending order
        [neighbouring_distance, neighbouring_index] = sort(distance,2);
        
        % Add neighbouring cell ID and distances to MS matrix
        MSC = [MSC(:,1:2), neighbouring_index(:,1:4), neighbouring_distance(:,1:4)];
        
        %plotMS(MSC(:,1),MSC(:,2));                         % Plot MS Deployment
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Traffic Generation
    
    % Probability of being in a call between [PcallMin, PcallMax]
    Pcall = abs(Pcall_average+Pcall_StDev*randn(1,1));
    
    % Generate random vector of 0s and 1s -> 1 = calling, 0 = not calling
    % (with a percentage of Pcall users calling)
    calling = rand(N_MS,1) <= Pcall;
    N_MScalling = sum(calling);                        % Number of calling MS
    MSC = [MSC(:,1:10), calling, zeros(N_MS,1)];       % Add call state column to MS matrix
    
    % Assign traffic type based on probabilities defined above
    % 1 = DOWNLINK
    % 2 = UPLINK
    traffic_kind = randsample(2,N_MScalling,true,[p_DL p_UL]);
    j = 1:N_MScalling;
    i = find(calling);
    MSC(i,12) = traffic_kind(j);                       % Add traffic type column to MS matrix
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Propagation
    
    % Inizialize first column with traffic type
    links = [];                                        % Inizialization of links for a new snapshot
    links(:,1) = MSC(:,12);
    temp_power = zeros(N_MS,4);
    
    % Compute received power for each link (using Beacon Signal)
    for i=1:N_MS
        switch links(i,1)
            case 1      % Downlink
                temp_power(i,:) = propagation(Ptmax_BS_dBm,fc,hBS,sigmadB,MSC(i,7:10));
            case 2      % Uplink
                temp_power(i,:) = propagation(Ptmax_MS_dBm,fc,hBS,sigmadB,MSC(i,7:10));
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
    links = [links(:,1), temp_BSid, temp_power];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% RUs Assignment
    
    N_RU_cellDL = round((N_RU/K)/2);           % Number of RRUs for DL per cell
    N_RU_cellUL = round((N_RU/K)/2);           % Number of RRUs for UL per cell
    
    % Add number of available RUs DL and RUs UL to BSC matrix
    BSC = [BSC(:,1:2), N_RU_cellDL*ones(N_BS,1), N_RU_cellUL*ones(N_BS,1)];
    N_RU_tot = sum(BSC(:,3)) + sum(BSC(:,4));
    N_RU_cell = N_RU_cellDL + N_RU_cellUL;
    
    % Admission Control - Directed Retry
    
    % Scheduling: First-Come-First-Served (FCFS)
    % MSid ordered by arrival
    
    % Links column 10:
    %    0 -> RU not requested
    %   -1 -> RU not assigned (blocked)
    % BSid -> RU assigned by BSid
    
    % Links column 11:
    %      0 -> RU not assigned
    % RUs ID -> number of the RU
    
    % Links column 12:
    % 0 -> MS not connected
    % 1 -> MS connected to the first best server
    % 2 -> MS connected to the second best server
    % 3 -> MS connected to the third best server
    % 4 -> MS connected to the fourth best server
    
    links = [links(:,1:9), zeros(N_MS,1), zeros(N_MS,1), zeros(N_MS,1)];
    for i = 1:N_MS
        switch links(i,1)
            case 1      % Downlink
                if(retry==0)
                    if((links(i,6)>Prmin_MS_dBm) && ((BSC(links(i,2),3)) > 0))
                        BSC(links(i,2),3) = BSC(links(i,2),3) - 1;
                        links(i,10) = links(i,2);
                        links(i,11) = N_RU_cellDL - BSC(links(i,2),3);
                        links(i,12) = 1;
                    else
                        links(i,10) = -1;   % Blocked MS (No RU assigned)
                    end
                else
                    if((links(i,6)>Prmin_MS_dBm) && ((BSC(links(i,2),3)) > 0))
                        BSC(links(i,2),3) = BSC(links(i,2),3) - 1;
                        links(i,10) = links(i,2);
                        links(i,11) = N_RU_cellDL - BSC(links(i,2),3);
                        links(i,12) = 1;
                    elseif((links(i,7)>Prmin_MS_dBm) && ((BSC(links(i,3),3)) > 0))
                        BSC(links(i,3),3) = BSC(links(i,3),3) - 1;
                        links(i,10) = links(i,3);
                        links(i,11) = N_RU_cellDL - BSC(links(i,3),3);
                        links(i,12) = 2;
                    elseif((links(i,8)>Prmin_MS_dBm) && ((BSC(links(i,4),3)) > 0))
                        BSC(links(i,4),3) = BSC(links(i,4),3) - 1;
                        links(i,10) = links(i,4);
                        links(i,11) = N_RU_cellDL - BSC(links(i,4),3);
                        links(i,12) = 3;
                    elseif((links(i,9)>Prmin_MS_dBm) && ((BSC(links(i,5),3)) > 0))
                        BSC(links(i,5),3) = BSC(links(i,5),3) - 1;
                        links(i,10) = links(i,5);
                        links(i,11) = N_RU_cellDL - BSC(links(i,5),3);
                        links(i,12) = 4;
                    else
                        links(i,10) = -1;    % Blocked MS (No RU assigned)
                    end
                end
            case 2      % Uplink
                if(retry==0)
                    if((links(i,6)>Prmin_BS_dBm) && ((BSC(links(i,2),4)) > 0))
                        BSC(links(i,2),4) = BSC(links(i,2),4) - 1;
                        links(i,10) = links(i,2);
                        links(i,11) = N_RU_cellUL - BSC(links(i,2),4);
                        links(i,12) = 1;
                    else
                        links(i,10) = -1;   % Blocked MS (No RU assigned)
                    end
                else
                    if((links(i,6)>Prmin_BS_dBm) && ((BSC(links(i,2),4)) > 0))
                        BSC(links(i,2),4) = BSC(links(i,2),4) - 1;
                        links(i,10) = links(i,2);
                        links(i,11) = N_RU_cellUL - BSC(links(i,2),4);
                        links(i,12) = 1;
                    elseif((links(i,7)>Prmin_BS_dBm) && ((BSC(links(i,3),4)) > 0))
                        BSC(links(i,3),4) = BSC(links(i,3),4) - 1;
                        links(i,10) = links(i,3);
                        links(i,11) = N_RU_cellUL - BSC(links(i,3),4);
                        links(i,12) = 2;
                    elseif((links(i,8)>Prmin_BS_dBm) && ((BSC(links(i,4),4)) > 0))
                        BSC(links(i,4),4) = BSC(links(i,4),4) - 1;
                        links(i,10) = links(i,4);
                        links(i,11) = N_RU_cellUL - BSC(links(i,4),4);
                        links(i,12) = 3;
                    elseif((links(i,9)>Prmin_BS_dBm) && ((BSC(links(i,5),4)) > 0))
                        BSC(links(i,5),4) = BSC(links(i,5),4) - 1;
                        links(i,10) = links(i,5);
                        links(i,11) = N_RU_cellUL - BSC(links(i,5),4);
                        links(i,12) = 4;
                    else
                        links(i,10) = -1;    % Blocked MS (No RU assigned)
                    end
                end
            otherwise   % Inactive (Not Calling)
                links(i,10) = 0;         % RU not requested
        end
    end
    
    N_RU_available = sum(BSC(:,3)) + sum(BSC(:,4));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Power Control
    
    connected_MSid = find(links(:,11)>0);
    Pt_dBm = zeros(N_MS,1);
    Pr_dBm = zeros(N_MS,1);
    for i = 1:N_MS
        switch links(i,12)
            case 1      % MS connected to the first best server
                if(links(i,1)==1)   % Downlink
                    Pt_dBm(i,1) = Ptmax_BS_dBm;
                else                % Uplink
                    Pt_dBm(i,1) = Ptmax_MS_dBm;
                end
                Pr_dBm(i,1) = links(i,6);
            case 2      % MS connected to the second best server
                if(links(i,1)==1)   % Downlink
                    Pt_dBm(i,1) = Ptmax_BS_dBm;
                else                % Uplink
                    Pt_dBm(i,1) = Ptmax_MS_dBm;
                end
                Pr_dBm(i,1) = links(i,7);
            case 3      % MS connected to the third best server
                if(links(i,1)==1)   % Downlink
                    Pt_dBm(i,1) = Ptmax_BS_dBm;
                else                % Uplink
                    Pt_dBm(i,1) = Ptmax_MS_dBm;
                end
                Pr_dBm(i,1) = links(i,8);
            case 4      % MS connected to the fourth best server
                if(links(i,1)==1)   % Downlink
                    Pt_dBm(i,1) = Ptmax_BS_dBm;
                else                % Uplink
                    Pt_dBm(i,1) = Ptmax_MS_dBm;
                end
                Pr_dBm(i,1) = links(i,9);
        end
    end
    
    connection_type = links(connected_MSid,1);
    connected_BSid = links(connected_MSid,10);
    RUid = links(connected_MSid,11);
    Pt_dBm = Pt_dBm(connected_MSid,1);
    Pr_dBm = Pr_dBm(connected_MSid,1);
    
    % Signal Based Power Control (PC)
    PtPC_dBm = zeros(length(connected_MSid),1);
    PrPC_dBm = zeros(length(connected_MSid),1);
    for i=1:length(connected_MSid)
        [PtPC_dBm(i,1), PrPC_dBm(i,1)] = powercontrol(connection_type(i,1),Pt_dBm(i,1),Pr_dBm(i,1),Prmin_BS_dBm,Prmin_MS_dBm,Ptmax_BS_dBm,Ptmin_BS_dBm,Ptmax_MS_dBm,Ptmin_MS_dBm,PCmargin_dB,delta);
    end
    
    % Create connected_links matrix with the following columns
    connected_links = [connection_type, connected_MSid, connected_BSid, RUid, PtPC_dBm, PrPC_dBm];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SNR Computation
    
    % Compute SNR for each connected link
    SNR_dB = zeros(length(connected_MSid),1);
    for i=1:length(connected_MSid)
        switch connected_links(i,1)
            case 1      % Downlink
                SNR_dB(i,1) = computeSNR(connected_links(i,6),F_MS,Rb);
            case 2      % Uplink
                SNR_dB(i,1) = computeSNR(connected_links(i,6),F_BS,Rb);
        end
    end
    
    % Add SNR (dB) column in connected_links matrix
    connected_links = [connected_links(:,1:6), SNR_dB];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SIR Computation
    
    % Add state column
    % Active = 1
    % Silent = 0
    state = rand(length(connected_links(:,1)),1) <= 0.45;
    connected_links = [connected_links(:,1:7), state];
    
    % Search cell ID of two interfering tiers
    interfering_cellid(:,1) = K+1:K:N_BS;
    
    % SIR DOWNLINK
    % Create MS_DLrefCell matrix: [connected_MSid, RUid, SIR_dB, SNR_dB]
    MS_DLrefCellindex = find(connected_links(:,3)==1 & connected_links(:,1)==1 & connected_links(:,8)==1);
    MS_DLrefCell = [connected_links(MS_DLrefCellindex,2),connected_links(MS_DLrefCellindex,4),zeros(length(MS_DLrefCellindex),1),zeros(length(MS_DLrefCellindex),1)];
    
    for i = 1:length(MS_DLrefCellindex)
        I = 0;
        for j = 1:length(interfering_cellid)
            interferinglinks_index = find(connected_links(:,3)==interfering_cellid(j,1) & connected_links(:,1)==1 & connected_links(:,4)==MS_DLrefCell(i,2) & connected_links(:,8)==1);
        end
        for k = 1:length(interferinglinks_index)
            interferingDL_distance = computeDistance(MSC(MS_DLrefCell(i,1),1:2),BSC(connected_links(interferinglinks_index(k,1)),1:2));
            I_dBm = propagation(connected_links(interferinglinks_index(k,1),5),fc,hBS,sigmadB,interferingDL_distance*scale);
            I = I + 10^((I_dBm-30)/10);
        end
        I_dBm = 10*log10(I);
        C_dBm = connected_links(MS_DLrefCellindex(i,1),6);
        SIR_dB = C_dBm - I_dBm;
        MS_DLrefCell(i,3) = SIR_dB;
        MS_DLrefCell(i,4) = connected_links(MS_DLrefCellindex(i,1),7);
    end
    
    % SIR UPLINK
    % Create MS_ULrefCell matrix: [connected_MSid, RUid, SIR_dB, SNR_dB]
    MS_ULrefCellindex = find(connected_links(:,3)==1 & connected_links(:,1)==2 & connected_links(:,8)==1);
    MS_ULrefCell = [connected_links(MS_ULrefCellindex,2),connected_links(MS_ULrefCellindex,4),zeros(length(MS_ULrefCellindex),1),zeros(length(MS_ULrefCellindex),1)];
    
    for i = 1:length(MS_ULrefCellindex)
        I = 0;
        for j = 1:length(interfering_cellid)
            interferinglinks_index = find(connected_links(:,3)==interfering_cellid(j,1) & connected_links(:,1)==2 & connected_links(:,4)==MS_ULrefCell(i,2) & connected_links(:,8)==1);
        end
        for k = 1:length(interferinglinks_index)
            interferingUL_distance = computeDistance(BSC(1,1:2),MSC(connected_links(interferinglinks_index(k,1),2),1:2));
            I_dBm = propagation(connected_links(interferinglinks_index(k,1),5),fc,hBS,sigmadB,interferingUL_distance*scale);
            I = I + 10^((I_dBm-30)/10);
        end
        I_dBm = 10*log10(I);
        C_dBm = connected_links(MS_ULrefCellindex(i,1),6);
        SIR_dB = C_dBm - I_dBm;
        MS_ULrefCell(i,3) = SIR_dB;
        MS_ULrefCell(i,4) = connected_links(MS_ULrefCellindex(i,1),7);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% KPIs Computation
    
    % Network Load Percentage (Theoretical)
    network_load_th = (N_MS/N_RU_tot)*100;
    % Incremental sum of network_load for KPI computation
    network_load_th_TOT = network_load_th_TOT + network_load_th;
    
    % Network Load Percentage (based on assigned RUs - Effective)
    network_load = ((N_RU_tot - N_RU_available)/N_RU_tot)*100;
    % Incremental sum of network_load for KPI computation
    network_load_TOT = network_load_TOT + network_load;
    
    % Reference Cell Load Percentage (BSid=1 - Effective)
    refCell_load = ((N_RU_cell - (BSC(1,3) + BSC(1,4)))/N_RU_cell)*100;
    % Incremental sum of refCell_load for KPI computation
    refCell_load_TOT = refCell_load_TOT + refCell_load;
    
    % Remove border effect
    if(K==3)
        border_cellsID = [23,27,26,30,29,33,32,31,35,34,38,37,41,40,44,43,45,46,48,49,51,52,54,55,57,56,24];
    elseif(K==7)
        border_cellsID = [51,52,63,58,59,60,65,66,67,72,73,74,75,80,81,82,87,88,89,90,95,96,97,102,103,104,105,110,11,112,117,118,119,114,125,126,121,132,133,128,129,56];
    end
    
    for i = 1:length(border_cellsID)
        border_links_index = find(links(:,10)==border_cellsID(i));
        links(border_links_index,:) = [];
    end
    
    % Remove not requested links
    inactive_links_index = find(links(:,10)==0);
    links(inactive_links_index,:) = [];
    
    % Average number of retries per user
    retries = links(:,12) - 1;
    if(retry==0)
        retries(retries==-1)=0;
    else
        retries(retries==-1)=3; 
    end
    Avg_retries = sum(retries)/length(retries);
    Avg_retries_TOT = Avg_retries_TOT + Avg_retries;
    
    % Retry rate
    retry_rate = (length(find(retries>0))/length(retries))*100;
    retry_rate_TOT = retry_rate_TOT + retry_rate;
    
    % Blocking Rate (%)
    N_MSblocked = length(find(links(:,10)==-1));
    N_MSdownlink = length(find(links(:,1)==1));
    N_MSuplink = length(find(links(:,1)==2));
    blocking_rate = (N_MSblocked/(N_MSdownlink + N_MSuplink))*100;
    % Incremental sum of blocking_rate for KPI computation
    blocking_rate_TOT = blocking_rate_TOT + blocking_rate;
    
    % Reference Cell Blocking Rate (%) (BSid=1)
    if(retry==0)
        row = find(links(:,2)==1);
        index_temp = sort(row);
    else        
        [row,column] = find(links(:,2:5)==1);
        index_temp = sort(row);
    end
    RU_temp = links(index_temp,10);
    blocked_index_refCell = find(RU_temp(:,1)==-1);
    N_MSblocked_refCell = length(blocked_index_refCell);
    refCell_blocking_rate = (N_MSblocked_refCell / (N_MSblocked_refCell + (N_RU_cell - (BSC(1,3) + BSC(1,4)))))*100;
    % Incremental sum of refCell_blocking_rate for KPI computation
    refCell_blocking_rate_TOT = refCell_blocking_rate_TOT + refCell_blocking_rate;
    
    % Outage Rate (%) (BSid=1)
    N_MS_Out_DL = 0;
    for i = 1:length(MS_DLrefCell(:,1))
        if(MS_DLrefCell(i,4)<SNR_Out_Thr_DL || MS_DLrefCell(i,3)<SIR_Out_Thr)
            N_MS_Out_DL = N_MS_Out_DL+1;
        end
    end
    N_MS_Out_UL = 0;
    for i = 1:length(MS_ULrefCell(:,1))
        if(MS_ULrefCell(i,4)<SNR_Out_Thr_UL || MS_ULrefCell(i,3)<SIR_Out_Thr)
            N_MS_Out_UL = N_MS_Out_UL+1;
        end
    end
    
    outage_rate = ((N_MS_Out_DL+N_MS_Out_UL)/(length(MS_DLrefCell(:,1))+length(MS_ULrefCell(:,1))))*100;
    % Incremental sum of outage_rate for KPI computation
    outage_rate_TOT = outage_rate_TOT + outage_rate;
    
    % Forced Termination Rate (%) (BSid=1)
    N_MS_FT_DL = 0;
    for i = 1:length(MS_DLrefCell(:,1))
        if(MS_DLrefCell(i,4)<SNR_FT_Thr_DL || MS_DLrefCell(i,3)<SIR_FT_Thr)
            N_MS_FT_DL = N_MS_FT_DL+1;
        end
    end
    N_MS_FT_UL = 0;
    for i = 1:length(MS_ULrefCell(:,1))
        if(MS_ULrefCell(i,4)<SNR_FT_Thr_UL || MS_ULrefCell(i,3)<SIR_FT_Thr)
            N_MS_FT_UL = N_MS_FT_UL+1;
        end
    end
    
    forced_termination_rate = ((N_MS_FT_DL+N_MS_FT_UL)/(length(MS_DLrefCell(:,1))+length(MS_ULrefCell(:,1))))*100;
    % Incremental sum of outage_rate for KPI computation
    forced_termination_rate_TOT = forced_termination_rate_TOT + forced_termination_rate;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print SIR_dB and SNR_dB in a file

% fileID_SIR_DL = fopen('SIR_DL.txt','at');
% fprintf(fileID_SIR_DL,'%.4f\n',MS_DLrefCell(:,3));
% fclose(fileID_SIR_DL);
% 
% fileID_SNR_DL = fopen('SNR_DL.txt','at');
% fprintf(fileID_SNR_DL,'%.4f\n',MS_DLrefCell(:,4));
% fclose(fileID_SNR_DL);
% 
% fileID_SIR_UL = fopen('SIR_UL.txt','at');
% fprintf(fileID_SIR_UL,'%.4f\n',MS_ULrefCell(:,3));
% fclose(fileID_SIR_UL);
% 
% fileID_SNR_UL = fopen('SNR_UL.txt','at');
% fprintf(fileID_SNR_UL,'%.4f\n',MS_ULrefCell(:,4));
% fclose(fileID_SNR_UL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total KPIs computation & Print to file

network_load_th_TOT = network_load_th_TOT / snapshots;
network_load_TOT = network_load_TOT / snapshots;
refCell_load_TOT = refCell_load_TOT / snapshots;
Avg_retries_TOT = Avg_retries_TOT / snapshots;
retry_rate_TOT = retry_rate_TOT / snapshots;
blocking_rate_TOT = blocking_rate_TOT / snapshots;
refCell_blocking_rate_TOT = refCell_blocking_rate_TOT / snapshots;
outage_rate_TOT = outage_rate_TOT / snapshots;
forced_termination_rate_TOT = forced_termination_rate_TOT / snapshots;

filename = ['KPIs_K=' clustersize '.txt'];
fileID = fopen(filename,'at');
fprintf(fileID,'%.4f\t%.4f\t%.4f\t%.4f\n',delta,PCmargin_dB,outage_rate_TOT,forced_termination_rate_TOT);
fclose('all');

% Print to video
fprintf('Directed Retry: %s\nPC Parameters: Delta = %.4f\tMargin = %d dB\nNetwork Load: %.4f (Effective)\t%.4f (Theoretical)\nReference Cell Load: %.4f (Effective)\nRetry Rate: %.4f\nBlocking Rate: %.4f\nReference Cell Blocking Rate: %.4f\nOutage Rate: %.4f\nForced Termination Rate: %.4f\n',retryS,delta,PCmargin_dB,network_load_TOT,network_load_th_TOT,refCell_load_TOT,retry_rate_TOT,blocking_rate_TOT,refCell_blocking_rate_TOT,outage_rate_TOT,forced_termination_rate_TOT);

toc          % Stop stopwatch