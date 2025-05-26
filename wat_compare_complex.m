% Author: Tomasz Syrylo
%% Setup Parameters
% Carrier frequency (Hz) - Band 77 is around 4 GHz
fc = 4e9;
% FIRST CONFIG
% Base station parameters 
% bsPosition1 = [52.251387,20.905994];    % [Latitude, Longitude]
% bsArrayOrientation1 = [120 -15].';     % [Azimuth; Elevation] (degrees)
% bsPosition2 = [52.256474,20.901313];    % [Latitude, Longitude] 
% bsArrayOrientation2 = [-140 -15].';     % [Azimuth; Elevation] (degrees)
% bsPosition3 = [52.251076,20.899524];    % [Latitude, Longitude]
% bsArrayOrientation3 = [40 -15].';       % [Azimuth; Elevation] (degrees)

%(0 deg is East, 90 deg is North)

% SECOND CONFIG
bsPosition1 = [52.256474,20.901313];    % [Latitude, Longitude]
bsArrayOrientation1 = [-120 -10].';     % [Azimuth; Elevation] (degrees)
bsPosition2 = [52.251076,20.899524];    % [Latitude, Longitude]
bsArrayOrientation2 = [20 -10].';     % [Azimuth; Elevation] (degrees)

% Base station transmit power (dBm)
bsTxPower = 43;  % 5W = 37dBm 20W = 43dBm

% User Equipment (UE) 'mesuring points' positions 
uePositions = [
    52.252251, 20.900002; % UE-1
    52.251634, 20.902267; % UE-2
    52.255475, 20.900689; % UE-3
    52.252207, 20.905517; % UE-4
    52.256861, 20.899916; % UE-5
    52.255117, 20.903365; % UE-6
    52.255420, 20.902401; % UE-7
    52.255834, 20.899132; % UE-8
    52.253490, 20.902260; % UE-9
    52.253909, 20.901067; % UE-10
    52.254563, 20.901170; % UE-11
    52.254347, 20.904763; % UE-12
    52.254769, 20.901956; % UE-13
    52.255537, 20.899738; % UE-14
    52.251868, 20.906254; % UE-15
    52.251063, 20.905155; % UE-16
    52.251985, 20.904676; % UE-17
    52.249458, 20.904992; % UE-18
    52.249775, 20.901617; % UE-19
    52.249952, 20.904060; % UE-20
    52.252526, 20.902980; % UE-21
    52.251281, 20.903752; % UE-22
    52.252720, 20.901149; % UE-23
    52.253238, 20.901263; % UE-24
    52.253472, 20.903349; % UE-25
    52.253667, 20.904444; % UE-26
    52.254404, 20.903308; % UE-27
    52.251634, 20.900128; % UE-28
    52.252414, 20.901770; % UE-29
    52.254630, 20.899806; % UE-30
    52.252778, 20.902257; % UE-31
    52.253655, 20.905455; % UE-32
    52.250856, 20.901424; % UE-33
    52.252727, 20.903551; % UE-34
    52.250628, 20.904618; % UE-35
    52.251811, 20.905739; % UE-36
    52.249370, 20.901342; % UE-37
    52.254002, 20.905815; % UE-38
    52.254365, 20.904431; % UE-39
    52.252017, 20.898901; % UE-40
    52.253447, 20.904195; % UE-41
    52.255057, 20.904285; % UE-42
    52.252723, 20.904444; % UE-43
    52.257395, 20.899547; % UE-44
    52.254044, 20.900135; % UE-45
    52.250900, 20.900878; % UE-46
    52.252218, 20.900563; % UE-47
    52.253187, 20.899381; % UE-48
    52.250716, 20.903421; % UE-49
    52.251188, 20.902568  % UE-50
];
numUEs = size(uePositions, 1);
fprintf('Total number of UEs: %d\n', numUEs);

% Define UE heights and array orientations (same for all UEs)
ueHeight = 1.5; % meters
ueAntSize = [2 2];                      % Array size (rows, columns)
ueArrayOrientation = [180 45].';        % [Azimuth; Elevation] (degrees)

% Ray tracing configuration
reflectionsOrder = 4;                   % Number of reflections (0 for LOS only), max 10
diffractionOrder = 0;                   % max 1!! , default 0 <--

% OFDM/NR bandwidth configuration
SCS = 30;       % Subcarrier spacing (kHz) 15 
NRB = 273;       % Number of resource blocks 
% FR1: 10MHz=52 ,50MHz=133 ,100MHz=273
% FR2: n258 scs 60 nrb 132=100MHz, 264=200MHz

% Noise parameters
noiseFigure = 7;  % UE noise figure in dB
temperature = 290; % Temperature in Kelvin
kB = physconst('Boltzmann');  % Boltzmann constant
BW = NRB * 12 * SCS * 1000;   % Bandwidth in Hz (12 subcarriers per RB)
noiseFloor = 10*log10(kB*temperature*BW) + 30; % Noise floor in dBm (+30 for dBW to dBm)
noiseFloorWithNF = noiseFloor + noiseFigure;   % Adding noise figure
noiseFloor_mW = db2pow(noiseFloorWithNF);      % Convert to linear

% Antenna configuration selection
% Set antennaConfig to "64TRX_192AE", "32TRX_128AE", or "8TRX_64AE"
antennaConfig = "32TRX_128AE";  % Change this to switch antenna configurations

minRxPower = -120; % Minimum RX power for connection (dBm)

% Channel estimation error parameter (0 = perfect, higher values = more error)
channel_estimation_error = 0.15; % 15% error

% HARQ parameters
max_harq_retx = 3;  % Maximum number of HARQ retransmissions
target_bler = 0.1;  % Target block error rate (10%)

% For analysis with many UEs, define a subset of UEs to display in detail
detailed_analysis_ues = 1:min(10, numUEs); % Display detailed results for first 10 UEs

% Coordinated Beamforming parameters
coordBeamforming_enabled = true;    % Enable coordinated beamforming
coord_interference_threshold = -50; % dBm - threshold to trigger coordination

% Power Control and Cell Range Expansion parameters
powerControl_enabled = true;        % Enable dynamic power control
max_power_adjust = 3;               % Maximum power adjustment in dB

% Flag to enable/disable spatial multiplexing
spatial_multiplexing_enabled = true;  % Set to false for single-layer beamforming only

% Rank adaptation parameters
rank_adaptation_enabled = true;  % Dynamically select optimal number of layers
rank_sinr_threshold = 15;        % Minimum SINR (dB) required for additional layers


%% Viewer Initialization

% Create or clear the map viewer if it exists
if exist('viewer','var') && isvalid(viewer)
    clearMap(viewer);
else
    viewer = siteviewer("Basemap","openstreetmap","Buildings","wat_expanded.osm"); 
end

%% Define Antenna Array for Base Stations based on selected configuration

% Calculate wavelength
c = physconst('LightSpeed');
lambda = c / fc;

% Create antenna array based on selected configuration
if strcmpi(antennaConfig, "64TRX_192AE")
    % Create the 64 TRX 192 AE antenna array in 12×8×1×2 configuration
    % - 12 rows per panel
    % - 8 columns per panel
    % - 1 panel
    % - 2 polarizations per position
    bsArray = phased.NRRectangularPanelArray(...
        'Size', [12 8 1 1], ...
        'Spacing', 0.5 * lambda * [1 1 1 1]);
    
    % Set cross-polarized elements (-45° and +45° polarizations)
    bsArray.ElementSet = {phased.NRAntennaElement('PolarizationAngle', -45), ...
                         phased.NRAntennaElement('PolarizationAngle', 45)};
    
    % Calculate total number of antenna elements
    numAntElements = prod(bsArray.Size) * length(bsArray.ElementSet); % Should be 192
    
    % TRX configuration
    numTRX = 64;
    elementsPerTRX = numAntElements / numTRX; % Should be 3
    
    % Display configuration info
    fprintf('Using 64TRX 192AE Antenna Configuration (12×8×1×2)\n');
    
elseif strcmpi(antennaConfig, "32TRX_128AE")
    % Create the 32 TRX 128 AE antenna array in 8×8×1×2 configuration
    % - 8 rows per panel
    % - 8 columns per panel
    % - 1 panel
    % - 2 polarizations per position
    bsArray = phased.NRRectangularPanelArray(...
        'Size', [8 8 1 1], ...
        'Spacing', 0.5 * lambda * [1 1 1 1]);
    
    % Set cross-polarized elements (-45° and +45° polarizations)
    bsArray.ElementSet = {phased.NRAntennaElement('PolarizationAngle', -45), ...
                         phased.NRAntennaElement('PolarizationAngle', 45)};
    
    % Calculate total number of antenna elements
    numAntElements = prod(bsArray.Size) * length(bsArray.ElementSet); % Should be 128
    
    % TRX configuration
    numTRX = 32;
    elementsPerTRX = numAntElements / numTRX; % Should be 4
    
    % Display configuration info
    fprintf('Using 32TRX 128AE Antenna Configuration (8×8×1×2)\n');
    
elseif strcmpi(antennaConfig, "8TRX_64AE")
    % Create an 8 TRX 64 AE antenna array in 8×4×1×2 configuration
    % - 8 rows per panel
    % - 4 columns per panel
    % - 1 panel
    % - 2 polarizations per position
    bsArray = phased.NRRectangularPanelArray(...
        'Size', [8 4 1 1], ...
        'Spacing', 0.5 * lambda * [1 1 1 1]);
    
    % Set cross-polarized elements (-45° and +45° polarizations)
    bsArray.ElementSet = {phased.NRAntennaElement('PolarizationAngle', -45), ...
                         phased.NRAntennaElement('PolarizationAngle', 45)};
    
    % Calculate total number of antenna elements
    numAntElements = prod(bsArray.Size) * length(bsArray.ElementSet); % Should be 64
    
    % TRX configuration
    numTRX = 8;
    elementsPerTRX = numAntElements / numTRX; % Should be 8
    
    % Display configuration info
    fprintf('Using 8TRX 64AE Antenna Configuration (8×4×1×2)\n');
    
else
    error('Invalid antenna configuration. Choose either "64TRX_192AE", "32TRX_128AE", or "8TRX_64AE"');
end

% Print array information
fprintf('Base Station Antenna Array Properties:\n');
fprintf('- Total number of antenna elements: %d\n', numAntElements);
fprintf('- Number of TRX elements: %d\n', numTRX);
fprintf('- Each TRX connects to %d antenna elements\n', elementsPerTRX);

% Create TRX to antenna element mapping
TRX_mapping = createTRXtoAEMapping(numTRX, numAntElements);

% Create a simple UE antenna array (2x2 with single polarization)
ueArray = phased.NRRectangularPanelArray(...
    'Size', [ueAntSize(1:2) 1 1], ...
    'Spacing', 0.5 * lambda * [1 1 1 1]);

ueArray.ElementSet = {phased.IsotropicAntennaElement};
figure;
pattern(ueArray, fc, 'ShowArray', true);
title('Radiation Pattern of UE Antenna Array');
figure;
pattern(bsArray, fc, 'ShowArray', true);
title(['Radiation Pattern of Base Station Antenna Array - ', antennaConfig]);


% MIMO Spatial Multiplexing Configuration
% Set max_spatial_layers based on the antenna configuration
if strcmpi(antennaConfig, "64TRX_192AE")
    max_spatial_layers = 4;  % Up to 4 layers for 64TRX
elseif strcmpi(antennaConfig, "32TRX_128AE")
    max_spatial_layers = 3;  % Up to 3 layers for 32TRX
else % 8TRX_64AE
    max_spatial_layers = 2;  % Up to 2 layers for 8TRX
end
%% Define Base Station and UE Sites

% Base Station Sites - All operating on same frequency (Band 77)
bsSite1 = txsite("Name","Base station 1", ...
    "Latitude",bsPosition1(1), "Longitude",bsPosition1(2),...
    "AntennaAngle",bsArrayOrientation1(1:2),...
    "AntennaHeight",4,...  % Height in meters
    "TransmitterFrequency",fc, ...
    "TransmitterPower", db2pow(bsTxPower)/1000); % Convert dBm to watts

bsSite2 = txsite("Name","Base Station 2", ...
    "Latitude",bsPosition2(1), "Longitude",bsPosition2(2), ...  
    "AntennaAngle",bsArrayOrientation2(1:2), ...       
    "AntennaHeight",4, ...                            
    "TransmitterFrequency",fc, ...
    "TransmitterPower", db2pow(bsTxPower)/1000);
%FIRST CONFIG
% bsSite3 = txsite("Name","Base Station 3", ...
%     "Latitude",bsPosition3(1), "Longitude",bsPosition3(2), ...  
%     "AntennaAngle",bsArrayOrientation3(1:2), ...       
%     "AntennaHeight",4, ...                            
%     "TransmitterFrequency",fc, ... 
%     "TransmitterPower", db2pow(bsTxPower)/1000);



% Store all BS sites in an array for easier access
% Change this depending on 1st or 2nd config
bsSites = {bsSite1, bsSite2}; % , bsSite3
numBS = length(bsSites);

%Create UE Sites
ueSites = cell(numUEs, 1);
for i = 1:numUEs
    ueSites{i} = rxsite("Name", "UE-" + i, ...
        "Latitude", uePositions(i, 1), "Longitude", uePositions(i, 2), ...
        "AntennaHeight", ueHeight, ...
        "AntennaAngle", ueArrayOrientation(1:2));

    % Show UE site with default appearance
    show(ueSites{i});
end

% Display the BS sites on the map
for bs = 1:numBS
    show(bsSites{bs});
end


%% Create a propagation model with ray tracing (SBR method)
pm = propagationModel("raytracing", "Method", "sbr", ...
    "MaxNumReflections", reflectionsOrder, "MaxNumDiffractions", diffractionOrder, ...
    "TerrainMaterial","vegetation",...
    "BuildingsMaterial", "concrete");

pm = pm + propagationModel("rain") + propagationModel("gas");

%% Initialize data structures for storing ray tracing and signal results
% Initialize 3D arrays to store results: [BS, UE, metric]
pathLoss = zeros(numBS, numUEs);
rxPower = zeros(numBS, numUEs);
sinrValues = zeros(numBS, numUEs);
sinrValuesCapped = zeros(numBS, numUEs); 
throughput = zeros(numBS, numUEs);
validConnections = false(numBS, numUEs);  % to track valid connections

%Initialize structures for AMC and HARQ analysis
cqiValues = zeros(numBS, numUEs);
mcsValues = zeros(numBS, numUEs);
spectralEfficiency = zeros(numBS, numUEs);
effectiveThroughput = zeros(numBS, numUEs); % After HARQ
blerEstimates = zeros(numBS, numUEs);

% Create arrays to store all ray tracing results
allRays = cell(numBS, numUEs);            % [BS, UE]
allBeamformingGains = zeros(numBS, numUEs);
allPathToAs = cell(numBS, numUEs);
allAvgPathGains = cell(numBS, numUEs);
allPathAoDs = cell(numBS, numUEs);
allPathAoAs = cell(numBS, numUEs);
allIsLOS = zeros(numBS, numUEs);
allChannels = cell(numBS, numUEs);

% Initialize the rtChannels array to prevent the error
rtChannels = cell(numBS, numUEs);

% Obtain OFDM information based on NRB and SCS
ofdmInfo = nrOFDMInfo(NRB, SCS);

% Initialize structures for MIMO analysis
allNumLayers = cell(numBS, numUEs);  % Store number of spatial layers used
layerDistribution = zeros(1, max_spatial_layers); % Count UEs using each rank

%% Ray Tracing for all BS-UE pairs
% Loop through all BS-UE combinations
for bs = 1:numBS
    fprintf('\nProcessing Base Station %d:\n', bs);
    
    for ue = 1:numUEs
        fprintf('  Processing UE-%d...\n', ue);
        
        % Perform ray tracing from current BS to current UE
        rays = raytrace(bsSites{bs}, ueSites{ue}, pm, "Type", "pathloss");
        
        % Store the rays
        if ~isempty(rays) && ~isempty(rays{1})
            allRays{bs, ue} = rays{1};
            
            % Only plot the ray-tracing result for a subset of UEs to avoid clutter
            % Plot just the first few UEs for each BS for visibility
            if ismember(ue, detailed_analysis_ues)
                plot(rays{1});
            end
            
            % Process ray tracing results
            % Normalize propagation delays to start at 0 sec
            allPathToAs{bs, ue} = [allRays{bs, ue}.PropagationDelay] - min([allRays{bs, ue}.PropagationDelay]);
            
            % Extract average path gains (negated path loss values)
            allAvgPathGains{bs, ue} = -[allRays{bs, ue}.PathLoss];
            
            % Store the total path loss (minimum path loss from all rays)
            pathLoss(bs, ue) = min([allRays{bs, ue}.PathLoss]);
            
            % Extract angles of departure and arrival
            allPathAoDs{bs, ue} = [allRays{bs, ue}.AngleOfDeparture];  % [Azimuth; Elevation]
            allPathAoAs{bs, ue} = [allRays{bs, ue}.AngleOfArrival];
            
            % Check if there is a Line-Of-Sight component
            allIsLOS(bs, ue) = any([allRays{bs, ue}.LineOfSight]);
            
            % Calculate received power (dBm)
            rxPower(bs, ue) = bsTxPower - pathLoss(bs, ue);
            
            % Check if connection is valid (power above threshold)
            validConnections(bs, ue) = (rxPower(bs, ue) > minRxPower);
            
            if validConnections(bs, ue)
                fprintf('    Valid connection! RSRP: %.2f dBm\n', rxPower(bs, ue));
                
                % Configure CDL Channel for this BS-UE pair
                channel = nrCDLChannel;
                channel.DelayProfile = 'Custom';
                channel.PathDelays = allPathToAs{bs, ue};
                channel.AveragePathGains = allAvgPathGains{bs, ue};
                
                % Set angles (convert elevation to zenith angle: zenith = 90 - elevation)
                channel.AnglesAoD = allPathAoDs{bs, ue}(1,:);       
                channel.AnglesZoD = 90 - allPathAoDs{bs, ue}(2,:);    
                channel.AnglesAoA = allPathAoAs{bs, ue}(1,:);       
                channel.AnglesZoA = 90 - allPathAoAs{bs, ue}(2,:);    
                
                channel.HasLOSCluster = allIsLOS(bs, ue);
                channel.CarrierFrequency = fc;
                channel.NormalizeChannelOutputs = false;
                channel.NormalizePathGains = false;
                
                % Assign antenna arrays - using the configured array for base stations
                channel.ReceiveAntennaArray = ueArray;
                channel.TransmitAntennaArray = bsArray;
                
                % Set array orientations based on which base station we're processing
                switch bs
                    case 1
                        channel.TransmitArrayOrientation = [bsArrayOrientation1(1); -bsArrayOrientation1(2); 0];
                    case 2
                        channel.TransmitArrayOrientation = [bsArrayOrientation2(1); -bsArrayOrientation2(2); 0];
                    case 3
                        channel.TransmitArrayOrientation = [bsArrayOrientation3(1); -bsArrayOrientation3(2); 0];
                end
                
                % Set UE array orientation
                channel.ReceiveArrayOrientation = [ueArrayOrientation(1); -ueArrayOrientation(2); 0];
                
                % Set sample rate from OFDM parameters
                channel.SampleRate = ofdmInfo.SampleRate;
                channel.ChannelFiltering = false;
                
                % Store the configured channel
                allChannels{bs, ue} = channel;
                rtChannels{bs, ue} = channel; 
                
                % Generate channel samples for this BS-UE pair
                [pathGains, sampleTimes] = channel();
                
                % Rearrange path gains dimensions
                pg = permute(pathGains, [2 1 3 4]); 
                
                % If LOS is present, sum the first two paths (assumed to be the LOS ray)
                if allIsLOS(bs, ue)
                    if size(pg, 1) >= 2  % Make sure there are at least 2 paths
                        pg = [sum(pg(1:2,:,:,:)); pg(3:end,:,:,:)];
                    end
                end
                pg = abs(pg).^2; % Compute power gains
                
                % Channel Estimation
                pathFilters = getPathFilters(channel);
                nSlot = 0;  % Initial slot number
                
                % Estimate timing offset using perfect timing estimation
                [offset, ~] = nrPerfectTimingEstimate(pathGains, pathFilters);
                
                % Obtain channel estimate
                hest = nrPerfectChannelEstimate(pathGains, pathFilters, NRB, SCS, nSlot, offset, sampleTimes);
                
                % Calculate SINR approximation for rank adaptation before beamforming
                % First calculate SINR
                signal_power = rxPower(bs, ue);
                noise_power = noiseFloorWithNF;
                
                % Simple interference estimation from other cells
                 interference_power = 0;
                interference_suppression_factor = 0.7; % 0.7 means 30% interference reduction
                for other_bs = 1:numBS
                    if other_bs ~= bs && validConnections(other_bs, ue)
                                interference_power = interference_power + db2pow(rxPower(other_bs, ue) - 3) * interference_suppression_factor;
                    end
                end
                
                if interference_power > 0
                    interference_power_dB = pow2db(interference_power);
                else
                    interference_power_dB = -120;
                end
                
                % Estimated SINR
                est_sinr_dB = signal_power - max(interference_power_dB, noise_power);
                
                if rank_adaptation_enabled && spatial_multiplexing_enabled
                    % Determine optimal number of layers based on SINR and channel conditions
                    if est_sinr_dB >= rank_sinr_threshold
                        % Get channel matrix through SVD for rank analysis
                        [~, ~, R, P] = size(hest);
                        scNo = 1;
                        hEstSubset = hest(scNo:scNo + 11, :, :, :); % Use first RB
                        H = permute(mean(reshape(hEstSubset, [], R, P)), [2 3 1]);
                        
                        % Perform SVD to get singular values for rank determination
                        [~, S, ~] = svd(H);
                        S = diag(S);
                        
                        % Determine optimal rank based on singular value distribution
                        num_spatial_layers = 1; % Default is rank 1
                        for layer = 2:min(max_spatial_layers, length(S))
                            % Add layer if singular value ratio is good enough 
                            % (within 10dB of strongest layer)
                            singular_value_ratio_dB = 20*log10(S(layer)/S(1));
                            if singular_value_ratio_dB > -10
                                num_spatial_layers = layer;
                            else
                                break; % Stop if singular value too small
                            end
                        end
                        
                        % Cap based on SINR
                        % Higher SINR allows more layers
                        if est_sinr_dB >= 25 
                            % Allow up to max layers at excellent SINR
                            num_spatial_layers = min(num_spatial_layers, max_spatial_layers);
                        elseif est_sinr_dB >= 18
                            % Good SINR - up to 3 layers 
                            num_spatial_layers = min(num_spatial_layers, 3);
                        elseif est_sinr_dB >= rank_sinr_threshold
                            % Decent SINR - up to 2 layers
                            num_spatial_layers = min(num_spatial_layers, 2);
                        else
                            % Low SINR - single layer only
                            num_spatial_layers = 1;
                        end
                    else
                        % Use single layer for low SINR conditions
                        num_spatial_layers = 1;
                    end
                else
                    % If rank adaptation is disabled, use configured value
                    if spatial_multiplexing_enabled
                        num_spatial_layers = max_spatial_layers;
                    else
                        num_spatial_layers = 1;
                    end
                end

                % Store the number of layers used for this BS-UE pair
                allNumLayers{bs, ue} = num_spatial_layers;

                % Log rank adaptation results for detailed analysis UEs
                if ismember(ue, detailed_analysis_ues)
                    fprintf('    Estimated SINR: %.2f dB\n', est_sinr_dB);
                    
                    % Create adaptation text based on condition
                    if rank_adaptation_enabled
                        adaptation_text = '(rank-adapted)';
                    else
                        adaptation_text = '(fixed)';
                    end
                    
                    fprintf('    Using %d spatial layer(s) %s\n', num_spatial_layers, adaptation_text);
                end
                
                % Calculate beamforming weights for the TRX system with estimation error
                scOffset = 0;
                noRBs = 1;
                [wbs_full, wue, singular_values] = getBeamformingWeightsWithError(hest, num_spatial_layers, scOffset, noRBs, channel_estimation_error);
                
                % Apply coordinated beamforming if enabled
                if exist('coordBeamforming_enabled', 'var') && coordBeamforming_enabled
                    % Only coordinate if this UE is experiencing significant interference
                    needs_coordination = false;
                    
                    % Check if there's significant interference from other cells
                    for other_bs = 1:numBS
                        if other_bs ~= bs && validConnections(other_bs, ue)
                            if rxPower(other_bs, ue) > coord_interference_threshold
                                needs_coordination = true;
                                break;
                            end
                        end
                    end
                    
                    if needs_coordination
                        % Apply coordinated beamforming
                        wbs_full = calculateCoordinatedBeamWeights(hest, bs, ue, allChannels, validConnections, rtChannels, rxPower);
                        fprintf('    Applying coordinated beamforming for UE-%d at BS-%d\n', ue, bs);
                    end
                end
                                
                % Convert full antenna element weights to TRX control signals
                wbs_trx = convertBeamWeightsToTRX(wbs_full, TRX_mapping, numTRX, numAntElements);
                
                % Apply TRX-level beamforming (convert back to full antenna elements)
                wbs = applyTRXBeamWeights(wbs_trx, TRX_mapping, numAntElements);
                
                % Store beamforming weights
                bsSites{bs}.Antenna = clone(channel.TransmitAntennaArray);
                bsSites{bs}.Antenna.Taper = wbs;
                
                % Calculate beamforming gain with a more realistic model
                % Theoretical array gain: 10*log10(numAntElements)
                theoreticalGain = 10*log10(numAntElements);
                directivityBonus = 4; % Reduced from 4dB to 2dB? 2?
                
                % Apply implementation losses
                implementationLoss = 3; % 3dB implementation loss due to hardware imperfections 1.5?
                
                % Calculate realistic gain and cap at 26dB
                realGain = min(26, theoreticalGain + directivityBonus - implementationLoss);
                allBeamformingGains(bs, ue) = realGain;
                
                % % Add logging for debugging
                % fprintf('    Beamforming gain: %.2f dB (theoretical: %.2f dB, capped at 26 dB)\n', ...
                %         allBeamformingGains(bs, ue), theoreticalGain + directivityBonus);
                
                % Update received power with beamforming gain
                rxPower(bs, ue) = bsTxPower - pathLoss(bs, ue) + allBeamformingGains(bs, ue);
                
                % Only show antenna pattern for detailed analysis UEs
                if ismember(ue, detailed_analysis_ues)
                    pattern(bsSites{bs}, fc, "Size", 8);
                end
                
            else
                fprintf('    Connection below threshold (%.2f dBm)\n', rxPower(bs, ue));
            end
        else
            fprintf('    No paths found from BS-%d to UE-%d\n', bs, ue);
            validConnections(bs, ue) = false;
        end
    end
end

% Set UE antennas with beamforming for valid UEs
for ue = 1:numUEs
    % Find last valid BS for this UE
    lastValidBS = find(validConnections(:, ue), 1, 'last');
    if ~isempty(lastValidBS) && ismember(ue, detailed_analysis_ues)
        ueSites{ue}.Antenna = clone(ueArray);
        if size(wue, 1) > 1
            % For multi-layer case, use just the first layer weights
            ueSites{ue}.Antenna.Taper = wue(1,:);
            fprintf('    Note: Using only first layer weights for UE antenna pattern\n');
        else
            ueSites{ue}.Antenna.Taper = wue;
        end
        
        % Only show pattern for detailed analysis UEs
        pattern(ueSites{ue}, fc, "Size", 4);
    end
end

%% Apply Power Control
if exist('powerControl_enabled', 'var') && powerControl_enabled
    fprintf('\n*** Applying Power Control ***\n');
    
    % Adjust BS transmit power based on path loss statistics
    for bs = 1:numBS
        valid_ues = find(validConnections(bs, :));
        if ~isempty(valid_ues)
            % Calculate average path loss to valid UEs
            avg_pathLoss = mean(pathLoss(bs, valid_ues));
            
            % Calculate standard deviation
            std_pathLoss = std(pathLoss(bs, valid_ues));
            
            % Adjust power based on path loss statistics
            power_adjustment = 0;
            
            if avg_pathLoss > 105  % High path loss
                power_adjustment = max_power_adjust;  % Increase power
            elseif avg_pathLoss < 95  % Low path loss
                power_adjustment = -max_power_adjust;  % Decrease power
            elseif std_pathLoss > 15  % High variance in path loss
                power_adjustment = max_power_adjust/2;  % Moderate increase
            end
            
            % Apply power adjustment
            new_power_dBm = bsTxPower + power_adjustment;
            fprintf('  BS-%d: Avg PathLoss: %.2f dB, Power adjustment: %.1f dB, New TxPower: %.1f dBm\n', ...
                    bs, avg_pathLoss, power_adjustment, new_power_dBm);
            
            % Update BS transmit power
            bsSites{bs}.TransmitterPower = db2pow(new_power_dBm)/1000;
            
            % Update all rxPower values for this BS
            for ue = 1:numUEs
                if validConnections(bs, ue)
                    rxPower(bs, ue) = rxPower(bs, ue) + power_adjustment;
                end
            end
        end
    end
end

%% Calculate SINR and Throughput for each valid BS-UE pair with improved SINR calculation

% Calculate thermal noise power
kB = physconst('Boltzmann');  % Boltzmann constant
BW = NRB * 12 * SCS * 1000;   % Bandwidth in Hz (12 subcarriers per RB)
noiseFloor = 10*log10(kB*temperature*BW) + 30; % Noise floor in dBm (+30 for dBW to dBm)
noiseFloorWithNF = noiseFloor + noiseFigure;   % Adding noise figure

% Convert background interference to mW
background_interference_dBm = -95;
background_interference_mW = db2pow(background_interference_dBm);

fprintf('\n*** Interference Analysis ***\n');
fprintf('Noise Floor: %.2f dBm\n', noiseFloorWithNF);
fprintf('Background Interference: %.2f dBm\n', background_interference_dBm);

% Display received signal powers only for detailed analysis UEs
fprintf('\nReceived Signal Powers at Selected UEs:\n');
for ue = detailed_analysis_ues
    fprintf('UE-%d:\n', ue);
    for bs = 1:numBS
        if validConnections(bs, ue)
            fprintf('  From BS-%d: %.2f dBm (with %.2f dB beamforming gain)\n', ...
                bs, rxPower(bs, ue), allBeamformingGains(bs, ue));
        end
    end
end

% Create tracking matrix for coordinated beamforming
coordBeamforming_applied = false(numBS, numUEs);

% Calculate SINR for each valid BS-UE pair with improved interference modeling
% First, process all UEs for calculation
for ue = 1:numUEs
    % Skip detailed logging for most UEs to avoid console clutter
    if ismember(ue, detailed_analysis_ues)
        fprintf('\nSINR Analysis for UE-%d:\n', ue);
    end
    
    % Check if there are any valid connections for this UE
    if ~any(validConnections(:, ue))
        if ismember(ue, detailed_analysis_ues)
            fprintf('  No valid connections for UE-%d\n', ue);
        end
        continue;
    end
    
    for servingBS = 1:numBS
        % Skip if this is not a valid connection
        if ~validConnections(servingBS, ue)
            continue;
        end
        
        % Signal power from serving cell (mW)
        signal_power_mW = db2pow(rxPower(servingBS, ue));
        
        % Enhanced interference calculation with dynamic ICIC
        [total_interference_mW, effective_icic_factors] = calculateEnhancedInterference(rxPower, servingBS, ue, validConnections, coordBeamforming_applied, sinrValuesCapped);        
        
        % Noise power (mW)
        noise_power_mW = db2pow(noiseFloorWithNF);
        
        % Calculate SINR
        sinr_mW = signal_power_mW / (total_interference_mW + noise_power_mW);
        sinrValues(servingBS, ue) = pow2db(sinr_mW);
        
        % Apply SINR capping to realistic values
        sinrValuesCapped(servingBS, ue) = capSINR(sinrValues(servingBS, ue));
        
        % Get the SINR value for this connection
        sinr_dB = sinrValuesCapped(servingBS, ue);
        
        % Determine modulation efficiency based on SINR regions
        % These thresholds align with 3GPP CQI to modulation mapping
        if sinr_dB >= 22.7 % 256-QAM region (CQI 15-16)
            % At excellent SINR, we can get closer to Shannon limit
            spectral_efficiency_factor = 0.85;
        elseif sinr_dB >= 14.1 % 64-QAM region (CQI 10-14)
            % Medium-high SINR
            spectral_efficiency_factor = 0.75;
        elseif sinr_dB >= 5.9 % 16-QAM region (CQI 7-9)
            % Medium SINR
            spectral_efficiency_factor = 0.65;
        elseif sinr_dB >= -2.3 % QPSK region (CQI 3-6)
            % Low SINR
            spectral_efficiency_factor = 0.55;
        else % Lowest QPSK region (CQI 1-2)
            % Very low SINR, very far from Shannon limit
            spectral_efficiency_factor = 0.45;
        end
        
        % Apply implementation loss scaling (larger arrays have greater implementation losses)
        implementation_loss_dB = 2 + log10(numAntElements) * 0.5; 
        
        % Apply implementation loss to SINR
        effective_sinr_dB = sinr_dB - implementation_loss_dB;
        effective_sinr_dB = max(effective_sinr_dB, -10); % Ensure SINR doesn't go below minimum
        
        % Calculate spectral efficiency using the variable efficiency factor
        shannon_limit = log2(1 + db2pow(effective_sinr_dB));
        spectral_efficiency_value = spectral_efficiency_factor * shannon_limit;
        
        % Account for system overhead (control signals, pilots, guard bands)
        % High modulation orders have less relative overhead
        if sinr_dB >= 14.1 % 64-QAM and 256-QAM
            overhead_factor = 0.88; % 12% overhead
        else
            overhead_factor = 0.85; % 15% overhead for lower modulation schemes
        end
        
        spectral_efficiency_value = spectral_efficiency_value * overhead_factor;
        
        % Get CQI, modulation order and theoretical spectral efficiency for this SINR
        [cqi, ~, modOrder, ~, theoretical_se] = getLinkAdaptationParameters(sinr_dB);
        
        % Apply a final upper bound based on the modulation order
        % Maximum spectral efficiency for each modulation:
        % QPSK: 2 bps/Hz, 16-QAM: 4 bps/Hz, 64-QAM: 6 bps/Hz, 256-QAM: 8 bps/Hz
        max_se_for_modulation = modOrder;
        
        % Ensure we don't exceed the theoretical maximum for the modulation
        spectral_efficiency_value = min(spectral_efficiency_value, max_se_for_modulation * 0.95);
        
        % Store the final realistic spectral efficiency
        spectralEfficiency(servingBS, ue) = spectral_efficiency_value;
        
        % Calculate raw throughput with spatial multiplexing
        if spatial_multiplexing_enabled && allNumLayers{servingBS, ue} > 1
            num_spatial_layers = allNumLayers{servingBS, ue};
            
            % Base throughput for one layer
            base_throughput = spectral_efficiency_value * BW / 1e6;
            
            % Calculate layer efficiencies
            layer_efficiency = zeros(1, num_spatial_layers);
            layer_efficiency(1) = 1.0;
            
            for l = 2:num_spatial_layers
                layer_efficiency(l) = max(0.3, 1.0 - (l-1)*0.15);
            end
            
            % Calculate total throughput with proper multiplexing gain
            spatial_gain = sum(layer_efficiency);
            max_throughput = base_throughput * spatial_gain;
            
            % fprintf('    DEBUG: MIMO throughput calculation:\n');
            % fprintf('    DEBUG: Base throughput (1 layer): %.2f Mbps\n', base_throughput);
            % fprintf('    DEBUG: Layer efficiencies: [%s]\n', sprintf('%.2f ', layer_efficiency));
            % fprintf('    DEBUG: Spatial multiplexing gain: %.2fx\n', spatial_gain);
            % fprintf('    DEBUG: Final throughput: %.2f Mbps\n', max_throughput);
            
            if ismember(ue, detailed_analysis_ues)
                fprintf('    Using %d spatial layers with efficiencies: [%s]\n',...
                        num_spatial_layers, sprintf('%.2f ', layer_efficiency));
                fprintf('    Spatial multiplexing gain: %.2f x\n', spatial_gain);
            end
        else
            % Single-layer throughput
            max_throughput = spectral_efficiency_value * BW / 1e6;
        end
        % Implement AMC - select CQI and MCS based on SINR
        [cqi, mcs, modOrder, codeRate, spectral_eff_cqi] = getLinkAdaptationParameters(sinrValuesCapped(servingBS, ue));
        cqiValues(servingBS, ue) = cqi;
        mcsValues(servingBS, ue) = mcs;
        
        % Calculate BLER estimate for the selected MCS
        [~, bler] = linkAdaptation(sinrValuesCapped(servingBS, ue), target_bler);
        blerEstimates(servingBS, ue) = bler;
        
        % Calculate effective throughput with improved HARQ model (Chase combining)
        effective_tput = calculateImprovedHARQThroughput(max_throughput, sinrValuesCapped(servingBS, ue), bler, max_harq_retx);
        effectiveThroughput(servingBS, ue) = effective_tput;
        
        % Only print detailed SINR analysis for selected UEs
        if ismember(ue, detailed_analysis_ues)
            % Calculate SIR (Signal to Interference Ratio)
            % Use max to avoid division by zero
            if total_interference_mW > 1e-10
                sir_mW = signal_power_mW / total_interference_mW;
            else
                sir_mW = signal_power_mW / 1e-10;
            end
        
            % Calculate SNR (Signal to Noise Ratio)
            snr_mW = signal_power_mW / noise_power_mW;
            
            % Print detailed ICIC information
            fprintf('  If BS-%d is serving cell:\n', servingBS);
            fprintf('    Signal Power: %.2f dBm\n', rxPower(servingBS, ue));
            fprintf('    Total Interference Power: %.2f dBm\n', pow2db(total_interference_mW));
            
            % Print ICIC factors applied
            fprintf('    ICIC factors applied: ');
            icic_applied = false;
            for bs = 1:numBS
                if bs ~= servingBS && effective_icic_factors(bs) > 0
                    fprintf('BS-%d: %.2f  ', bs, effective_icic_factors(bs));
                    icic_applied = true;
                end
            end
            if ~icic_applied
                fprintf('None');
            end
            fprintf('\n');
            
            fprintf('    Signal-to-Interference Ratio (SIR): %.2f dB\n', pow2db(sir_mW));
            fprintf('    Signal-to-Noise Ratio (SNR): %.2f dB\n', pow2db(snr_mW));
            fprintf('    Uncapped SINR: %.2f dB\n', sinrValues(servingBS, ue));
            fprintf('    Capped SINR: %.2f dB\n', sinrValuesCapped(servingBS, ue));
            fprintf('    CQI: %d (Modulation Order: %d-QAM, Code Rate: %.4f)\n', cqi, 2^modOrder, codeRate);
            fprintf('    Estimated BLER: %.4f\n', bler);
            fprintf('    Spectral Efficiency: %.2f bps/Hz\n', spectral_efficiency_value);
            fprintf('    Raw Throughput: %.2f Mbps\n', max_throughput);
            fprintf('    Effective Throughput with HARQ: %.2f Mbps\n', effective_tput);
        end
    end
end

%% Find the best serving cell for each UE based on effective throughput

% Determine serving cells based on throughput
fprintf('\n*** Best Serving Cell Analysis ***\n');
best_serving_bs = zeros(numUEs, 1);

for ue = 1:numUEs
    best_val = -inf;
    best_idx = 0;
    for bs = 1:numBS
        if validConnections(bs, ue) && effectiveThroughput(bs, ue) > best_val
            best_val = effectiveThroughput(bs, ue);
            best_idx = bs;
        end
    end
    if best_idx > 0
        best_serving_bs(ue) = best_idx;
    end
end

% Count connections per UE and BS
bs_ue_counts = histcounts(best_serving_bs(best_serving_bs > 0), 1:numBS+1);
ue_connection_count = sum(validConnections, 1);
bs_connection_count = sum(validConnections, 2);

fprintf('\nConnection statistics:\n');
fprintf('- Number of UEs with valid connections: %d out of %d (%.1f%%)\n', ...
    sum(ue_connection_count > 0), numUEs, sum(ue_connection_count > 0)/numUEs*100);
fprintf('- Average number of BSs connected per UE: %.2f\n', mean(ue_connection_count(ue_connection_count > 0)));
fprintf('- Average number of UEs connected per BS: %.2f\n', mean(bs_connection_count));
% 3-BS print statement
if numBS == 2
    fprintf('- UE distribution by serving BS: BS1=%d, BS2=%d\n', bs_ue_counts(1), bs_ue_counts(2));
elseif numBS == 3
    fprintf('- UE distribution by serving BS: BS1=%d, BS2=%d, BS3=%d\n', bs_ue_counts(1), bs_ue_counts(2), bs_ue_counts(3));
else
    % For any other number of BS (1 or more than 3)
    str = '- UE distribution by serving BS: ';
    for bs = 1:numBS
        if bs < numBS
            str = str + sprintf('BS%d=%d, ', bs, bs_ue_counts(bs));
        else
            str = str + sprintf('BS%d=%d\n', bs, bs_ue_counts(bs));
        end
    end
    fprintf(str);
end

%% Visualizations - Grouped into 3 logical categories
% Define a peak theoretical throughput
peak_theoretical_rate = 1000; % Mbps

% Create color map for UEs
ueColorMap = jet(numUEs);

% Pre-calculate data needed for multiple plots
% Number of UEs served by each BS (only serving connections)
ues_per_bs = bs_ue_counts;

% Average SINR per BS (only serving connections)
avg_sinr_per_bs = zeros(numBS, 1);
for bs = 1:numBS
    valid_serving_ues = find(best_serving_bs == bs);
    if ~isempty(valid_serving_ues)
        avg_sinr_per_bs(bs) = mean(sinrValuesCapped(bs, valid_serving_ues));
    end
end

% Throughput metrics (only serving connections)
total_throughput_per_bs = zeros(numBS, 1);
avg_throughput_per_ue_per_bs = zeros(numBS, 1);
for bs = 1:numBS
    valid_serving_ues = find(best_serving_bs == bs);
    if ~isempty(valid_serving_ues)
        total_throughput_per_bs(bs) = sum(effectiveThroughput(bs, valid_serving_ues));
        avg_throughput_per_ue_per_bs(bs) = mean(effectiveThroughput(bs, valid_serving_ues));
    end
end

% Spectral efficiency per BS (only serving connections)
network_spectral_efficiency = zeros(numBS, 1);
for bs = 1:numBS
    bs_efficiency = 0;
    valid_serving_ues = find(best_serving_bs == bs);
    for ue = valid_serving_ues'
        bs_efficiency = bs_efficiency + spectralEfficiency(bs, ue);
    end
    network_spectral_efficiency(bs) = bs_efficiency;
end

% All throughputs and CQIs for distribution analysis (only serving connections)
all_throughputs = [];
all_cqis = [];
valid_sinrs = [];
valid_tputs = [];
for ue = 1:numUEs
    if best_serving_bs(ue) > 0
        bs = best_serving_bs(ue);
        all_throughputs = [all_throughputs; effectiveThroughput(bs, ue)];
        all_cqis = [all_cqis; cqiValues(bs, ue)];
        valid_sinrs = [valid_sinrs, sinrValuesCapped(bs, ue)];
        valid_tputs = [valid_tputs, effectiveThroughput(bs, ue)];
    end
end

% Modulation types distribution (only serving connections)
mod_types = {'QPSK', '16-QAM', '64-QAM', '256-QAM'};
mod_counts = zeros(1, 4);
for ue = 1:numUEs
    if best_serving_bs(ue) > 0
        bs = best_serving_bs(ue);
        cqi = cqiValues(bs, ue);
        if cqi <= 6
            mod_counts(1) = mod_counts(1) + 1; % QPSK
        elseif cqi <= 9
            mod_counts(2) = mod_counts(2) + 1; % 16-QAM
        elseif cqi <= 14
            mod_counts(3) = mod_counts(3) + 1; % 64-QAM
        else
            mod_counts(4) = mod_counts(4) + 1; % 256-QAM
        end
    end
end
non_zero_idx = mod_counts > 0;

% Calculate network capacity (only serving connections)
total_capacity_per_bs = zeros(numBS, 1);
for bs = 1:numBS
    valid_serving_ues = find(best_serving_bs == bs);
    total_capacity_per_bs(bs) = sum(effectiveThroughput(bs, valid_serving_ues));
end
total_network_capacity = sum(total_capacity_per_bs);

%% Group 1: Base Station Performance Dashboard
figure('Name', 'BS Performance Dashboard - Serving Connections Only', 'Position', [50, 50, 1000, 800]);

% Number of UEs served by each BS
subplot(2, 3, 1);
bar(ues_per_bs);
grid on;
xlabel('Base Station ID');
ylabel('Number of UEs Served');
title('UEs Connected to Each BS (Serving Only)');
xticks(1:numBS);

% Average SINR per BS
subplot(2, 3, 2);
bar(avg_sinr_per_bs);
grid on;
xlabel('Base Station ID');
ylabel('Average SINR (dB)');
title('Average SINR per BS');
xticks(1:numBS);

% Total throughput per BS
subplot(2, 3, 3);
bar(total_throughput_per_bs);
grid on;
xlabel('Base Station ID');
ylabel('Total Throughput (Mbps)');
title('Total Throughput per BS');
xticks(1:numBS);

% Average throughput per UE for each BS
subplot(2, 3, 4);
bar(avg_throughput_per_ue_per_bs);
grid on;
xlabel('Base Station ID');
ylabel('Avg Throughput per UE (Mbps)');
title('Average UE Throughput per BS');
xticks(1:numBS);

% Spectral efficiency per BS
subplot(2, 3, 5);
bar(network_spectral_efficiency);
grid on;
xlabel('Base Station ID');
ylabel('Spectral Efficiency (bps/Hz)');
title('Total Network Spectral Efficiency per BS');
xticks(1:numBS);

% UE distribution by serving BS
subplot(2, 3, 6);
if sum(bs_ue_counts) > 0
    pie(bs_ue_counts);
    legend(arrayfun(@(x) ['BS-' num2str(x)], 1:numBS, 'UniformOutput', false), 'Location', 'eastoutside');
else
    text(0.5, 0.5, 'No valid connections', 'HorizontalAlignment', 'center');
    axis off;
end
hold off;

% Add overall title
sgtitle(['Link Adaptation & AMC - ' char(antennaConfig) ' Configuration (Serving Only)'], 'FontSize', 16, 'FontWeight', 'bold');

%% Helper Functions

function sinr_capped = capSINR(sinr_dB)
    % Cap SINR between practical limits
    sinr_min = -10; % dB
    sinr_max = 30;  % dB
    sinr_capped = max(min(sinr_dB, sinr_max), sinr_min);
end

function [interference_power_mW, effective_icic_factors] = calculateEnhancedInterference(rxPower, servingBS, ue, validConnections, coordBeamforming_applied, sinrValuesCapped)
    % Optimized interference calculation for small private networks with variable BS count
    % Inputs remain the same as before
    
    % Initialize interference and tracking of ICIC factors
    interference_power_mW = 0;
    numBS = size(rxPower, 1);
    effective_icic_factors = zeros(1, numBS);
    
    % Get serving power for reference
    serving_power = rxPower(servingBS, ue);
    
    % Create an interference priority list to identify strongest interferers
    interference_list = [];
    for bs = 1:numBS
        if bs ~= servingBS && validConnections(bs, ue)
            % Store BS index and power
            interference_list = [interference_list; bs, rxPower(bs, ue)];
        end
    end
    
    % Sort by descending interference power
    if ~isempty(interference_list)
        interference_list = sortrows(interference_list, -2);
    end
    
    % Process interferers in order of significance
    for i = 1:size(interference_list, 1)
        interferingBS = interference_list(i, 1);
        interferer_power = interference_list(i, 2);
        
        % Convert to mW for calculations
        interferer_power_mW = db2pow(interferer_power);
        
        % Calculate interference severity (dB)
        power_difference = interferer_power - serving_power;
        
        % Set ICIC factor based on severity - using optimized values for private networks
        if power_difference > 3
            % Dominant interferer - very aggressive ICIC
            icic_factor = 0.05; % -13dB interference reduction
        elseif power_difference > -3
            % Strong interferer 
            icic_factor = 0.10; % -10dB interference reduction
        elseif power_difference > -10
            % Moderate interferer
            icic_factor = 0.20; % -7dB interference reduction
        else
            % Weak interferer
            icic_factor = 0.35; % -4.6dB interference reduction
        end
        
        % Apply relative ranking adjustment - first interferers get stronger ICIC
        % Small networks benefit from more aggressive treatment of primary interferers
        if i == 1 && numBS <= 4  % Primary interferer in small network
            icic_factor = icic_factor * 0.7;  % Even stronger reduction
        elseif i > (numBS/2)  % Less significant interferers
            icic_factor = icic_factor * 1.3;  % Less aggressive
        end
        
        % Apply coordinated beamforming effect if enabled
        if coordBeamforming_applied(interferingBS, ue)
            % More effective in small networks with fewer signal paths
            beam_factor = 0.4 + (0.1 * (i-1));  % 0.4 for primary, increasing for others
            icic_factor = icic_factor * beam_factor;
        end
        
        % Apply minimum ICIC factor based on network size
        if numBS <= 3
            min_icic_factor = 0.02;  % -17dB maximum for very small networks
        else
            min_icic_factor = 0.03;  % -15.2dB for larger networks
        end
        icic_factor = max(min_icic_factor, icic_factor);
        
        % Calculate interference contribution
        interference_contribution = interferer_power_mW * icic_factor;
        interference_power_mW = interference_power_mW + interference_contribution;
        
        % Store effective ICIC factor
        effective_icic_factors(interferingBS) = icic_factor;
    end
    
    % Adjust background interference based on network size
    background_interference_dBm = -95;
    background_interference_mW = db2pow(background_interference_dBm);
    interference_power_mW = interference_power_mW + background_interference_mW;
end


function [cqi, mcs, modOrder, codeRate, spectralEfficiency] = getLinkAdaptationParameters(sinr_dB)
    % Map SINR to CQI and determine modulation and coding parameters
    
    % Define SINR thresholds for CQI values (dB) - 3GPP based
    sinr_thresholds = [-6.7, -4.7, -2.3, 0.2, 2.4, 4.3, 5.9, 8.1, 10.3, 11.7, ...
                       14.1, 16.3, 18.7, 21.0, 22.7, 24.2];
    
    % Define MCS table including modulation order, code rate, and spectral efficiency
    % Format: [CQI, Modulation Order, Code Rate, Spectral Efficiency]
    mcs_table = [
        1,   2, 78/1024,  0.1523;  % QPSK, 78/1024
        2,   2, 120/1024, 0.2344;  % QPSK, 120/1024
        3,   2, 193/1024, 0.3770;  % QPSK, 193/1024
        4,   2, 308/1024, 0.6016;  % QPSK, 308/1024
        5,   2, 449/1024, 0.8770;  % QPSK, 449/1024
        6,   2, 602/1024, 1.1758;  % QPSK, 602/1024
        7,   4, 378/1024, 1.4766;  % 16QAM, 378/1024
        8,   4, 490/1024, 1.9141;  % 16QAM, 490/1024
        9,   4, 616/1024, 2.4063;  % 16QAM, 616/1024
        10,  6, 466/1024, 2.7305;  % 64QAM, 466/1024
        11,  6, 567/1024, 3.3223;  % 64QAM, 567/1024
        12,  6, 666/1024, 3.9023;  % 64QAM, 666/1024
        13,  6, 772/1024, 4.5234;  % 64QAM, 772/1024
        14,  6, 873/1024, 5.1152;  % 64QAM, 873/1024
        15,  8, 711/1024, 5.5547;  % 256QAM, 711/1024
        16,  8, 797/1024, 6.2266]; % 256QAM, 797/1024
    
    % Find the highest CQI that the current SINR can support
    cqi = 1; % Default to the lowest CQI
    for i = 1:length(sinr_thresholds)
        if sinr_dB >= sinr_thresholds(i)
            cqi = i;
        else
            break;
        end
    end
    
    % Cap CQI at 16 (maximum defined value)
    cqi = min(cqi, 16);
    
    % Get MCS parameters from the table
    mcs = mcs_table(cqi, 1);
    modOrder = mcs_table(cqi, 2);
    codeRate = mcs_table(cqi, 3);
    spectralEfficiency = mcs_table(cqi, 4);
end

function [selected_mcs, predicted_bler] = linkAdaptation(sinr_dB, target_bler)
    % Enhanced link adaptation with MCS selection based on target BLER
    
    % MCS options and their SINR thresholds for 10% BLER
    mcs_options = 0:28; % 5G NR MCS indices
    
    % SINR thresholds for 10% BLER (example values, would be calibrated in practice)
    % These values should increase with higher MCS as higher-order modulations need better SINR
    sinr_thresholds = [-7.5, -5.5, -3.5, -1.5, 0.5, 2.0, 3.5, 5.0, 6.5, 8.0, ...
                      9.5, 11.0, 12.5, 14.0, 15.5, 16.5, 17.5, 18.5, 19.5, ...
                      20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5];
    
    % Find highest MCS that meets SINR requirement
    selected_mcs = 0;
    for i = 1:length(sinr_thresholds)
        if sinr_dB >= sinr_thresholds(i)
            selected_mcs = mcs_options(i);
        else
            break;
        end
    end
    
    % Calculate predicted BLER for selected MCS
    if selected_mcs > 0
        mcs_idx = find(mcs_options == selected_mcs);
        sinr_diff = sinr_dB - sinr_thresholds(mcs_idx);
        % Simple BLER model: BLER decreases by factor of 10 for each 10dB SINR increase
        predicted_bler = target_bler * 10^(-sinr_diff/10);
        predicted_bler = min(1, max(0.001, predicted_bler));
    else
        predicted_bler = 1; % 100% BLER if no suitable MCS found
    end
end

function effective_throughput = calculateHARQThroughput(raw_throughput, sinr_dB, bler, max_retx)
    % Calculate effective throughput with HARQ retransmissions
    
    % Calculate success probability after all retransmissions
    success_prob = 0;
    for retx = 0:max_retx
        % Probability of success at this retransmission
        % For simplicity, we assume each retransmission has the same BLER
        % In practice, combining of retransmissions would improve BLER
        success_prob = success_prob + (1-success_prob) * (1-bler^(retx+1));
    end
    
    % Calculate effective throughput
    effective_throughput = raw_throughput * success_prob;
end

function [wtx, wrx, D] = getBeamformingWeights(hEst, nLayers, scOffset, noRBs)
% Enhanced function that computes beamforming weights using SVD
% with better support for multiple spatial layers
%
% hEst    : Channel estimate (matrix)
% nLayers : Number of layers to extract (spatial streams)
% scOffset: Offset from the first subcarrier
% noRBs   : Number of resource blocks over which to average

    % Get dimensions of channel estimate
    [~, ~, R, P] = size(hEst);

    % Select the relevant part of the channel estimate based on subcarrier offset
    scNo = scOffset + 1;
    hEstSubset = hEst(scNo:scNo + (12 * noRBs - 1), :, :, :);

    % Average the channel estimate over the selected subcarriers
    H = permute(mean(reshape(hEstSubset, [], R, P)), [2 3 1]);

    % Perform singular value decomposition
    [U, D, V] = svd(H);
    
    % Ensure D is a diagonal matrix with the right dimensions
    if ~isequal(size(D), [min(R,P), min(R,P)])
        D_temp = zeros(min(R,P), min(R,P));
        for i = 1:size(D,1)
            for j = 1:size(D,2)
                D_temp(i,j) = D(i,j);
            end
        end
        D = D_temp;
    end
    
    % Ensure we're not asking for more layers than are available
    actual_layers = min(nLayers, rank(H));
    if actual_layers < nLayers
        fprintf('Warning: Channel only supports %d layers, requested %d\n', actual_layers, nLayers);
    end
    
    % Extract beamforming weights for transmit and receive sides
    % For multiple layers, we use columns of V and U corresponding to the largest singular values
    wtx = V(:, 1:actual_layers).';
    wrx = U(:, 1:actual_layers).';
    
    % Return singular values for condition number analysis
    D = diag(D);
end

function optimal_layers = determineOptimalRank(hEst, max_layers, rank_sinr_threshold)
% Determines the optimal number of spatial layers based on channel conditions
%
% hEst               : Channel estimate (matrix)
% max_layers         : Maximum number of layers supported by hardware
% rank_sinr_threshold: SINR threshold (dB) for adding additional layers

    % Get channel matrix through SVD
    [~, ~, R, P] = size(hEst);
    scNo = 1;
    hEstSubset = hEst(scNo:scNo + 11, :, :, :); % Use first RB
    H = permute(mean(reshape(hEstSubset, [], R, P)), [2 3 1]);
    
    % Perform SVD to get singular values
    [~, S, ~] = svd(H);
    S = diag(S);
    
    % Calculate condition number for different ranks
    optimal_layers = 1; % Start with rank 1
    
    % Estimate SINR for each layer using singular values
    for layer = 2:min(max_layers, length(S))
        % Calculate layer power relative to first layer
        layer_power_db = 20*log10(S(layer)/S(1));
        
        % If the layer has sufficient power (SINR), include it
        if layer_power_db > -rank_sinr_threshold
            optimal_layers = layer;
        else
            break; % Stop at first layer that doesn't meet the threshold
        end
    end
    
    fprintf('Optimal number of layers based on channel conditions: %d\n', optimal_layers);
end

function TRX_mapping = createTRXtoAEMapping(numTRX, numAntElements)
% Creates a mapping matrix relating TRX units to antenna elements
%
% numTRX: Number of TRX units (64, 32 or 8)
% numAntElements: Number of antenna elements (192, 128 or 64)
%
% Returns: 
% TRX_mapping: A numTRX x numAntElements matrix where each row represents
%              a TRX and has 1's in positions corresponding to the antenna
%              elements it controls

    % Initialize mapping matrix with the correct dimensions
    TRX_mapping = zeros(numTRX, numAntElements);
    
    % Each TRX controls elementsPerTRX antenna elements
    elementsPerTRX = numAntElements / numTRX;
    
    % Assign each TRX to control consecutive elements
    for trx = 1:numTRX
        startElement = floor((trx-1) * elementsPerTRX) + 1;
        endElement = floor(trx * elementsPerTRX);
        TRX_mapping(trx, startElement:endElement) = 1;
    end

end

function wbs_trx = convertBeamWeightsToTRX(wbs_full, TRX_mapping, numTRX, numAntElements)
% Converts full antenna element weights to TRX control signals
%
% wbs_full: Full antenna element weight vector (1 x numAntElements)
% TRX_mapping: Mapping matrix (numTRX x numAntElements)
% numTRX: Number of TRX units (64, 32, or 8)
% numAntElements: Number of antenna elements (192, 128, or 64)
%
% Returns:
% wbs_trx: TRX-level weights (1 x numTRX)

    % Check that TRX_mapping has the correct dimensions
    [mappingRows, mappingCols] = size(TRX_mapping);
    if mappingRows < numTRX || mappingCols < numAntElements
        error(['TRX_mapping matrix dimensions (%d,%d) do not match ' ...
               'expected dimensions (%d,%d)'], ...
               mappingRows, mappingCols, numTRX, numAntElements);
    end

    % Initialize TRX weights
    wbs_trx = zeros(1, numTRX);
    
    % For each TRX, calculate average weight of its controlled elements
    for trx = 1:numTRX
        % Find antenna elements controlled by this TRX
        elements = find(TRX_mapping(trx, :));
        
        % Calculate average weight (both amplitude and phase)
        if ~isempty(elements)
            trx_elements_weights = wbs_full(elements);
            wbs_trx(trx) = mean(trx_elements_weights);
        end
    end
end

function wbs = applyTRXBeamWeights(wbs_trx, TRX_mapping, numAntElements)
% Converts TRX control signals back to full antenna element weights
%
% wbs_trx: TRX-level weights (1 x numTRX)
% TRX_mapping: Mapping matrix (numTRX x numAntElements)
% numAntElements: Number of antenna elements
%
% Returns:
% wbs: Full antenna element weight vector (1 x numAntElements)

    % Check that TRX_mapping has appropriate dimensions
    [numTRX, mappingCols] = size(TRX_mapping);
    if mappingCols < numAntElements
        error(['TRX_mapping matrix columns (%d) are less than ' ...
               'expected numAntElements (%d)'], ...
               mappingCols, numAntElements);
    end
    
    if length(wbs_trx) ~= numTRX
        error(['Length of wbs_trx (%d) does not match ' ...
               'number of TRX units in mapping (%d)'], ...
               length(wbs_trx), numTRX);
    end

    % Initialize full weight vector
    wbs = zeros(1, numAntElements);
    
    % For each antenna element, apply the weight from its controlling TRX
    for ae = 1:numAntElements
        % Find which TRX controls this element
        [trx, ~] = find(TRX_mapping(:, ae));
        
        % Apply TRX weight to this element
        if ~isempty(trx)
            wbs(ae) = wbs_trx(trx(1));  % Use first TRX if multiple (shouldn't happen)
        end
    end
end

function [wbs, wue, D] = getBeamformingWeightsWithError(hest, nLayers, scOffset, noRBs, error_factor)
    % Get the original weights
    [wbs_perfect, wue_perfect, D] = getBeamformingWeights(hest, nLayers, scOffset, noRBs);
    
   
    % Apply estimation error to the weights
    if nargin > 4 && error_factor > 0
        % Generate complex errors for each layer
        err_tx = complex(randn(size(wbs_perfect)) * error_factor, randn(size(wbs_perfect)) * error_factor);
        err_rx = complex(randn(size(wue_perfect)) * error_factor, randn(size(wue_perfect)) * error_factor);
        
        % Apply errors to weights
        wbs = wbs_perfect + err_tx;
        wue = wue_perfect + err_rx;
        
        % Normalize weights for each layer
        for i = 1:size(wbs, 1)
            wbs(i,:) = wbs(i,:) / norm(wbs(i,:));
        end
        for i = 1:size(wue, 1)
            wue(i,:) = wue(i,:) / norm(wue(i,:));
        end
    else
        wbs = wbs_perfect;
        wue = wue_perfect;
    end
end

function effective_throughput = calculateImprovedHARQThroughput(raw_throughput, sinr_dB, initial_bler, max_retx)
    % Calculate effective throughput with HARQ retransmissions using Chase combining
    
    % Initial BLER
    bler = initial_bler;
    
    % Calculate success probability after all retransmissions
    success_prob = 0;
    cumulative_success = 0;
    
    % For Chase combining, each retransmission improves SINR by ~3dB
    sinr_gain_per_retx = 3; % dB
    
    for retx = 0:max_retx
        if retx > 0
            % Improve SINR for retransmissions due to soft combining
            effective_sinr = sinr_dB + retx * sinr_gain_per_retx;
            
            % Recalculate BLER based on improved SINR
            [~, bler] = linkAdaptation(effective_sinr, initial_bler);
        end
        
        % Probability of success at this retransmission
        tx_success_prob = (1 - cumulative_success) * (1 - bler);
        success_prob = success_prob + tx_success_prob;
        cumulative_success = success_prob;
    end
    
    % Calculate effective throughput, accounting for retransmission overhead
    avg_tx_count = 0;
    for retx = 0:max_retx
        if retx < max_retx
            prob_exactly_retx = initial_bler^retx * (1-initial_bler);
        else
            prob_exactly_retx = initial_bler^max_retx;
        end
        avg_tx_count = avg_tx_count + prob_exactly_retx * (retx + 1);
    end
    
    % Throughput reduction factor due to retransmissions
    retx_overhead = 1 / avg_tx_count;
    
    % Final effective throughput
    effective_throughput = raw_throughput * success_prob * retx_overhead;
end


function [coordinated_weights] = calculateCoordinatedBeamWeights(hest, serving_bs, ue, allChannels, validConnections, rtChannels, rxPower)
    % This function calculates coordinated beamforming weights that create nulls toward other UEs
    
    % Get the number of spatial layers this UE was assigned
    if exist('allNumLayers', 'var') && ~isempty(allNumLayers{serving_bs, ue})
        nLayers = allNumLayers{serving_bs, ue};
    else
        % If not available, determine from channel dimensions
        [~, ~, R, P] = size(hest);
        nLayers = min(R, P);
        nLayers = min(nLayers, 4); % Cap at 4 layers max
    end
    
    scOffset = 0;
    noRBs = 1;
    [wtx_init, ~, ~] = getBeamformingWeights(hest, nLayers, scOffset, noRBs);
    
    % Initialize coordinated weights with initial weights
    coordinated_weights = wtx_init;
    
    % Get the number of BS and UEs
    numBS = size(validConnections, 1);
    numUEs = size(validConnections, 2);
    
    % These are the UEs that we want to nullify interference towards
    % (UEs that have stronger connection to other BSs)
    interfering_ues = [];
    
    % Find UEs that are better served by other BSs
    for other_ue = 1:numUEs
        if other_ue ~= ue && validConnections(serving_bs, other_ue)
            % Check if this UE is better served by another BS
            is_better_elsewhere = false;
            
            for bs = 1:numBS
                if bs ~= serving_bs && validConnections(bs, other_ue) && rxPower(bs, other_ue) > rxPower(serving_bs, other_ue)
                    is_better_elsewhere = true;
                    break;
                end
            end
            
            if is_better_elsewhere
                interfering_ues = [interfering_ues; other_ue];
            end
        end
    end
    
    % If we found UEs to protect from interference
    if ~isempty(interfering_ues)
        % Get dimensions from initial weights to ensure compatibility
        [rows, cols] = size(coordinated_weights);
        
        % Use a simple but effective interference reduction approach
        % that doesn't depend on details of the channel structure
        
        % For each UE we want to protect
        for i = 1:length(interfering_ues)
            int_ue = interfering_ues(i);
            
            % Simple interference mitigation: 
            % Just slightly tilt the beam in a different direction
            
            % We're making a small 5% modification to the weights in a deterministic way
            tilt_factor = 0.05;
            phase_shift = pi/4 * (i/length(interfering_ues));
            
            % Apply a deterministic phase shift based on UE index
            for r = 1:rows
                for c = 1:cols
                    old_weight = coordinated_weights(r,c);
                    if old_weight ~= 0
                        angle_shift = phase_shift * (int_ue / numUEs);
                        old_mag = abs(old_weight);
                        old_phase = angle(old_weight);
                        new_phase = old_phase + angle_shift;
                        coordinated_weights(r,c) = old_mag * exp(1i * new_phase);
                    end
                end
            end
        end
        
        % Normalize the weights per layer
        for r = 1:rows
            if norm(coordinated_weights(r,:)) > 0
                coordinated_weights(r,:) = coordinated_weights(r,:) / norm(coordinated_weights(r,:));
            end
        end
    end
end

function [selected_mcs, predicted_bler] = linkAdaptationWithCQI(sinr_dB, target_bler, cqi)
    % Enhanced BLER calculation for different CQI values
    % Apply different SINR offsets based on CQI
    switch cqi
        case 1  % QPSK, lowest rate
            sinr_offset = 0;
        case 7  % 16-QAM, medium rate
            sinr_offset = 8;
        case 10 % 64-QAM, high rate
            sinr_offset = 13;
        case 15 % 256-QAM, highest rate
            sinr_offset = 18;
        otherwise
            sinr_offset = 0;
    end
    
    % Apply offset to simulate different curves for different CQIs
    effective_sinr = sinr_dB - sinr_offset;
    
    % Simple BLER model
    predicted_bler = 1 / (1 + exp(effective_sinr/1.5));
    predicted_bler = min(1, max(0.001, predicted_bler));
    
    % Set MCS based on CQI (simplified mapping)
    selected_mcs = cqi * 2;
end



%% Group 2: Network Performance Metrics
figure('Name', 'Network Performance Metrics - Serving Connections Only', 'Position', [50, 50, 1000, 800]);

% CDF of throughput values (serving connections only)
subplot(2, 3, 1);
if ~isempty(all_throughputs)
    [f, x] = ecdf(all_throughputs);
    plot(x, f, 'LineWidth', 2);
    grid on;
    xlabel('Throughput (Mbps)');
    ylabel('Cumulative Probability');
    title('CDF of Throughput');
    
    % Add vertical line for median throughput
    hold on;
    median_tput = median(all_throughputs);
    plot([median_tput, median_tput], [0, 0.5], 'r--', 'LineWidth', 1.5);
    plot([0, median_tput], [0.5, 0.5], 'r--', 'LineWidth', 1.5);
    text(median_tput*1.1, 0.52, ['Median: ' num2str(median_tput, '%.1f') ' Mbps'], 'Color', 'r');
    hold off;
else
    text(0.5, 0.5, 'No valid connections', 'HorizontalAlignment', 'center');
    axis off;
end

% CQI distribution histogram (serving connections only)
subplot(2, 3, 2);
if ~isempty(all_cqis)
    histogram(all_cqis, 1:16, 'FaceColor', 'b', 'EdgeColor', 'k');
    grid on;
    xlabel('CQI Value');
    ylabel('Number of Serving Links');
    title('CQI Distribution');
    xticks(1:15);
else
    text(0.5, 0.5, 'No valid connections', 'HorizontalAlignment', 'center');
    axis off;
end

% Modulation distribution pie chart (serving connections only)
subplot(2, 3, 3);
if sum(mod_counts) > 0
    pie(mod_counts(non_zero_idx), mod_types(non_zero_idx));
    title('Modulation Distribution');
else
    text(0.5, 0.5, 'No valid connections', 'HorizontalAlignment', 'center');
    axis off;
end

% Network capacity analysis
subplot(2, 3, 4);
bar(total_capacity_per_bs);
hold on;
yline(total_network_capacity/numBS, 'r--', 'LineWidth', 2);
text(0.5, total_network_capacity/numBS*1.1, ['Avg: ' num2str(total_network_capacity/numBS, '%.1f') ' Mbps']);
hold off;
xlabel('Base Station ID');
ylabel('Total Capacity (Mbps)');
title(['Network Capacity: ' num2str(total_network_capacity, '%.1f') ' Mbps']);
grid on;
xticks(1:numBS);


% Load balancing visualization
subplot(2, 3, 6);
% Calculate statistics
ue_per_bs_std = std(ues_per_bs);
max_ue_per_bs = max(ues_per_bs);
min_ue_per_bs = min(ues_per_bs);
% Create a visual representation of load balance
bar(ues_per_bs);
hold on;
yline(mean(ues_per_bs), 'r--', 'LineWidth', 2);
text(0.5, mean(ues_per_bs)*1.1, ['Avg: ' num2str(mean(ues_per_bs), '%.1f') ' UEs']);
hold off;
title(['Load Balance (StdDev: ' num2str(ue_per_bs_std, '%.2f') ' UEs)']);
xlabel('Base Station ID');
ylabel('Number of UEs');
grid on;
xticks(1:numBS);

% Add overall title
sgtitle(['Network Performance - ' char(antennaConfig) ' Configuration (Serving Only)'], 'FontSize', 16, 'FontWeight', 'bold');

%% Group 3: Link Adaptation and AMC Performance
figure('Name', 'Link Adaptation & AMC Performance - Serving Connections Only', 'Position', [50, 50, 1000, 800]);

% CQI vs SINR scatter plot with modulation regions (serving connections only)
subplot(2, 2, 1);
hold on;

% Plot data points for each serving UE
for ue = 1:numUEs
    if best_serving_bs(ue) > 0
        bs = best_serving_bs(ue);
        scatter(sinrValuesCapped(bs, ue), cqiValues(bs, ue), 50, 'filled');
    end
end

% Add labels for modulation regions
text(-5, 3, 'QPSK', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'blue');
text(7, 8, '16-QAM', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'green');
text(16, 12, '64-QAM', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
text(25, 15, '256-QAM', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'magenta');

grid on;
xlabel('SINR (dB)');
ylabel('CQI Value');
title('CQI vs SINR - Serving Connections');
xlim([-10 30]);
ylim([0 17]);

% Add colorful background regions for different modulations
patch([-10 -10 5.9 5.9], [0 6.5 6.5 0], 'blue', 'FaceAlpha', 0.1);  % QPSK
patch([5.9 5.9 14.1 14.1], [0 9.5 9.5 0], 'green', 'FaceAlpha', 0.1);  % 16-QAM
patch([14.1 14.1 22.7 22.7], [0 14.5 14.5 0], 'red', 'FaceAlpha', 0.1);  % 64-QAM
patch([22.7 22.7 30 30], [0 16.5 16.5 0], 'magenta', 'FaceAlpha', 0.1);  % 256-QAM

hold off;

% Throughput vs SINR curve with actual data points (serving connections only)
subplot(2, 2, 2);
hold on;

% In the "Throughput vs SINR curve with actual data points" subplot section (subplot 2,2,2)

% Plot theoretical Shannon capacity curves for different layers
sinr_range = -10:0.1:30;  % dB
shannon_capacity = log2(1 + 10.^(sinr_range/10));  % bps/Hz

% Define layer efficiencies - these determine the scaling for each layer
layer_efficiency = [1.0, 0.85, 0.7, 0.55]; % Efficiency factors for layers 1-4
layer_colors = {'k-', 'b--', 'g:', 'r-.'};  % Different line styles for each layer
legend_entries = {};

% Plot Shannon capacity for each layer
for layer = 1:min(max_spatial_layers, 4)
    if layer == 1
        % First layer - base capacity
        cumulative_eff = layer_efficiency(1);
        shannon_throughput = shannon_capacity * BW / 1e6 * cumulative_eff;
        plot(sinr_range, shannon_throughput, layer_colors{1}, 'LineWidth', 2);
        legend_entries{end+1} = 'Shannon (1 layer)';
    else
        % Additional layers - show cumulative capacity
        cumulative_eff = sum(layer_efficiency(1:layer));
        shannon_throughput = shannon_capacity * BW / 1e6 * cumulative_eff;
        plot(sinr_range, shannon_throughput, layer_colors{layer}, 'LineWidth', 2);
        legend_entries{end+1} = sprintf('Shannon (%d layers)', layer);
    end
end

% Throughput vs SINR curve with actual data points (serving connections only)

valid_sinrs = [];
valid_tputs = [];
cqi_data = [];

% Only collect data for UEs with valid connections to their best serving BS
for ue = 1:numUEs
    if best_serving_bs(ue) > 0
        bs = best_serving_bs(ue);
        % Make sure this is a valid connection
        if validConnections(bs, ue)
            valid_sinrs = [valid_sinrs, sinrValuesCapped(bs, ue)];
            valid_tputs = [valid_tputs, effectiveThroughput(bs, ue)];
            cqi_data = [cqi_data, cqiValues(bs, ue)];
        end
    end
end

% Check if we have data to plot
if ~isempty(valid_sinrs) && length(valid_sinrs) == length(valid_tputs) && length(valid_tputs) == length(cqi_data)
    % Create scatter plot with CQI-based coloring
    scatter(valid_sinrs, valid_tputs, 80, cqi_data, 'filled', 'MarkerEdgeColor', 'k');
    cb = colorbar;
    ylabel(cb, 'CQI Value');
    
    % Apply colormap safely to current axes only
    try
        colormap(gca, jet);
    catch
        % If jet fails, MATLAB's default colormap
        try
            colormap(gca, parula);
        catch
            % Last resort - don't change colormap
            fprintf('Warning: Could not apply custom colormap\n');
        end
    end
    
    % Add UEs to legend entries
    legend_entries{end+1} = 'Serving UEs (color = CQI)';
else
    % No valid data to plot
    text(10, max(shannon_throughput)/2, 'No valid UE connections to display', 'FontSize', 12, 'HorizontalAlignment', 'center');
    legend_entries{end+1} = 'No valid UE data';
end

legend(legend_entries);
title('Throughput vs SINR - Serving Connections');
xlabel('SINR (dB)');
ylabel('Throughput (Mbps)');
grid on;
hold off;

% BLER vs SINR analysis
subplot(2, 2, 3);
hold on;
% Plot BLER for different CQI values
sinr_range = -10:0.5:30;
colors = {'b', 'g', 'r', 'm'}; % Colors for CQI 1, 7, 10, 15
cqi_values = [1, 7, 10, 15];  % Representative CQI values

for idx = 1:length(cqi_values)
    cqi = cqi_values(idx);
    bler_values = zeros(size(sinr_range));
    
    for i = 1:length(sinr_range)
        [~, bler_values(i)] = linkAdaptationWithCQI(sinr_range(i), target_bler, cqi);
    end
    
    plot(sinr_range, bler_values, colors{idx}, 'LineWidth', 2);
end

set(gca, 'YScale', 'log');
ylim([0.001 1]);
grid on;
xlabel('SINR (dB)');
ylabel('Block Error Rate (BLER)');
title('BLER vs SINR for Different CQI Values');
legend('CQI 1 (QPSK)', 'CQI 7 (16-QAM)', 'CQI 10 (64-QAM)', 'CQI 15 (256-QAM)');
hold off;

% Spectral efficiency analysis (serving connections only)
subplot(2, 2, 4);
hold on;
% Scatter plot of spectral efficiency vs SINR
se_values = [];
for ue = 1:numUEs
    if best_serving_bs(ue) > 0
        bs = best_serving_bs(ue);
        se_values = [se_values; sinrValuesCapped(bs, ue), spectralEfficiency(bs, ue)];
    end
end
if ~isempty(se_values)
    scatter(se_values(:,1), se_values(:,2), 60, 'filled');
    
    % Add Shannon capacity curve for reference
    sinr_range = -10:0.5:30;
    shannon_se = log2(1 + 10.^(sinr_range/10));
    plot(sinr_range, shannon_se, 'r-', 'LineWidth', 2);
    
    xlabel('SINR (dB)');
    ylabel('Spectral Efficiency (bps/Hz)');
    title('Spectral Efficiency vs SINR - Serving');
    legend('Serving UEs', 'Shannon Limit');
    grid on;
else
    text(0.5, 0.5, 'No valid connections', 'HorizontalAlignment', 'center');
    axis off;
end

%% Save Simulation Results
% Create a structure with all results
resultsToSave = struct();
resultsToSave.antennaConfig = antennaConfig;
resultsToSave.sinrValuesCapped = sinrValuesCapped;
resultsToSave.effectiveThroughput = effectiveThroughput;
resultsToSave.cqiValues = cqiValues;
resultsToSave.blerEstimates = blerEstimates;
resultsToSave.rxPower = rxPower;
resultsToSave.allBeamformingGains = allBeamformingGains;
resultsToSave.validConnections = validConnections;
resultsToSave.numBS = numBS;
resultsToSave.numUEs = numUEs;
resultsToSave.allNumLayers = allNumLayers;
resultsToSave.max_spatial_layers = max_spatial_layers;
resultsToSave.spatial_multiplexing_enabled = spatial_multiplexing_enabled;

% Save the results with the configuration name in the filename
filename = ['results_' char(antennaConfig) '.mat'];
save(filename, 'resultsToSave');
fprintf('Results saved to %s\n', filename);
