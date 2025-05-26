close all
clear variables

% Author: Tomasz Syrylo
%% Setup Parameters
% Carrier frequency (Hz)
fc = 25e9;  % 4 GHz

% Ray tracing configuration 
reflectionsOrder = 2;                   % Set number of reflections here (0 for LOS only)
diffractionOrder = 0;                   % max 1!! , default 0

% FIRST CONFIG
% Base station parameters
bsPosition1 = [52.251387,20.905994];    % [Latitude, Longitude]
bsArrayOrientation1 = [120 -10].';     % [Azimuth; Elevation] (degrees)
bsPosition2 = [52.256474,20.901313];    % [Latitude, Longitude]
bsArrayOrientation2 = [-140 -5].';     % [Azimuth; Elevation] (degrees)
bsPosition3 = [52.251076,20.899524];    % [Latitude, Longitude]
bsArrayOrientation3 = [40 -5].';       % [Azimuth; Elevation] (degrees)

% SECOND CONFIG
% bsPosition1 = [52.256474,20.901313];    % [Latitude, Longitude]
% bsArrayOrientation1 = [-120 -5].';     % [Azimuth; Elevation] (degrees)
% bsPosition2 = [52.251076,20.899524];    % [Latitude, Longitude]
% bsArrayOrientation2 = [-260 -5].';     % [Azimuth; Elevation] (degrees)

% Base station transmit power (dBm)
bsTxPower = 43; % 5W = 37dBm 20W = 43dBm

% Antenna configuration selection
% Set antennaConfig to "32TRX_128AE", or "8TRX_64AE"
antennaConfig = "32TRX_128AE"; 


% Calculate wavelength for antenna elements
c = physconst('LightSpeed');
lambda = c / fc;

% Create antenna array based on selected configuration
if strcmpi(antennaConfig, "32TRX_128AE")
    bsArray = phased.NRRectangularPanelArray(...
        'Size', [8 8 1 1], ...
        'Spacing', 0.5 * lambda * [1 1 1 1]);
    bsArray.ElementSet = {phased.NRAntennaElement('PolarizationAngle', -45), ...
                         phased.NRAntennaElement('PolarizationAngle', 45)};
    numAntElements = prod(bsArray.Size) * length(bsArray.ElementSet); % 128
    numTRX = 32;
    beamformingGain = 10*log10(numAntElements) - 10; % Reduced gain by increasing loss
    
elseif strcmpi(antennaConfig, "8TRX_64AE")
    bsArray = phased.NRRectangularPanelArray(...
        'Size', [8 4 1 1], ...
        'Spacing', 0.5 * lambda * [1 1 1 1]);
    bsArray.ElementSet = {phased.NRAntennaElement('PolarizationAngle', -45), ...
                         phased.NRAntennaElement('PolarizationAngle', 45)};
    numAntElements = prod(bsArray.Size) * length(bsArray.ElementSet); % 64
    numTRX = 8;
    beamformingGain = 10*log10(numAntElements) - 10; % Reduced gain by increasing loss
end

% Calculate effective transmit power including beamforming gain
effectiveTxPower = bsTxPower;% + beamformingGain;
effectiveTxPower_w = db2pow(effectiveTxPower - 30);
%% Viewer Initialization
% Create the map viewer
viewer = siteviewer("Basemap","openstreetmap","Buildings","wat_expanded.osm"); 

%% Define Base Station Transmitter Sites
bs1 = txsite("Name", "Base Station 1", ...
    "Latitude", bsPosition1(1), ...
    "Longitude", bsPosition1(2), ...
    "AntennaHeight", 4, ...
    "AntennaAngle", bsArrayOrientation1(1:2), ...
    "TransmitterPower",effectiveTxPower_w , ... 
    "TransmitterFrequency", fc);

bs2 = txsite("Name", "Base Station 2", ...
    "Latitude", bsPosition2(1), ...
    "Longitude", bsPosition2(2), ...
    "AntennaHeight", 4, ...
    "AntennaAngle", bsArrayOrientation2(1:2), ...
    "TransmitterPower", effectiveTxPower_w, ...
    "TransmitterFrequency", fc);

bs3 = txsite("Name", "Base Station 3", ...
    "Latitude", bsPosition3(1), ...
    "Longitude", bsPosition3(2), ...
    "AntennaHeight", 4, ...
    "AntennaAngle", bsArrayOrientation3(1:2), ...
    "TransmitterPower", effectiveTxPower_w, ...
    "TransmitterFrequency", fc);

% Store all transmitters in array
bss = [bs1, bs2, bs3]; %, bs3

% Show all transmitters
for i = 1:length(bss)
    show(bss(i));
end

%% Create a propagation model with ray tracing (SBR method)
rtpm = propagationModel("raytracing", ...
    "Method", "sbr", ...
    "MaxNumReflections", reflectionsOrder, ...
    "MaxNumDiffractions", diffractionOrder, ...
    "BuildingsMaterial", "concrete", ...
    "TerrainMaterial", "concrete");

% Add weather effects
rtPlusWeather = rtpm   + propagationModel("rain") + propagationModel("gas"); % about 1.5 dB of loss

%% Coverage Parameters
signalStrengths = -120:-5;  % Signal strength range
maxRange = 500;  % meters
resolution =3;  % meters

% Display Coverage Map
fprintf('Displaying coverage map with %d reflections...\n', reflectionsOrder);
coverage(bss, rtPlusWeather, ...
    "SignalStrengths", signalStrengths, ...
    "MaxRange", maxRange, ...
    "Resolution", resolution, ...
    "Transparency", 0.6);

% Wait for user input before continuing
fprintf('Coverage map displayed. Press any key to continue to SINR map...\n');
pause;

%% Display SINR Map
% Clear only the coverage map, not the entire viewer
hide(bss(1)); % Temporarily hide and show to refresh without clearing
show(bss(1));

% Calculate and display SINR map
sinr(bss, rtPlusWeather, ...
    "MaxRange", maxRange, ...
    "Resolution", resolution, ...
    "Transparency", 0.6);

fprintf('SINR map displayed.\n');

%% Display System Information
fprintf('\n*** System Configuration ***\n');
fprintf('Number of Reflections: %d\n', reflectionsOrder);
fprintf('Antenna Configuration: %s\n', antennaConfig);
fprintf('Number of Antenna Elements: %d\n', numAntElements);
fprintf('Beamforming Gain: %.1f dB\n', beamformingGain);
fprintf('Effective Transmit Power: %.1f dBm\n', effectiveTxPower);