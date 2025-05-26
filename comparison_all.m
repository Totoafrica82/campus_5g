% Author: Tomasz Syrylo
%% MIMO Configuration Comparison Script
% This script loads and compares results from different antenna configurations
% Modified to only display values for serving UE-BS pairs

% Load the saved results
results_8TRX = load('results_8TRX_64AE.mat').resultsToSave;
results_32TRX = load('results_32TRX_128AE.mat').resultsToSave;
results_64TRX = load('results_64TRX_192AE.mat').resultsToSave;

% Define colors for consistent visualization
colors_8TRX = [0.3 0.5 0.8];  % Blue
colors_32TRX = [0.8 0.3 0.3]; % Red
colors_64TRX = [0.3 0.8 0.3]; % Green

% Create serving pair masks for each configuration
serving_8TRX = false(size(results_8TRX.validConnections));
serving_32TRX = false(size(results_32TRX.validConnections));
serving_64TRX = false(size(results_64TRX.validConnections));

% Build serving pairs masks based on best performance
% For a more accurate representation, we'll use effective throughput to determine the serving BS
serving_8TRX = false(size(results_8TRX.validConnections));
serving_32TRX = false(size(results_32TRX.validConnections));
serving_64TRX = false(size(results_64TRX.validConnections));

% For each UE, find which BS provides the highest throughput
for ue = 1:results_8TRX.numUEs
    [max_val, best_bs] = max(results_8TRX.effectiveThroughput(:, ue));
    if max_val > 0 && results_8TRX.validConnections(best_bs, ue)
        serving_8TRX(best_bs, ue) = true;
    end
end

for ue = 1:results_32TRX.numUEs
    [max_val, best_bs] = max(results_32TRX.effectiveThroughput(:, ue));
    if max_val > 0 && results_32TRX.validConnections(best_bs, ue)
        serving_32TRX(best_bs, ue) = true;
    end
end

for ue = 1:results_64TRX.numUEs
    [max_val, best_bs] = max(results_64TRX.effectiveThroughput(:, ue));
    if max_val > 0 && results_64TRX.validConnections(best_bs, ue)
        serving_64TRX(best_bs, ue) = true;
    end
end

% Create a comparison figure for throughput analysis
figure('Name', 'Throughput Performance Comparison', 'Position', [100, 100, 1200, 800]);

%% 1. Peak Throughput Comparison
subplot(2, 2, 1);

% Extract valid throughput values for serving pairs only
tput_8TRX = [];
tput_32TRX = [];
tput_64TRX = [];

% Extract throughput values from serving BS-UE pairs
for bs = 1:results_8TRX.numBS
    for ue = 1:results_8TRX.numUEs
        if serving_8TRX(bs, ue)
            tput_8TRX = [tput_8TRX; results_8TRX.effectiveThroughput(bs, ue)];
        end
        
        if bs <= results_32TRX.numBS && ue <= results_32TRX.numUEs && serving_32TRX(bs, ue)
            tput_32TRX = [tput_32TRX; results_32TRX.effectiveThroughput(bs, ue)];
        end
        
        if bs <= results_64TRX.numBS && ue <= results_64TRX.numUEs && serving_64TRX(bs, ue)
            tput_64TRX = [tput_64TRX; results_64TRX.effectiveThroughput(bs, ue)];
        end
    end
end

% Calculate peak throughput for each configuration
peak_tput_8TRX = max(tput_8TRX);
peak_tput_32TRX = max(tput_32TRX);
peak_tput_64TRX = max(tput_64TRX);

% Create bar chart for peak throughput
bar_handle = bar([peak_tput_8TRX, peak_tput_32TRX, peak_tput_64TRX]);
set(bar_handle, 'FaceColor', 'flat');
bar_handle.CData(1,:) = colors_8TRX;
bar_handle.CData(2,:) = colors_32TRX;
bar_handle.CData(3,:) = colors_64TRX;

set(gca, 'XTickLabel', {'8TRX_64AE', '32TRX_128AE', '64TRX_192AE'});
ylabel('Peak Throughput (Mbps)');
title('Peak Throughput Comparison');
grid on;

% Add text with values
text(1, peak_tput_8TRX + 10, sprintf('%.1f Mbps', peak_tput_8TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(2, peak_tput_32TRX + 10, sprintf('%.1f Mbps', peak_tput_32TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(3, peak_tput_64TRX + 10, sprintf('%.1f Mbps', peak_tput_64TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');

% Calculate improvement percentages
improvement_32vs8 = ((peak_tput_32TRX - peak_tput_8TRX) / peak_tput_8TRX) * 100;
improvement_64vs32 = ((peak_tput_64TRX - peak_tput_32TRX) / peak_tput_32TRX) * 100;
improvement_64vs8 = ((peak_tput_64TRX - peak_tput_8TRX) / peak_tput_8TRX) * 100;

% Add improvement percentages as annotations
if abs(improvement_32vs8) > 0.1
    text(1.5, max([peak_tput_8TRX, peak_tput_32TRX]) * 0.85, sprintf('%.1f%%', improvement_32vs8), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10, ...
        'BackgroundColor', [0.9 0.9 0.9]);
end

if abs(improvement_64vs32) > 0.1
    text(2.5, max([peak_tput_32TRX, peak_tput_64TRX]) * 0.85, sprintf('%.1f%%', improvement_64vs32), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10, ...
        'BackgroundColor', [0.9 0.9 0.9]);
end

if abs(improvement_64vs8) > 0.1
    text(2, max([peak_tput_8TRX, peak_tput_64TRX]) * 0.7, sprintf('64TRX vs 8TRX: %.1f%%', improvement_64vs8), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10, ...
        'BackgroundColor', [0.9 0.9 0.9]);
end

%% 2. Average Throughput Comparison
subplot(2, 2, 2);

% Calculate average throughput for each configuration
avg_tput_8TRX = mean(tput_8TRX);
avg_tput_32TRX = mean(tput_32TRX);
avg_tput_64TRX = mean(tput_64TRX);

% Create bar chart for average throughput
bar_handle = bar([avg_tput_8TRX, avg_tput_32TRX, avg_tput_64TRX]);
set(bar_handle, 'FaceColor', 'flat');
bar_handle.CData(1,:) = colors_8TRX;
bar_handle.CData(2,:) = colors_32TRX;
bar_handle.CData(3,:) = colors_64TRX;

set(gca, 'XTickLabel', {'8TRX_64AE', '32TRX_128AE', '64TRX_192AE'});
ylabel('Average Throughput (Mbps)');
title('Average Throughput Comparison');
grid on;

% Add text with values
text(1, avg_tput_8TRX + 5, sprintf('%.1f Mbps', avg_tput_8TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(2, avg_tput_32TRX + 5, sprintf('%.1f Mbps', avg_tput_32TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(3, avg_tput_64TRX + 5, sprintf('%.1f Mbps', avg_tput_64TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');

% Calculate improvement percentages
avg_improvement_32vs8 = ((avg_tput_32TRX - avg_tput_8TRX) / avg_tput_8TRX) * 100;
avg_improvement_64vs32 = ((avg_tput_64TRX - avg_tput_32TRX) / avg_tput_32TRX) * 100;
avg_improvement_64vs8 = ((avg_tput_64TRX - avg_tput_8TRX) / avg_tput_8TRX) * 100;

% Add improvement percentages as annotations
if abs(avg_improvement_32vs8) > 0.1
    text(1.5, max([avg_tput_8TRX, avg_tput_32TRX]) * 0.85, sprintf('%.1f%%', avg_improvement_32vs8), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10, ...
        'BackgroundColor', [0.9 0.9 0.9]);
end

if abs(avg_improvement_64vs32) > 0.1
    text(2.5, max([avg_tput_32TRX, avg_tput_64TRX]) * 0.85, sprintf('%.1f%%', avg_improvement_64vs32), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10, ...
        'BackgroundColor', [0.9 0.9 0.9]);
end

if abs(avg_improvement_64vs8) > 0.1
    text(2, max([avg_tput_8TRX, avg_tput_64TRX]) * 0.7, sprintf('64TRX vs 8TRX: %.1f%%', avg_improvement_64vs8), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10, ...
        'BackgroundColor', [0.9 0.9 0.9]);
end

%% 3. Throughput Distribution
subplot(2, 2, 3);

% Create histogram with equal bins for fair comparison
bin_edges = linspace(0, max([tput_8TRX; tput_32TRX; tput_64TRX]) * 1.05, 15);

% Create histograms using same bins
h1 = histogram(tput_8TRX, bin_edges, 'FaceColor', colors_8TRX, 'FaceAlpha', 0.7);
hold on;
h2 = histogram(tput_32TRX, bin_edges, 'FaceColor', colors_32TRX, 'FaceAlpha', 0.7);
h3 = histogram(tput_64TRX, bin_edges, 'FaceColor', colors_64TRX, 'FaceAlpha', 0.7);

% Calculate percentage of connections above 500 Mbps
percent_high_tput_8TRX = 100 * sum(tput_8TRX > 500) / length(tput_8TRX);
percent_high_tput_32TRX = 100 * sum(tput_32TRX > 500) / length(tput_32TRX);
percent_high_tput_64TRX = 100 * sum(tput_64TRX > 500) / length(tput_64TRX);

% Calculate percentage of connections below 50 Mbps
percent_low_tput_8TRX = 100 * sum(tput_8TRX < 50) / length(tput_8TRX);
percent_low_tput_32TRX = 100 * sum(tput_32TRX < 50) / length(tput_32TRX);
percent_low_tput_64TRX = 100 * sum(tput_64TRX < 50) / length(tput_64TRX);

% Add legend with percentages
legend(['8TRX: ' num2str(percent_high_tput_8TRX, '%.1f') '% > 500 Mbps'], ...
       ['32TRX: ' num2str(percent_high_tput_32TRX, '%.1f') '% > 500 Mbps'], ...
       ['64TRX: ' num2str(percent_high_tput_64TRX, '%.1f') '% > 500 Mbps'], ...
       'Location', 'northwest');

grid on;
xlabel('Throughput (Mbps)');
ylabel('Number of BS-UE Connections');
title('Throughput Distribution');

%% 4. SINR Comparison with Cell-Edge Focus
subplot(2, 2, 4);

% Extract SINR values for serving pairs only
sinr_8TRX = [];
sinr_32TRX = [];
sinr_64TRX = [];

% Extract SINR values from serving BS-UE pairs
for bs = 1:results_8TRX.numBS
    for ue = 1:results_8TRX.numUEs
        if serving_8TRX(bs, ue)
            sinr_8TRX = [sinr_8TRX; results_8TRX.sinrValuesCapped(bs, ue)];
        end
        
        if bs <= results_32TRX.numBS && ue <= results_32TRX.numUEs && serving_32TRX(bs, ue)
            sinr_32TRX = [sinr_32TRX; results_32TRX.sinrValuesCapped(bs, ue)];
        end
        
        if bs <= results_64TRX.numBS && ue <= results_64TRX.numUEs && serving_64TRX(bs, ue)
            sinr_64TRX = [sinr_64TRX; results_64TRX.sinrValuesCapped(bs, ue)];
        end
    end
end

% Create a CDF plot (empirical cumulative distribution function)
[f_8TRX, x_8TRX] = ecdf(sinr_8TRX);
[f_32TRX, x_32TRX] = ecdf(sinr_32TRX);
[f_64TRX, x_64TRX] = ecdf(sinr_64TRX);

plot(x_8TRX, f_8TRX, 'Color', colors_8TRX, 'LineWidth', 2);
hold on;
plot(x_32TRX, f_32TRX, 'Color', colors_32TRX, 'LineWidth', 2);
plot(x_64TRX, f_64TRX, 'Color', colors_64TRX, 'LineWidth', 2);

% Calculate cell-edge SINR (10th percentile)
edge_sinr_8TRX = prctile(sinr_8TRX, 10);
edge_sinr_32TRX = prctile(sinr_32TRX, 10);
edge_sinr_64TRX = prctile(sinr_64TRX, 10);

% Calculate median SINR (50th percentile)
median_sinr_8TRX = prctile(sinr_8TRX, 50);
median_sinr_32TRX = prctile(sinr_32TRX, 50);
median_sinr_64TRX = prctile(sinr_64TRX, 50);

% Add vertical lines for cell-edge SINR
xline(edge_sinr_8TRX, 'Color', colors_8TRX, 'LineStyle', '--', 'LineWidth', 1.5, ...
      'Label', [' 8TRX Edge: ' num2str(edge_sinr_8TRX, '%.1f') ' dB']);
xline(edge_sinr_32TRX, 'Color', colors_32TRX, 'LineStyle', '--', 'LineWidth', 1.5, ...
      'Label', [' 32TRX Edge: ' num2str(edge_sinr_32TRX, '%.1f') ' dB']);
xline(edge_sinr_64TRX, 'Color', colors_64TRX, 'LineStyle', '--', 'LineWidth', 1.5, ...
      'Label', [' 64TRX Edge: ' num2str(edge_sinr_64TRX, '%.1f') ' dB']);

% Add horizontal lines at 10% and 50% for easier reading
yline(0.1, 'k:', '10% (Cell Edge)');
yline(0.5, 'k:', '50% (Median)');

grid on;
xlabel('SINR (dB)');
ylabel('CDF (Probability)');
title('SINR Cumulative Distribution Function (Serving Pairs)');
legend('8TRX_64AE', '32TRX_128AE', '64TRX_192AE', 'Location', 'northwest');

% Add overall title
sgtitle('MIMO Configuration Performance Comparison (Serving Pairs Only)', 'FontSize', 14, 'FontWeight', 'bold');

%% Create a second figure focusing on signal power and beamforming
figure('Name', 'Signal Power and Beamforming Comparison', 'Position', [100, 100, 1200, 400]);

%% 1. Beamforming Gain Comparison
subplot(1, 3, 1);

% Extract valid beamforming gain values for serving pairs only
bf_gain_8TRX = [];
bf_gain_32TRX = [];
bf_gain_64TRX = [];

% Extract beamforming gains from serving BS-UE pairs
for bs = 1:results_8TRX.numBS
    for ue = 1:results_8TRX.numUEs
        if serving_8TRX(bs, ue)
            bf_gain_8TRX = [bf_gain_8TRX; results_8TRX.allBeamformingGains(bs, ue)];
        end
        
        if bs <= results_32TRX.numBS && ue <= results_32TRX.numUEs && serving_32TRX(bs, ue)
            bf_gain_32TRX = [bf_gain_32TRX; results_32TRX.allBeamformingGains(bs, ue)];
        end
        
        if bs <= results_64TRX.numBS && ue <= results_64TRX.numUEs && serving_64TRX(bs, ue)
            bf_gain_64TRX = [bf_gain_64TRX; results_64TRX.allBeamformingGains(bs, ue)];
        end
    end
end

% Calculate average beamforming gain
avg_bf_gain_8TRX = mean(bf_gain_8TRX);
avg_bf_gain_32TRX = mean(bf_gain_32TRX);
avg_bf_gain_64TRX = mean(bf_gain_64TRX);

% Create bar chart
bar_handle = bar([avg_bf_gain_8TRX, avg_bf_gain_32TRX, avg_bf_gain_64TRX]);
set(bar_handle, 'FaceColor', 'flat');
bar_handle.CData(1,:) = colors_8TRX;
bar_handle.CData(2,:) = colors_32TRX;
bar_handle.CData(3,:) = colors_64TRX;

set(gca, 'XTickLabel', {'8TRX_64AE', '32TRX_128AE', '64TRX_192AE'});
ylabel('Beamforming Gain (dB)');
title('Beamforming Gain Comparison');
grid on;

% Add text with values
text(1, avg_bf_gain_8TRX + 0.2, sprintf('%.2f dB', avg_bf_gain_8TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(2, avg_bf_gain_32TRX + 0.2, sprintf('%.2f dB', avg_bf_gain_32TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(3, avg_bf_gain_64TRX + 0.2, sprintf('%.2f dB', avg_bf_gain_64TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');

% Calculate and show improvements
gain_difference_32vs8 = avg_bf_gain_32TRX - avg_bf_gain_8TRX;
gain_difference_64vs32 = avg_bf_gain_64TRX - avg_bf_gain_32TRX;
gain_difference_64vs8 = avg_bf_gain_64TRX - avg_bf_gain_8TRX;

text(1.5, (avg_bf_gain_8TRX + avg_bf_gain_32TRX)/2, sprintf('%.2f dB', gain_difference_32vs8), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9, ...
    'BackgroundColor', [0.9 0.9 0.9]);
text(2.5, (avg_bf_gain_32TRX + avg_bf_gain_64TRX)/2, sprintf('%.2f dB', gain_difference_64vs32), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9, ...
    'BackgroundColor', [0.9 0.9 0.9]);

%% 2. Signal Power (RSRP) Comparison
subplot(1, 3, 2);

% Extract valid RSRP values for serving pairs only
rsrp_8TRX = [];
rsrp_32TRX = [];
rsrp_64TRX = [];

% Extract RSRP values from serving BS-UE pairs
for bs = 1:results_8TRX.numBS
    for ue = 1:results_8TRX.numUEs
        if serving_8TRX(bs, ue)
            rsrp_8TRX = [rsrp_8TRX; results_8TRX.rxPower(bs, ue)];
        end
        
        if bs <= results_32TRX.numBS && ue <= results_32TRX.numUEs && serving_32TRX(bs, ue)
            rsrp_32TRX = [rsrp_32TRX; results_32TRX.rxPower(bs, ue)];
        end
        
        if bs <= results_64TRX.numBS && ue <= results_64TRX.numUEs && serving_64TRX(bs, ue)
            rsrp_64TRX = [rsrp_64TRX; results_64TRX.rxPower(bs, ue)];
        end
    end
end

% Calculate average RSRP
avg_rsrp_8TRX = mean(rsrp_8TRX);
avg_rsrp_32TRX = mean(rsrp_32TRX);
avg_rsrp_64TRX = mean(rsrp_64TRX);

% Create bar chart for average RSRP
bar_handle = bar([avg_rsrp_8TRX, avg_rsrp_32TRX, avg_rsrp_64TRX]);
set(bar_handle, 'FaceColor', 'flat');
bar_handle.CData(1,:) = colors_8TRX;
bar_handle.CData(2,:) = colors_32TRX;
bar_handle.CData(3,:) = colors_64TRX;

set(gca, 'XTickLabel', {'8TRX_64AE', '32TRX_128AE', '64TRX_192AE'});
ylabel('Average RSRP (dBm)');
title('Received Signal Power Comparison');
grid on;

% Add text with values
text(1, avg_rsrp_8TRX + 1, sprintf('%.2f dBm', avg_rsrp_8TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(2, avg_rsrp_32TRX + 1, sprintf('%.2f dBm', avg_rsrp_32TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
text(3, avg_rsrp_64TRX + 1, sprintf('%.2f dBm', avg_rsrp_64TRX), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');

% Calculate and show improvements
rsrp_improvement_32vs8 = avg_rsrp_32TRX - avg_rsrp_8TRX;
rsrp_improvement_64vs32 = avg_rsrp_64TRX - avg_rsrp_32TRX;
rsrp_improvement_64vs8 = avg_rsrp_64TRX - avg_rsrp_8TRX;

text(1.5, (avg_rsrp_8TRX + avg_rsrp_32TRX)/2, sprintf('%.2f dB', rsrp_improvement_32vs8), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9, ...
    'BackgroundColor', [0.9 0.9 0.9]);
text(2.5, (avg_rsrp_32TRX + avg_rsrp_64TRX)/2, sprintf('%.2f dB', rsrp_improvement_64vs32), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9, ...
    'BackgroundColor', [0.9 0.9 0.9]);

%% 3. Modulation Order Distribution
subplot(1, 3, 3);

% Extract CQI values for serving pairs to determine modulation orders
cqi_8TRX = [];
cqi_32TRX = [];
cqi_64TRX = [];

% Extract CQI values from serving BS-UE pairs
for bs = 1:results_8TRX.numBS
    for ue = 1:results_8TRX.numUEs
        if serving_8TRX(bs, ue)
            cqi_8TRX = [cqi_8TRX; results_8TRX.cqiValues(bs, ue)];
        end
        
        if bs <= results_32TRX.numBS && ue <= results_32TRX.numUEs && serving_32TRX(bs, ue)
            cqi_32TRX = [cqi_32TRX; results_32TRX.cqiValues(bs, ue)];
        end
        
        if bs <= results_64TRX.numBS && ue <= results_64TRX.numUEs && serving_64TRX(bs, ue)
            cqi_64TRX = [cqi_64TRX; results_64TRX.cqiValues(bs, ue)];
        end
    end
end

% Calculate modulation order distribution
mod_orders_8TRX = zeros(1, 4); % QPSK, 16QAM, 64QAM, 256QAM
mod_orders_32TRX = zeros(1, 4);
mod_orders_64TRX = zeros(1, 4);

% Count modulation orders based on CQI values
for cqi = cqi_8TRX'
    if cqi <= 6
        mod_orders_8TRX(1) = mod_orders_8TRX(1) + 1; % QPSK
    elseif cqi <= 9
        mod_orders_8TRX(2) = mod_orders_8TRX(2) + 1; % 16QAM
    elseif cqi <= 14
        mod_orders_8TRX(3) = mod_orders_8TRX(3) + 1; % 64QAM
    else
        mod_orders_8TRX(4) = mod_orders_8TRX(4) + 1; % 256QAM
    end
end

for cqi = cqi_32TRX'
    if cqi <= 6
        mod_orders_32TRX(1) = mod_orders_32TRX(1) + 1; % QPSK
    elseif cqi <= 9
        mod_orders_32TRX(2) = mod_orders_32TRX(2) + 1; % 16QAM
    elseif cqi <= 14
        mod_orders_32TRX(3) = mod_orders_32TRX(3) + 1; % 64QAM
    else
        mod_orders_32TRX(4) = mod_orders_32TRX(4) + 1; % 256QAM
    end
end

for cqi = cqi_64TRX'
    if cqi <= 6
        mod_orders_64TRX(1) = mod_orders_64TRX(1) + 1; % QPSK
    elseif cqi <= 9
        mod_orders_64TRX(2) = mod_orders_64TRX(2) + 1; % 16QAM
    elseif cqi <= 14
        mod_orders_64TRX(3) = mod_orders_64TRX(3) + 1; % 64QAM
    else
        mod_orders_64TRX(4) = mod_orders_64TRX(4) + 1; % 256QAM
    end
end

% Convert to percentages
if sum(mod_orders_8TRX) > 0
    mod_orders_8TRX = (mod_orders_8TRX / sum(mod_orders_8TRX)) * 100;
end

if sum(mod_orders_32TRX) > 0
    mod_orders_32TRX = (mod_orders_32TRX / sum(mod_orders_32TRX)) * 100;
end

if sum(mod_orders_64TRX) > 0
    mod_orders_64TRX = (mod_orders_64TRX / sum(mod_orders_64TRX)) * 100;
end

% Create grouped bar chart
mod_data = [mod_orders_8TRX; mod_orders_32TRX; mod_orders_64TRX]';
b = bar(mod_data, 'grouped');
set(b(1), 'FaceColor', colors_8TRX);
set(b(2), 'FaceColor', colors_32TRX);
set(b(3), 'FaceColor', colors_64TRX);

% Add labels
modulation_labels = {'QPSK', '16QAM', '64QAM', '256QAM'};
set(gca, 'XTickLabel', modulation_labels);
ylabel('Percentage (%)');
title('Modulation Order Distribution');
legend('8TRX_64AE', '32TRX_128AE', '64TRX_192AE');
grid on;

% Add percentage values on top of each bar (only for values > 5% to avoid clutter)
for i = 1:length(mod_orders_8TRX)
    if mod_orders_8TRX(i) > 5
        text(i-0.25, mod_orders_8TRX(i)+1, [num2str(mod_orders_8TRX(i), '%.0f') '%'], ...
            'HorizontalAlignment', 'center', 'FontSize', 7);
    end
    if mod_orders_32TRX(i) > 5
        text(i, mod_orders_32TRX(i)+1, [num2str(mod_orders_32TRX(i), '%.0f') '%'], ...
            'HorizontalAlignment', 'center', 'FontSize', 7);
    end
    if mod_orders_64TRX(i) > 5
        text(i+0.25, mod_orders_64TRX(i)+1, [num2str(mod_orders_64TRX(i), '%.0f') '%'], ...
            'HorizontalAlignment', 'center', 'FontSize', 7);
    end
end

sgtitle('Signal Power and Modulation Analysis (Serving Pairs Only)', 'FontSize', 14, 'FontWeight', 'bold');

%% Create a third figure to examine connection quality and CQI distribution
figure('Name', 'Connection Quality Analysis', 'Position', [100, 100, 1200, 400]);

%% 1. Valid Connections Comparison
subplot(1, 2, 1);

% Instead of counting percentages, let's show the distribution of UEs across BSs
% Count how many UEs are served by each BS
bs_ue_counts_8TRX = zeros(results_8TRX.numBS, 1);
bs_ue_counts_32TRX = zeros(results_32TRX.numBS, 1);
bs_ue_counts_64TRX = zeros(results_64TRX.numBS, 1);

% For each UE, find which BS serves it
for ue = 1:results_8TRX.numUEs
    for bs = 1:results_8TRX.numBS
        if serving_8TRX(bs, ue)
            bs_ue_counts_8TRX(bs) = bs_ue_counts_8TRX(bs) + 1;
            break; % Each UE should only have one serving BS
        end
    end
    
    if ue <= results_32TRX.numUEs
        for bs = 1:results_32TRX.numBS
            if serving_32TRX(bs, ue)
                bs_ue_counts_32TRX(bs) = bs_ue_counts_32TRX(bs) + 1;
                break;
            end
        end
    end
    
    if ue <= results_64TRX.numUEs
        for bs = 1:results_64TRX.numBS
            if serving_64TRX(bs, ue)
                bs_ue_counts_64TRX(bs) = bs_ue_counts_64TRX(bs) + 1;
                break;
            end
        end
    end
end

% Create a grouped bar chart for BS distribution
bs_data = [bs_ue_counts_8TRX, bs_ue_counts_32TRX, bs_ue_counts_64TRX];
b = bar(bs_data, 'grouped');
set(b(1), 'FaceColor', colors_8TRX);
set(b(2), 'FaceColor', colors_32TRX);
set(b(3), 'FaceColor', colors_64TRX);

xlabel('Base Station ID');
ylabel('Number of UEs Served');
title('UE Distribution Across Base Stations');
legend('8TRX_64AE', '32TRX_128AE', '64TRX_192AE', 'Location', 'best');
grid on;
xticks(1:max([results_8TRX.numBS, results_32TRX.numBS, results_64TRX.numBS]));

% Add text with values on top of each bar
for bs = 1:results_8TRX.numBS
    text(bs-0.25, bs_ue_counts_8TRX(bs)+0.5, num2str(bs_ue_counts_8TRX(bs)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
    
    if bs <= results_32TRX.numBS
        text(bs, bs_ue_counts_32TRX(bs)+0.5, num2str(bs_ue_counts_32TRX(bs)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    if bs <= results_64TRX.numBS
        text(bs+0.25, bs_ue_counts_64TRX(bs)+0.5, num2str(bs_ue_counts_64TRX(bs)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end

%% 2. CQI Distribution
subplot(1, 2, 2);

% Calculate CQI distribution
cqi_bins = 1:16;
cqi_dist_8TRX = histcounts(cqi_8TRX, [cqi_bins, 17]) / length(cqi_8TRX) * 100;
cqi_dist_32TRX = histcounts(cqi_32TRX, [cqi_bins, 17]) / length(cqi_32TRX) * 100;
cqi_dist_64TRX = histcounts(cqi_64TRX, [cqi_bins, 17]) / length(cqi_64TRX) * 100;

% Calculate average CQI
avg_cqi_8TRX = mean(cqi_8TRX);
avg_cqi_32TRX = mean(cqi_32TRX);
avg_cqi_64TRX = mean(cqi_64TRX);

% Plot CQI distribution
plot(cqi_bins, cqi_dist_8TRX, 'Color', colors_8TRX, 'LineWidth', 2, 'Marker', 'o');
hold on;
plot(cqi_bins, cqi_dist_32TRX, 'Color', colors_32TRX, 'LineWidth', 2, 'Marker', 's');
plot(cqi_bins, cqi_dist_64TRX, 'Color', colors_64TRX, 'LineWidth', 2, 'Marker', 'd');

grid on;
xlabel('CQI Value');
ylabel('Percentage (%)');
title('CQI Distribution (Serving Pairs)');
legend(['8TRX (Avg: ' num2str(avg_cqi_8TRX, '%.1f') ')'], ...
       ['32TRX (Avg: ' num2str(avg_cqi_32TRX, '%.1f') ')'], ...
       ['64TRX (Avg: ' num2str(avg_cqi_64TRX, '%.1f') ')'], ...
       'Location', 'best');

% Add vertical lines for average CQI values
xline(avg_cqi_8TRX, 'Color', colors_8TRX, 'LineStyle', '--');
xline(avg_cqi_32TRX, 'Color', colors_32TRX, 'LineStyle', '--');
xline(avg_cqi_64TRX, 'Color', colors_64TRX, 'LineStyle', '--');

sgtitle('Connection Quality Analysis (Serving Pairs Only)', 'FontSize', 14, 'FontWeight', 'bold');

%% Summary Box for All Configurations
% Define the peak improvement variable
peak_improvement_32vs8 = improvement_32vs8;
peak_improvement_64vs32 = improvement_64vs32;
peak_improvement_64vs8 = improvement_64vs8;

% Print summary to the command window for all configuration comparisons
fprintf('\n----- MIMO Configuration Comparison Summary -----\n');

% 32TRX vs 8TRX comparison
fprintf('32TRX_128AE vs 8TRX_64AE Comparison:\n');
fprintf('- Beamforming gain improvement: %.2f dB\n', gain_difference_32vs8);
fprintf('- Signal power improvement: %.2f dB\n', rsrp_improvement_32vs8);
fprintf('- Average throughput improvement: %.1f%%\n', avg_improvement_32vs8);
fprintf('- Peak throughput improvement: %.1f%%\n', improvement_32vs8);
fprintf('\n');

% 64TRX vs 32TRX comparison
fprintf('64TRX_192AE vs 32TRX_128AE Comparison:\n');
fprintf('- Beamforming gain improvement: %.2f dB\n', gain_difference_64vs32);
fprintf('- Signal power improvement: %.2f dB\n', rsrp_improvement_64vs32);
fprintf('- Average throughput improvement: %.1f%%\n', avg_improvement_64vs32);
fprintf('- Peak throughput improvement: %.1f%%\n', improvement_64vs32);
fprintf('\n');

% 64TRX vs 8TRX comparison (overall improvement)
fprintf('64TRX_192AE vs 8TRX_64AE Comparison (Overall):\n');
fprintf('- Beamforming gain improvement: %.2f dB\n', gain_difference_64vs8);
fprintf('- Signal power improvement: %.2f dB\n', rsrp_improvement_64vs8);
fprintf('- Average throughput improvement: %.1f%%\n', avg_improvement_64vs8);
fprintf('- Peak throughput improvement: %.1f%%\n', peak_improvement_64vs8);
fprintf('\n');

fprintf('--------------------------------------------------\n\n');

%% Create a simplified figure for MIMO spatial multiplexing gains analysis
figure('Name', 'MIMO Spatial Multiplexing Gains', 'Position', [100, 100, 900, 400]);

%% Check if MIMO data is available
has_mimo_data = isfield(results_8TRX, 'allNumLayers') && ...
                isfield(results_32TRX, 'allNumLayers') && ...
                isfield(results_64TRX, 'allNumLayers');

if ~has_mimo_data
    % Display message if no MIMO data available
    text(0.5, 0.5, 'MIMO layer data not available in results', ...
        'HorizontalAlignment', 'center');
    axis off;
    return;
end

%% 1. Distribution of number of layers used
subplot(1, 2, 1);

% Extract layer information for each configuration
layers_8TRX = [];
layers_32TRX = [];
layers_64TRX = [];

% Extract layers from serving BS-UE pairs
for bs = 1:results_8TRX.numBS
    for ue = 1:results_8TRX.numUEs
        if serving_8TRX(bs, ue) && ~isempty(results_8TRX.allNumLayers{bs, ue})
            layers_8TRX = [layers_8TRX; results_8TRX.allNumLayers{bs, ue}];
        end
    end
end

for bs = 1:results_32TRX.numBS
    for ue = 1:results_32TRX.numUEs
        if serving_32TRX(bs, ue) && ~isempty(results_32TRX.allNumLayers{bs, ue})
            layers_32TRX = [layers_32TRX; results_32TRX.allNumLayers{bs, ue}];
        end
    end
end

for bs = 1:results_64TRX.numBS
    for ue = 1:results_64TRX.numUEs
        if serving_64TRX(bs, ue) && ~isempty(results_64TRX.allNumLayers{bs, ue})
            layers_64TRX = [layers_64TRX; results_64TRX.allNumLayers{bs, ue}];
        end
    end
end

% Calculate the maximum layer count across all configurations
max_layers = max([max(layers_8TRX), max(layers_32TRX), max(layers_64TRX)]);
if isempty(max_layers) || max_layers < 1
    max_layers = 1;
end

% Calculate distribution of layers for each configuration
layer_counts_8TRX = histcounts(layers_8TRX, 0.5:1:(max_layers+0.5));
layer_counts_32TRX = histcounts(layers_32TRX, 0.5:1:(max_layers+0.5));
layer_counts_64TRX = histcounts(layers_64TRX, 0.5:1:(max_layers+0.5));

% Convert to percentages
if ~isempty(layer_counts_8TRX) && sum(layer_counts_8TRX) > 0
    layer_percent_8TRX = (layer_counts_8TRX / sum(layer_counts_8TRX)) * 100;
else
    layer_percent_8TRX = zeros(1, max_layers);
end

if ~isempty(layer_counts_32TRX) && sum(layer_counts_32TRX) > 0
    layer_percent_32TRX = (layer_counts_32TRX / sum(layer_counts_32TRX)) * 100;
else
    layer_percent_32TRX = zeros(1, max_layers);
end

if ~isempty(layer_counts_64TRX) && sum(layer_counts_64TRX) > 0
    layer_percent_64TRX = (layer_counts_64TRX / sum(layer_counts_64TRX)) * 100;
else
    layer_percent_64TRX = zeros(1, max_layers);
end

% Plot distribution as individual bars
bar_width = 0.25;
hold on;

% Create one bar of each color for legend purposes
h1 = bar(NaN, NaN, bar_width, 'FaceColor', colors_8TRX, 'EdgeColor', 'k');
h2 = bar(NaN, NaN, bar_width, 'FaceColor', colors_32TRX, 'EdgeColor', 'k');
h3 = bar(NaN, NaN, bar_width, 'FaceColor', colors_64TRX, 'EdgeColor', 'k');

% Plot actual data
for i = 1:max_layers
    % Plot bars with slight offset for each configuration
    if i <= length(layer_percent_8TRX) && layer_percent_8TRX(i) > 0
        bar(i-bar_width, layer_percent_8TRX(i), bar_width, 'FaceColor', colors_8TRX, 'EdgeColor', 'k');
        % Add text label for values
        if layer_percent_8TRX(i) > 5
            text(i-bar_width, layer_percent_8TRX(i)+2, sprintf('%.0f%%', layer_percent_8TRX(i)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
    
    if i <= length(layer_percent_32TRX) && layer_percent_32TRX(i) > 0
        bar(i, layer_percent_32TRX(i), bar_width, 'FaceColor', colors_32TRX, 'EdgeColor', 'k');
        % Add text label for values
        if layer_percent_32TRX(i) > 5
            text(i, layer_percent_32TRX(i)+2, sprintf('%.0f%%', layer_percent_32TRX(i)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
    
    if i <= length(layer_percent_64TRX) && layer_percent_64TRX(i) > 0
        bar(i+bar_width, layer_percent_64TRX(i), bar_width, 'FaceColor', colors_64TRX, 'EdgeColor', 'k');
        % Add text label for values
        if layer_percent_64TRX(i) > 5
            text(i+bar_width, layer_percent_64TRX(i)+2, sprintf('%.0f%%', layer_percent_64TRX(i)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
end
hold off;

% Add labels
xlabel('Number of Spatial Layers (Rank)');
ylabel('Percentage of UEs (%)');
title('Distribution of Spatial Layers Used');
legend([h1, h2, h3], {'8TRX_64AE', '32TRX_128AE', '64TRX_192AE'}, 'Location', 'best');
grid on;
xticks(1:max_layers);
xlim([0.5, max_layers+0.5]);

%% 2. Average Throughput by Number of Layers
subplot(1, 2, 2);

% Calculate average throughput per number of layers for each configuration
avg_tput_by_layer_8TRX = zeros(1, max_layers);
avg_tput_by_layer_32TRX = zeros(1, max_layers);
avg_tput_by_layer_64TRX = zeros(1, max_layers);

% Fill in zeros for counts to avoid division by zero
count_by_layer_8TRX = zeros(1, max_layers);
count_by_layer_32TRX = zeros(1, max_layers);
count_by_layer_64TRX = zeros(1, max_layers);

% For 8TRX
for bs = 1:results_8TRX.numBS
    for ue = 1:results_8TRX.numUEs
        if serving_8TRX(bs, ue) && ~isempty(results_8TRX.allNumLayers{bs, ue})
            layer = results_8TRX.allNumLayers{bs, ue};
            if layer <= max_layers
                avg_tput_by_layer_8TRX(layer) = avg_tput_by_layer_8TRX(layer) + results_8TRX.effectiveThroughput(bs, ue);
                count_by_layer_8TRX(layer) = count_by_layer_8TRX(layer) + 1;
            end
        end
    end
end

% For 32TRX
for bs = 1:results_32TRX.numBS
    for ue = 1:results_32TRX.numUEs
        if serving_32TRX(bs, ue) && ~isempty(results_32TRX.allNumLayers{bs, ue})
            layer = results_32TRX.allNumLayers{bs, ue};
            if layer <= max_layers
                avg_tput_by_layer_32TRX(layer) = avg_tput_by_layer_32TRX(layer) + results_32TRX.effectiveThroughput(bs, ue);
                count_by_layer_32TRX(layer) = count_by_layer_32TRX(layer) + 1;
            end
        end
    end
end

% For 64TRX
for bs = 1:results_64TRX.numBS
    for ue = 1:results_64TRX.numUEs
        if serving_64TRX(bs, ue) && ~isempty(results_64TRX.allNumLayers{bs, ue})
            layer = results_64TRX.allNumLayers{bs, ue};
            if layer <= max_layers
                avg_tput_by_layer_64TRX(layer) = avg_tput_by_layer_64TRX(layer) + results_64TRX.effectiveThroughput(bs, ue);
                count_by_layer_64TRX(layer) = count_by_layer_64TRX(layer) + 1;
            end
        end
    end
end

% Calculate averages, handling division by zero
for layer = 1:max_layers
    if count_by_layer_8TRX(layer) > 0
        avg_tput_by_layer_8TRX(layer) = avg_tput_by_layer_8TRX(layer) / count_by_layer_8TRX(layer);
    else
        avg_tput_by_layer_8TRX(layer) = 0;
    end
    
    if count_by_layer_32TRX(layer) > 0
        avg_tput_by_layer_32TRX(layer) = avg_tput_by_layer_32TRX(layer) / count_by_layer_32TRX(layer);
    else
        avg_tput_by_layer_32TRX(layer) = 0;
    end
    
    if count_by_layer_64TRX(layer) > 0
        avg_tput_by_layer_64TRX(layer) = avg_tput_by_layer_64TRX(layer) / count_by_layer_64TRX(layer);
    else
        avg_tput_by_layer_64TRX(layer) = 0;
    end
end

% Plot with a simplified legend approach
hold on;
% Create bars for legend only
h1 = bar(NaN, NaN, bar_width, 'FaceColor', colors_8TRX, 'EdgeColor', 'k');
h2 = bar(NaN, NaN, bar_width, 'FaceColor', colors_32TRX, 'EdgeColor', 'k');
h3 = bar(NaN, NaN, bar_width, 'FaceColor', colors_64TRX, 'EdgeColor', 'k');

% Plot actual data
for i = 1:max_layers
    % Plot bars with slight offset for each configuration
    if avg_tput_by_layer_8TRX(i) > 0
        bar(i-bar_width, avg_tput_by_layer_8TRX(i), bar_width, 'FaceColor', colors_8TRX, 'EdgeColor', 'k');
        
        % Add text with value above bar
        text(i-bar_width, avg_tput_by_layer_8TRX(i)+5, num2str(avg_tput_by_layer_8TRX(i), '%.0f'), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    if avg_tput_by_layer_32TRX(i) > 0
        bar(i, avg_tput_by_layer_32TRX(i), bar_width, 'FaceColor', colors_32TRX, 'EdgeColor', 'k');
        
        % Add text with value above bar
        text(i, avg_tput_by_layer_32TRX(i)+5, num2str(avg_tput_by_layer_32TRX(i), '%.0f'), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    if avg_tput_by_layer_64TRX(i) > 0
        bar(i+bar_width, avg_tput_by_layer_64TRX(i), bar_width, 'FaceColor', colors_64TRX, 'EdgeColor', 'k');
        
        % Add text with value above bar
        text(i+bar_width, avg_tput_by_layer_64TRX(i)+5, num2str(avg_tput_by_layer_64TRX(i), '%.0f'), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
hold off;

% Add labels
xlabel('Number of Spatial Layers (Rank)');
ylabel('Average Throughput (Mbps)');
title('Average Throughput by Number of Spatial Layers');
legend([h1, h2, h3], {'8TRX_64AE', '32TRX_128AE', '64TRX_192AE'}, 'Location', 'best');
grid on;
xticks(1:max_layers);
xlim([0.5, max_layers+0.5]);

% Add overall title
sgtitle('MIMO Spatial Multiplexing Performance Analysis', 'FontSize', 14, 'FontWeight', 'bold');

% Add overall title
sgtitle('MIMO Spatial Multiplexing Performance Analysis', 'FontSize', 14, 'FontWeight', 'bold');
