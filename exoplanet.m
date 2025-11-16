% Sean Hordines, Lee McKinstry, Name 3
% Digital Information Processing
% Yong Wei

filename = 'data/Proxima_Cen_(1).tbl';
fileID = fopen(filename, 'r');

if fileID == -1
    error('File could not be opened. Check if the file exists and the name is correct.');
end

RadialVelocityData = [];

% Read in the star ID
tline = fgetl(fileID);
starID = extractBetween(tline, '"', '"');
fprintf('Star ID: %s\n', starID{1});

% Read in the data
while ischar(tline)
    if ~startsWith(tline, '|') && ~startsWith(tline, '\')
        dataPoints = sscanf(tline, '%f', 3);
        RadialVelocityData = [RadialVelocityData; dataPoints'];
    end
    tline = fgetl(fileID);
end
fclose(fileID);

% Plot the data
figure;
errorbar(RadialVelocityData(:, 1), RadialVelocityData(:, 2), RadialVelocityData(:, 3), 'o-', 'Color', 'b', 'MarkerFaceColor', 'r');
xlabel('Time (days)');
ylabel('Radial Velocity');
title([starID, ' - Radial Velocity vs Time']);
grid on;

% Extract data
time = RadialVelocityData(:, 1);
velocity = RadialVelocityData(:, 2);
sigma = RadialVelocityData(:, 3);

% Normalize time
time = time - min(time);

% Signal-to-noise ratio
signalPower = mean(velocity.^2);
noisePower = mean(sigma.^2);
SNR = 10 * log10(signalPower / noisePower);
fprintf('Signal to Noise Ratio (SNR): %.2f dB\n', SNR);

% Frequency range
minFreq = 1 / max(time);
maxFreq = 1;
freq = linspace(minFreq, maxFreq, 10000);

% Compute weighted Lomb-Scargle periodogram
[power, freq] = weighted_lomb_scargle(time, velocity, sigma, freq);

% Convert to periods and reverse for plotting
periods = 1 ./ freq;
power = flip(power);
periods = flip(periods);

% Bandpass filter
minPeriod = 1;
maxPeriod = 250;
bandpassIdx = periods >= minPeriod & periods <= maxPeriod;
filteredPeriods = periods(bandpassIdx);
filteredPower = power(bandpassIdx);

% Plot periodogram
figure;
plot(filteredPeriods, filteredPower, 'b');
xlabel('Period (days)');
ylabel('Weighted Lomb-Scargle Power');
title([starID, ' - Weighted LS Periodogram']);
grid on;

% Peak detection
smoothedPower = smooth(filteredPower, 50);
threshold = 0.2 * max(smoothedPower);
[peakVals, peakLocs] = findpeaks(smoothedPower, filteredPeriods, ...
                                 'MinPeakHeight', threshold, ...
                                 'MinPeakDistance', 5);

% Plot the weighted Lomb-Scargle periodogram with peaks
figure;
plot(filteredPeriods, filteredPower, 'b-', 'LineWidth', 1.5);
xlabel('Period (days)');
ylabel('Weighted Lomb-Scargle Power');
title([starID, ' - Weighted LS Periodogram']);
set(gca, 'XScale', 'log');
grid on;
hold on;

% Overlay detected peaks
plot(peakLocs, peakVals, 'ro', 'MarkerFaceColor', 'r');

% Annotate each peak with its period
for i = 1:length(peakLocs)
    text(peakLocs(i) + 0.5, peakVals(i), sprintf('%.1f d', peakLocs(i)), ...
         'Color', 'red', 'FontSize', 8);
end

hold off;

% --- Weighted Lomb-Scargle Periodogram Function ---
function [power, freq] = weighted_lomb_scargle(time, velocity, sigma, freq)
    % Ensure column vectors
    time = time(:);
    velocity = velocity(:);
    sigma = sigma(:);
    freq = freq(:);

    % Compute weights from error bars
    w = 1 ./ sigma.^2;
    W = sum(w);

    % Weighted mean of the velocity
    y_mean = sum(w .* velocity) / W;

    % Center the data
    y = velocity - y_mean;

    % Initialize power array
    power = zeros(size(freq));

    % Loop over each frequency
    for i = 1:length(freq)
        omega = 2 * pi * freq(i);

        % Compute phase offset tau
        tan_2omega_tau = sum(w .* sin(2 * omega * time)) / sum(w .* cos(2 * omega * time));
        tau = atan(tan_2omega_tau) / (2 * omega);

        % Compute sine and cosine terms
        cos_term = cos(omega * (time - tau));
        sin_term = sin(omega * (time - tau));

        % Weighted projections
        C = sum(w .* y .* cos_term);
        S = sum(w .* y .* sin_term);
        CC = sum(w .* cos_term.^2);
        SS = sum(w .* sin_term.^2);

        % Compute normalized power
        power(i) = (C^2 / CC + S^2 / SS) / (2 * sum(w .* y.^2));
    end
end