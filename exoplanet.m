% Sean Hordines, Lee McKinstry, Name 3
% Digital Information Processing
% Yong Wei

filename = 'data/Proxima_Cen_(1).tbl';
fileID = fopen(filename, 'r');

if fileID == -1
    error('Files could not be opened. Check if the file exists and the name is correct.');
end

% Read in the star ID
tline = fgetl(fileID);
starID = extractBetween(tline, '"', '"'); % Extract the star ID using extractBetween
fprintf('Star ID: %s\n', starID{1}); % Access the first element of the cell array

% Read in the data from file
RadialVelocityData = [];
while ischar(tline)
    if ~startsWith(tline, '|') && ~startsWith(tline, '\')
        dataPoints = sscanf(tline, '%f', 3);
        RadialVelocityData = [RadialVelocityData; dataPoints'];
    end
    tline = fgetl(fileID);
end

% Preview the data
format long g;
% disp(size(RadialVelocityData));
% disp(RadialVelocityData);

% Exoplanet Detection via Lomb-Scargle Periodogram with Peak Detection
time = RadialVelocityData(:, 1);
velocity = RadialVelocityData(:, 2);
sigma = RadialVelocityData(:, 3);

% Normalize Time (makes sure time starts at zero for numerical stability)
time = time - min(time);

% Calculate the signal to noise ratio (SNR)
signalPower = mean(velocity.^2);
noisePower = mean(sigma.^2);
SNR = 10 * log10(signalPower / noisePower);
fprintf('Signal to Noise Ratio (SNR): %.2f dB\n', SNR);

% Plot the data with error bars
figure;
errorbar(time, velocity, sigma, 'o-', 'Color', 'b', 'MarkerFaceColor', 'r');
hold on;
errorbar(time, velocity, sigma, 'o', 'Color', 'c');
hold off;
xlabel('Time (days)');
ylabel('Radial Velocity');
title([starID, '- Radial Velocity vs Time']);
grid on;

% Define Frequency Range (Sets up Frequency sweep to detect periodic signals)
minFreq = 1 / max(time);
maxFreq = 1;
freq = linspace(minFreq, maxFreq, 10000); % create 10,000 evenly spaced freq samples

% Compute Lomb-Scargle periodogram
[power, ~] = plomb(velocity, time, freq, 'normalized');

% Convert freq values to orbital periods (in days)
periods = 1 ./ freq;

% reverse for peak detection
power = flip(power);
periods = flip(periods);

% Apply bandpass filter
minPeriod = 1.15;    % days
maxPeriod = 250;  % days
bandpassIdx = periods >= minPeriod & periods <= maxPeriod;
filteredPeriods = periods(bandpassIdx);
filteredPower = power(bandpassIdx);

% Plot periodogram
figure;
plot(filteredPeriods, filteredPower, 'b');      % Plot power vs. period in blue
xlabel('Period (days)');
ylabel('Lomb-Scargle Power');
title([starID, ' - LS Periodogram']);
grid on;

% Peak Detection
smoothedPower = smooth(filteredPower, 50);  % Optional smoothing with moving average (window size = 50)
threshold = 0.5 * max(filteredPower);  % Only consider peaks above 20% of the strongest signal
[peakVals, peakLocs] = findpeaks(filteredPower, filteredPeriods, ...
                                 'MinPeakHeight', threshold, ...
                                 'MinPeakDistance', 50); % Detects logical maxima in smoothed power spectrum

% Annotate peaks
set(gca, 'XScale', 'log');
hold on;
plot(peakLocs, peakVals, 'ro', 'MarkerFaceColor', 'r');  % Red dots at peak locations
text(peakLocs + 0.5, peakVals, compose('%.1f d', peakLocs), ...
     'Color', 'red', 'FontSize', 8);
hold off;

% Add vertical lines for periods of 11.2 days and its multiples (0.5, 2, 3, and 4)
hold on;
multiples = [1, 2, 3, 4, 5];
truePeriod = 11.2;
for i = multiples
    % xline(i * truePeriod, 'c--', 'LineWidth', 1); % Dashed lines for specified multiples
end
hold off;

% Print detected periods
fprintf('Detected candidate orbital periods:\n');
for i = 1:length(peakLocs)
    fprintf('  %.2f days (Power = %.3f)\n', peakLocs(i), peakVals(i));
end

% Empirical False Alarm Probability via permutation
rng(314159265);
numPerms = 5000;

% Map each detected period to its index in filteredPeriods
peakIdx = zeros(size(peakLocs));
for i = 1:length(peakLocs)
    [~, peakIdx(i)] = min(abs(filteredPeriods - peakLocs(i)));
end

% Preallocate counters
exceedCount = zeros(size(peakIdx));

fprintf('\nRunning %d permutations for empirical FAP (may take time)...\n', numPerms);
maxPowers = zeros(numPerms,1);

for b = 1:numPerms
    % Permute velocity values
    permVel = velocity(randperm(length(velocity)));
    % recompute LS power on same frequency grid (normalized)
    p_boot = plomb(permVel, time, 1./filteredPeriods, 'normalized');
    maxPowers(b) = max(p_boot);
end

% Compute FAP per peak
FAP = zeros(length(peakIdx),1);
for i = 1:length(peakIdx)
    z_obs = filteredPower(peakIdx(i));
    exceedCount(i) = sum(maxPowers >= z_obs);
    FAP(i) = (exceedCount(i) + 1) / (numPerms + 1);
    fprintf('Period %.2f: FAP ≈ %.4f ( %d / %d )\n', ...
            peakLocs(i), FAP(i), exceedCount(i), numPerms);
end