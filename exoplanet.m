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
starID = extractBetween(tline, '"', '"'); % Extract the star ID using extractBetween
fprintf('Star ID: %s\n', starID{1}); % Access the first element of the cell array

% Read in the data from file line by line
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

% Plot the data with error bars
figure;
errorbar(RadialVelocityData(:, 1), RadialVelocityData(:, 2), RadialVelocityData(:, 3), 'o-', 'Color', 'b', 'MarkerFaceColor', 'r');
hold on;
errorbar(RadialVelocityData(:, 1), RadialVelocityData(:, 2), RadialVelocityData(:, 3), 'o', 'Color', 'c');
hold off;
xlabel('Days since Start');
ylabel('Radial Velocity');
title([starID, '- Radial Velocity vs Time']);
grid on;

fclose(fileID);

% Exoplanet Detection via Lomb-Scargle Periodogram with Peak Detection
time = RadialVelocityData(:, 1);
velocity = RadialVelocityData(:, 2);
sigma = RadialVelocityData(:, 3);

% Normalize Time (makes sure time starts at zero for numerical stability)
time = time - min(time);

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
minPeriod = 1;    % days
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
threshold = 0.2 * max(smoothedPower);  % Only consider peaks above 20% of the strongest signal
[peakVals, peakLocs] = findpeaks(smoothedPower, filteredPeriods, ...
                                 'MinPeakHeight', threshold, ...
                                 'MinPeakDistance', 5); % Detects logical maxima in smoothed power spectrum

% Annotate peaks
set(gca, 'XScale', 'log');
hold on;
plot(peakLocs, peakVals, 'ro', 'MarkerFaceColor', 'r');  % Red dots at peak locations
text(peakLocs + 0.5, peakVals, compose('%.1f d', peakLocs), ...
     'Color', 'red', 'FontSize', 8);
hold off;

% Print detected periods
fprintf('Detected candidate orbital periods:\n');
for i = 1:length(peakLocs)
    fprintf('  %.2f days (Power = %.3f)\n', peakLocs(i), peakVals(i));
end