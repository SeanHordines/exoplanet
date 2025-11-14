% Sean Hordines, Name 2, Name 3
% Digital Information Processing
% Yong Wei

filename = 'data/Proxima_Cen_(1).tbl';
fileID = fopen(filename, 'r');

if fileID == -1
    error('File could not be opened. Check if the file exists and the name is correct.');
end

RadialVelocityData = [];

% Read in the data from file line by line
tline = fgetl(fileID);
while ischar(tline)
    if ~startsWith(tline, '|') && ~startsWith(tline, '\')
        dataPoints = sscanf(tline, '%f', 3);
        RadialVelocityData = [RadialVelocityData; dataPoints'];
    end
    tline = fgetl(fileID);
end

% Convert Julian Date to days since start
JulianDate = RadialVelocityData(:, 1);
Days = JulianDate - RadialVelocityData(1, 1);
RadialVelocityData(:, 1) = Days;

% Preview the data
disp(size(RadialVelocityData));
format long g;
disp(RadialVelocityData);

% Plot the data with error bars
figure;
errorbar(RadialVelocityData(:, 1), RadialVelocityData(:, 2), RadialVelocityData(:, 3), 'o-', 'Color', 'b', 'MarkerFaceColor', 'c');
hold on;
errorbar(RadialVelocityData(:, 1), RadialVelocityData(:, 2), RadialVelocityData(:, 3), 'o', 'Color', 'r');
hold off;
xlabel('Days since Start');
ylabel('Radial Velocity');
title('Proxima Centauri'); % Change to match file
grid on;

fclose(fileID);