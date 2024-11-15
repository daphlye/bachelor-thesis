function [xRange,yRange] = best_legend_pos(timeSeries)
%BEST_LEGEND_POS finds the best position to put the legend in in a plot.
%   works as an input for the irf_legend legend_position argument. 
%   NOT YET WORKING

% Data-based dynamic legend placement
data = timeSeries.data;
time = timeSeries.time;

% Find location with minimum data density
xRange = [time(1) time(end)];
yRange = [min(data), max(data)];
legendPosition = [0.1, 0.9]; % default

if mean(data) < (max(data) - min(data))/2 % Adjust to place it on less crowded side
    legendPosition = [0.88, 0.10];
end

end