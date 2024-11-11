function [avg_value, med_value, std_value] = get_stats(afData)
% GET_STATS Return mean, median, and std of inputted data array
%   Turns out current arrays contain nan-values
    avg_value = mean(afData, "omitnan"); %Computes average of inputted data array
    med_value = median(afData, "omitnan"); %Computes median of inputteed data array
    std_value = std(afData, "omitnan");
end
