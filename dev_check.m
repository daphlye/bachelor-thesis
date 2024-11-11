function [res] = dev_check(inputArr,max)
%STDDEV_PERC tests if all values in the INPUTARR are within a range of MAX
%around the mean (check the units of your data!)
%   INPUT: inputArr - array you want to test 
%          max - Maximum difference between the value in your array and the
%          mean of all values in that array
%
%   OUTPUT: res - returns 1 if all values are within that range, 0 if not. 
%
    devVal = abs(mean(inputArr)-inputArr);
    res = ~any(devVal(~isnan(devVal)) > max); % returns 0 if at least one element is larger than max and ignores NaN-vals
end

