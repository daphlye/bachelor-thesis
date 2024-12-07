function [structure] = jr_data_ana(jr_data, R_gsm_data)
%JR_DATA_ANA Analyze if inputted current points in the right direction.
%   Remember to change the name of the structure you want to put this data
%   in sufficiently, so that you dont overwrite!
%
%   also please remember to INPUT jr blablabla .DATA!!! it wont work with
%   the TSeries as an input!
%
%    Output is going to be:
%   STRUCTURE(:,1)=MEAN
%   STRUCTURE(:,2/3)=MEAN_a0/b0
%   STRUCTURE(:,4/5)=STD_a0/b0
%   STRUCTURE(:,6)=PERCENT
%   
%   ======================================================================
structure(:,1)= mean(jr_data, 'omitnan');

structure(:,2) = mean(jr_data(jr_data>0)); % data above 0
structure(:,3)= abs(mean(jr_data(jr_data<0))); % data below 0
structure(:,4) = std(jr_data(jr_data>0));
structure(:,5) = std(jr_data(jr_data<0));
% calculate percentage of peaks that point in right direction:
if mean(R_gsm_data(:,2))> 0 % dawnside +y, J*r should point in +r, away from earth
    structure(:,6) = numel(jr_data(jr_data>0))./numel(jr_data);    % percentage in right direction
else
    structure(:,6) = (numel(jr_data(jr_data<0))./numel(jr_data))*100; % J should be negative, 
end

end