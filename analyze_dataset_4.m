% make tables 
tmp.filename = fullfile('events/tables');
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end

% Put information from arrays into readable time format
tmp.adjusted_time = good_time(:,1) + (15 * 60);
tmp.adjusted_time = datetime(tmp.adjusted_time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
tmp.end_time = datetime(good_time(:,2), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
values.t_start = tmp.adjusted_time; % Assign datetime array
values.t_end = tmp.end_time;        % Assign datetime array
values{:,3:end} = val_arr;
current{:, :} = curr_arr;
values_time = array2table(zeros(height(good_time),2));
values_time.Properties.VariableNames = {'t_start', 't_end'};
values_time = array2table(values{:,1:2});
valtime_ltx = table2latex(values_time);


values_ltx = table2latex(values);
fid = fopen('events/tables/values_ltx.txt', 'w');  % Open the file for writing
fprintf(fid, '%s\n', values_ltx);           % Write each line of text
fclose(fid);  
current_ltx = array2latex(current);
fid = fopen('events/tables/current_ltx.txt', 'w');  % Open the file for writing
fprintf(fid, '%s\n', current_ltx);           % Write each line of text
fclose(fid);  

% list event numbers tints
    event_tints = array2table(zeros(height(good_time),3));
    event_tints.Properties.VariableNames = {'Event','Start','End'};
    event_tints{:,1}= (1:length(good_time))';
    event_tints{:, 2} = string(values{:,1});
    event_tints{:, 3} = string(values{:,2});
    % make table
    times_ltx = table2latex(event_tints);
    fid = fopen('events/tables/event_tint.txt', 'w');  % Open the file for writing
    fprintf(fid, '%s\n', times_ltx);           % Write each line of text
    fclose(fid);  


% calculate Jpks_>0.25/Jpks
parts_ab25 = curr_arr(:,7)./curr_arr(:,5);
newcurr_arr = zeros(length(curr_arr), 6);
newcurr_arr(:,1:5) = curr_arr(:,1:5);
newcurr_arr(:,6) = parts_ab25*100;

% create matrix for color 
comat_arr = newcurr_arr(:,2:end);      
minusy = find(newcurr_arr(:,1)==-1);
comat_arr(minusy,3:4) = newcurr_arr(minusy,4:5)*(-1); % changing sign of values on -y side for unified colormap

colorMatrix = array2table(comat_arr);
eventMatrix = newcurr_arr(:,2:end);   % Matrix for displaying numbers (with correct sign of y)
criteriaLabels = {'(a) $(\mathbf{\hat{J}}\cdot\mathbf{\hat{r}})_{all} [\%]$', '(b) $(\mathbf{\hat{J}}\cdot\mathbf{\hat{r}})_{peaks} [\%]$', '(c) $(\mathbf{\vec{J}}\cdot\mathbf{\hat{r}})_{all} [\mu A/ m^2]$', ...
    '(d) $(\mathbf{\vec{J}}\cdot\mathbf{\hat{r}})_{peaks} [\mu A/ m^2]$', '(e) $(\frac{\mathbf{\vec{J}}\cdot\mathbf{\hat{r}}_{> 0.25}}{\mathbf{\vec{J}}\cdot\mathbf{\hat{r}}})_{peaks} [\%]$'};
eventLabels = 1:length(good_time);

%manually exclude event 13 because current values exceed 
comatwo13 = comat_arr;
comatwo13(13, :) = 0;

% create colorbar limits for each column
limit_perc = [0,100];
limit_all  = [max(abs(comatwo13(:,3)))*-1,max(abs(comatwo13(:,3)))];
limit_pks  = [max(abs(comatwo13(:,4)))*-1,max(abs(comatwo13(:,4)))];
limit_perc_special = [-100,100]; % includes negative percentage because 
% currents above 0.25 might point in different direction than mean

colorlimitsArray = {limit_perc, limit_perc, limit_all, limit_pks, limit_perc_special}; % limit_perc, limit_perc, limit_all, limit_pks, limit_alla25, limit_pksa25
% make correlation table
correlation_matrix_pmy(colorMatrix, eventMatrix, eventLabels, criteriaLabels, colorlimitsArray)
