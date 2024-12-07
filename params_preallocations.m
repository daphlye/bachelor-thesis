% define parameters and preallocate 

% preallocation for "values" (size of number of tints)
val_arr = zeros(height(good_time),7);
values = array2table(zeros(height(good_time),9));
values.Properties.VariableNames = ...
    {'t_start', 't_end', 'cone_a', 'nslv', 'Rx', 'Ry', 'Rz', 'tetr_E', 'tetr_P'};
values.t_start = NaT(height(good_time), 1); % Preallocate with NaT for datetime
values.t_end = NaT(height(good_time), 1);   % Preallocate with NaT for datetime

% preallocation for "current_info" (size of number of tints)
curr_arr = zeros(height(good_time),7);
current = array2table(zeros(height(good_time),7));
current.Properties.VariableNames = ...
    {'right dir.', 'j_uv_all %', 'j_uv_pks %', 'tot_curr_all', ...
    'tot_curr_pks', 'tot_curr_a25_all', 'tot_curr_a25_pks'};