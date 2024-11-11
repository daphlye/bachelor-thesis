% define parameters and preallocate "event" struct

% preallocation (size of number of tints)
values = array2table(zeros(height(good_time),16));
values.Properties.VariableNames = ...
    {'t_start', 't_end', 'cone_a', 'cone_a_ok', 'j_maguv_all', 'j_maguv_pks','j_uv_all',...
    'j_uv_pks', 'j_ok', 'nslv', 'Rx', 'Ry', 'Rz', 'tetr_E', 'tetr_P','tetr_ok'};