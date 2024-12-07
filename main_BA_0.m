% main document.
irf; 
addpath('plot');
addpath('my_functions');
mms.db_init('local_file_db','/mnt/mms'); % initialize path to mms data
% access to MMS data is needed for this program. Either mount to a server
% or download the raw data from NASA or use HAPI for matlab to download. 

tic;
Param.SCall = [1,2,3,4]; % defining global parameters
Param.ic = 1;
Param.RE = 6378;
% make database
get_database_1;

% filter MMS 
% Q_param, tetr, amount of peaks are in right
% range
fprintf('\n Now checking for Q, tetr, cone_a... ');
ok_arr = filter_mms_data(good_time, Param);
good_time = good_time(find(all(ok_arr ~= 0, 2)),:);
fprintf('Finished! \n');
% preallocate arrays
params_preallocations;
save('events/time_intervals.mat', 'good_time');

%% stat. analyze extracted time intervals and plot functions
for i = 1:height(good_time)
    fprintf('\n Starting to analyze time interval i = %d\n', i);
    tmp.start = irf_time(good_time(i,1), 'epoch>epochtt');
    tmp.stop = irf_time(good_time(i,2), 'epoch>epochtt');
    tmp.tint = irf.tint(tmp.start, tmp.stop); 
    tmp.tint_beob = [tmp.tint(1)+(15*60),tmp.tint(2)]; % maybe change tmp.stop to tmp.start+25*60? - so that all tints have same length!
    % analyze tint
    [omni, mmsd, event, val_arr(i,:), curr_arr(i,:)] = load_analyze_tint_2(tmp.tint, tmp.tint_beob, Param, i);
    
    % save data to plot with python
    % omni, mmsd, event, tint, tint_beob, Param
    event.tint_string = irf_fname(tmp.tint_beob);
    tmp.filename = fullfile('events/data');
    if ~exist(tmp.filename, 'dir')
        mkdir(tmp.filename);
    end
    save(['events/data/', event.tint_string,'_omni_data.mat'], 'omni');
    save(['events/data/', event.tint_string,'_mms_data.mat'], 'event');
    
    % plot tint
    if 0
        plot_tint_3(omni, mmsd, event,tmp.tint, tmp.tint_beob, Param);
    end
    fprintf('\nTime interval number %d/%d: Finished. \n', i, length(good_time));
end

%% make tables
analyze_dataset_4;

counter.exe = toc;
%% --------------------------------------------------------------------------
fprintf('Execution time: %.2f minutes\n', counter.exe/60);
clearvars -except good_time Param values val_arr current curr_arr ok_arr;
