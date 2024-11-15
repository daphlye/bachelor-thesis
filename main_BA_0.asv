% main document.
irf; 
mms.db_init('local_file_db','/mnt/mms'); %initialize path to mms data
tic;
Param.SCall = [1,2,3,4]; % defining global parameters
Param.ic = 1;
Param.RE = 6378;
% make database
get_database_1;
% preallocate the events array
params_preallocations;
% analyze extracted time intervals and plot functions
counter.rd = 0; % counts events in which more than 50 percent of the data 
% points are pointing in the right direction

for i = 1:height(good_time)
    
    tmp.start = irf_time(good_time(i,1), 'epoch>epochtt');
    tmp.stop = irf_time(good_time(i,2), 'epoch>epochtt');
    tmp.tint = irf.tint(tmp.start, tmp.stop); 
    tmp.tint_beob = [tmp.tint(1)+(15*60),tmp.tint(2)]; % maybe change tmp.stop to tmp.start+25*60? - so that all tints have same length!
    values{i,1:2} = [good_time(i,1)+(15*60),good_time(i,2)];
    % analyze tint
    [omni, mmsd, event, val_arr] = load_analyze_tint_2(tmp.tint, tmp.tint_beob, Param, i);
    values{i,3:end} = val_arr;
    % if elongation, planarity are okay 
    % plot tint
    if 1 %values.cone_a_ok(i)==1 && values.j_ok(i)==1 && values.tetr_ok(i)==1
        plot_tint_3(omni, mmsd, event, tmp.tint_beob, tmp.tint, Param);
    end

    fprintf('\nTime interval number %d/%d: Finished. \n', i, length(good_time));
end
save('events/analysis_result.mat', 'values');
uitable('Data', values{:,:}, 'ColumnName', values.Properties.VariableNames);
saveas(gcf, 'events/analysis_result.mat.png');
% analyze val_arr by counting how many tints fullfilled different
% criteria FUNCTION
analyze_dataset_4;



counter.exe = toc;
%--------------------------------------------------------------------------
fprintf('Execution time: %.2f minutes\n', counter.exe/60);
clearvars -except counter event good_time mmsd omni Param tints values tints_ana_summary;

%% 
clearvars -except tints good_time counter Param;

%%
for i = 1:3
    test.structure(i,:)= mod(i, 2) == 0;
end