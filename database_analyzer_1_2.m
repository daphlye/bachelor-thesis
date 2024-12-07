%% database analyzer
% Confirmation and presenting results
% this part is not needed anymore!
if time_start_omni == (time_start.epoch(end)-(15*60))
    fprintf('All time intervals checked. \n ---- number of windows passed ---- \n Magnetosonic Mach Number: %d \n IMF: %d \n SW speed: %d \n Density: %d \n---------------------------------- \n', counter.Ms, counter.B, counter.v, counter.np)
else 
    fprintf('Something went wrong.\n')
end

% converting all time intervals as irf_tint data into a struct called
% tints:
for i = 1:length(good_time)
    tmp.fieldName = ['tint' num2str(i)];
    tmp.start = irf_time(good_time(i,1), 'epoch>epochtt');
    tmp.stop = irf_time(good_time(i,2), 'epoch>epochtt');
    tints.(tmp.fieldName) = irf.tint(tmp.start, tmp.stop); 
end

% check if all timeintervals are within magnetosheath tints in the "data"
% table
tmp.ntint = fieldnames(tints);
for i = 1:numel(tmp.ntint)
    tmp.tint = tints.(tmp.ntint{i});
    if length(tmp.tint)==2
        for j = 1:height(data)
            if all((tmp.tint+15*60)>=EpochUnix(posixtime(data{j,1})) & tmp.tint<=EpochUnix(posixtime(data{j,2})))
                counter.good_tint = counter.good_tint + 1;
            end
        end
        if counter.good_tint == counter.bad_tint+1
            counter.bad_tint = counter.bad_tint+1;
        elseif counter.good_tint == counter.bad_tint
            fprintf('Epoch data in field %s is NOT within any time ranges in the data table. \n', tmp.ntint{i});
        else
            fprintf("something I dont understand is happening. look at code in line 170ff.\n");
        end
    else
        fprintf('the tint %d in tints does not have two entries.\n', i)
    end
end
if counter.good_tint == counter.bad_tint
    fprintf("everything worked out fine:)\n"); 
    fprintf('Extracted %d time intervals that fit the criteria. \n', counter.good_tint);
else
    fprintf("sth went wrong. \n");
end

% workspace info
workspace_info = whos;
total_bytes = sum([workspace_info.bytes])*10^-6;
fprintf('Total memory used: %.2f MB ... Giving free workspace. \n', total_bytes); 

