% Step 1: Read the CSV file
filename = 'magnetosheath_region_list.csv';  % Replace with your CSV file path
opts = detectImportOptions(filename);
data = readtable(filename, opts); % now look at CSV file

% Step 2: Convert to Unix time
time_start  = EpochUnix(posixtime(data{:,1}));
time_stop  = EpochUnix(posixtime(data{:,2}));
%% Find longest time interval
for i = 1:9649
    tint = irf_time(time_stop(i), 'epochtt>tt')-irf_time(time_start(i), 'epochtt>tt');
    ref_tint=0;
    if tint>ref_tint
        ref_tint=tint;
    end
end
ref_tint./60
% RESULT: 217 minutes is longest time that scis in dayside MS
