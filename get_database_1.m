%% BIG_BABA_BA_SCRIPT
% ====================================
% this is the script that documents nicely and in the approximately correct
% order, how I filter OMNI data with the criteria that MMS is in the
% dayside magnetosheath and that the solar wind is as stable as possible
%
% ---> currently strong Bz and np criteria applied! for more results lower these
%      criteria in the if loops!
%
% ====================================
%% Load data 

% Input parameters

trange.yr1 = irf.tint('2015-09-01T00:00:00/2016-06-01T00:00:00');
trange.yr2 = irf.tint('2016-09-01T00:00:00/2017-06-01T00:00:00');
trange.yr3 = irf.tint('2017-09-01T00:00:00/2018-06-01T00:00:00');
trange.yr4 = irf.tint('2018-09-01T00:00:00/2019-06-01T00:00:00');
trange.yr5 = irf.tint('2019-09-01T00:00:00/2020-06-01T00:00:00');
trange.yr6 = irf.tint('2020-09-01T00:00:00/2021-07-01T00:00:00');
trange.yr7 = irf.tint('2021-09-01T00:00:00/2022-07-01T00:00:00');
trange.yr8 = irf.tint('2022-09-01T00:00:00/2023-07-01T00:00:00');
% irf.tint('2015-09-01T00:00:00/2023-07-01T00:00:00'); too long - not more
% than 600000 records allowed to download
% all data is included since no magnetosheath crossings between june/july
% and september

%-------------------------------------

% Step 1: Read the CSV file
filename = 'inp/magnetosheath_region_list.csv';  % find the Database at https://doi.org/10.5281/zenodo.10491878
opts = detectImportOptions(filename);
tmp.data = readtable(filename, opts); % now look at CSV file
data = tmp.data(1:2451,:); % all entries for when mms1 is in ms (entries below are mms2,3,4)

% Step 2: Convert to Unix time
time_start  = EpochUnix(posixtime(data{:,1}));
time_stop  = EpochUnix(posixtime(data{:,2}));

%=====================================
% Step 4: Apply the criteria function
% Load OMNI data for the extracted tints and save the timestamps of the
% ones that survive the criteria function
%=====================================
% Testvalues
years = fieldnames(trange);
counter.Ms = 0;
counter.B = 0;
counter.v = 0;
counter.np = 0;
counter.good_tint = 0;
counter.bad_tint = 0;
good_time = [];
for k = 1:length(years) % downloading and whole routine for all years
    fprintf('\n Loading Year %d data. \n',k);
    tint = trange.(years{k});
    tmp.B = irf_get_data_omni(tint, 'b,bx,byGSM,bzGSM', 'omni_min');% loads 1-min OMNI data flow speed velocity in gsm coordinates (OMNI v data is in gse) for selected tint
    if isempty(tmp.B)
        fprintf('\n \n ---- ATTENTION ---- \n Download of OMNI data from nasa.gov is not working. \n Don''t download more than 600000 entrys \n');
        return;
    end
    tmp.V = irf_gse2gsm(irf_get_data_omni(tint, 'v,vx,vy,vz', 'omni_min'));
    tmp.Ms_n = irf_get_data_omni(tint, 'Ms,n', 'omni_min'); % 1 AU IP Magnetosonic Mach number, ion dens
    tmp.IND = find(tmp.Ms_n(:,2)>99);    tmp.Ms_n(tmp.IND,2) = NaN; % removing fillvalues
    time = EpochUnix(tmp.V(:,1:1)); % convert downloaded time to unix data
    B_tot = TSeries(time, [tmp.B(:,2:2)], 'to', 1); % GSE [nT]
    Bxyz = TSeries(time, [tmp.B(:,3:5)], 'to', 1); % GSM [nT]
    V_tot = TSeries(time, [tmp.V(:,2:2)], 'to', 1); % GSM [km/s]
    Vxyz = TSeries(time, [tmp.V(:,3:5)], 'to', 1); % GSM [km/s]
    Ms = TSeries(time, [tmp.Ms_n(:,2:2)], 'to', 1); % [#]
    np = TSeries(time, [tmp.Ms_n(:,3:3)], 'to', 1); % [cc]
    clear tmp;

    fprintf('Download succeeded. Data is being processed... ');

    %function
    good_intervals.(years{k}) = [];
    for i = 1:length(time_start) % moves through all timeintervalls the paper mentioned
        time_start_omni = time_start.epoch(i)-(15*60); %omni_tint shall start 15mins before tint_ms
        time_stop_omni = time_stop.epoch(i);
        ind = find((time.epoch >= time_start_omni) & (time.epoch <= time_stop_omni)); % calculates indexes for data within timeframe
        counter.gtims = []; % counts tints that fulfill criteria within one tint_ms
        if ~isempty(ind) % checks if ms_tint is still within the downloaded omni_tint
            tint_omni = time(ind); % gives array of measurement times for whole time + 20mins before
    %--------------------------------------------------------------------------
            input = ind; window = 30; %window of 30mins
            n = length(input);
            wind.good_ind = zeros((n - window + 1), 2); % includes indices of the omni dataset that pass the criteria
            for j = 1:(n - window + 1) % loops until window looks at last frame of "input"
                subwind = input(j:j + window - 1);
    % Magnetosonic Mach Number
                wind.Ms = Ms.data(subwind);
        % Data gaps && M_MS < 5
                if (NaNalyzer(wind.Ms, 0.1) == 1) && all(wind.Ms(~isnan(wind.Ms))<5) % data gaps and Ms Machnumber 
                    counter.Ms = counter.Ms + 1;
    % Magnetic field 
                    wind.Bz = Bxyz.data(subwind,3:3);
                    wind.ca = atan2d(Bxyz.data(subwind,2:2),Bxyz.data(subwind,3:3)); % calc. clock angle between By,z
                    % takes into account, that B comp can be negative
        % Data gaps &&  Bz variies less than 6nT && Bz<-15nT && CA close to
        % 180 deg && CA doesnt vary more than 30 deg
        % varys less than 30 deg
                    if (NaNalyzer(wind.Bz, 0.1) == 1) && (dev_check(wind.Bz, 0.1*mean(wind.Bz, 'omitnan')) == 1) && all(wind.Bz(~isnan(wind.Bz))<(-15)) && all(abs(wind.ca(~isnan(wind.ca)))>135) && (dev_check(wind.ca, 30) == 1) % its ok if I only check Bz - Bx, By become NaN simultaneously
                        % instead of 0.1*mean(wind.Bz, 'omitnan') try
                        % 6*1e(-9)
                        % - variation of 10% (dev_check) doesnt make too much sense.
                        counter.B = counter.B + 1;
    % Solar Wind speed
                        wind.v = V_tot.data(subwind);
                        if (NaNalyzer(wind.v, 0.1) == 1) && (dev_check(wind.v, 50) == 1) && all(wind.v(~isnan(wind.v))<(500))
                            counter.v = counter.v + 1;
    % Particle density
                            wind.np = np.data(subwind);
                            if (NaNalyzer(wind.np, 0.1) == 1) && (dev_check(wind.np, 0.1*mean(wind.np, 'omitnan')) == 1)
                                counter.np = counter.np + 1;
    % Finding the indeces for time intervals that pass the test
                                wind.good_ind(j, :) = [subwind(1),subwind(end)];
    % Could include: quasi perp. BS, algorithm to load MMS data and
    % control if in tetrahedal shape, (more ideas here)
                            end
                        end
                    end
                end
            end
            if any(wind.good_ind ~= 0) % only calculates an index if tint passed the criteria
                wind.counter=1; % counts how many windows are right after each other
                for l = 1:length(wind.good_ind) %connects tints that are 1min apart from each other
                    if wind.good_ind(l)~=0 %if it is 0 the tint didnt pass
                        if (wind.good_ind(l)+1) == wind.good_ind(l+1) % equal to: if wind.good_ind(l+1)~=0
                            wind.counter = wind.counter+1;
                        else
                            counter.gtims = [counter.gtims;[wind.good_ind(l-wind.counter+1,1:1), wind.good_ind(l,2:2)]];
                            wind.counter = 1; % reset, so that a possible second tint can be found. 
                        end
                    end
                end
            end
        end
        if ~isempty(counter.gtims)
            %counter.gtims = counter.gtims(counter.gtims ~= 0); %transforms
            %array weirdly
            good_intervals.(years{k}) = [good_intervals.(years{k});counter.gtims]; 
% Filter matching tints
    % find SC_in_MS times that lie within the OMNI timestamps 
            good_time = [good_time;time.epoch(good_intervals.(years{k})(:,1)),time.epoch(good_intervals.(years{k})(:,2))];
            
        end
    end
    fprintf('Finished! \n \n')
end

% Confirmation and presenting results
if time_start_omni == (time_start.epoch(end)-(15*60))
    fprintf('All time intervals checked. \n ---- number of windows passed ---- \n Magnetosonic Mach Number: %d \n IMF: %d \n SW speed: %d \n Density: %d \n---------------------------------- \n', counter.Ms, counter.B, counter.v, counter.np)
else 
    fprintf('Something went wrong.\n')
end

% converting all time intervals as irf_tint data into a stiPDruct called
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

clearvars -except tints good_time counter Param;

