function ok_arr = filter_mms_data(good_time, Param)
    %% check current quality factor (divovercurl)
    
    % preallocation for ok_arr
    ok_arr = zeros(height(good_time),4);
    
    for i = 1:height(good_time)
        tmp.start = irf_time(good_time(i,1), 'epoch>epochtt');
        tmp.stop = irf_time(good_time(i,2), 'epoch>epochtt');
        tmp.tint = irf.tint(tmp.start, tmp.stop); 
        tint = [tmp.tint(1)+(15*60),tmp.tint(2)];

    % downloading needed daata
        c_eval('event.B?_gsm = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gsm_srvy_l2'',tint);',Param.SCall);
        % resampling on B since this has the highest resolution
        c_eval('event.B?_gsm = event.B?_gsm.resample(event.B1_gsm);',2:4);
        c_eval('event.R?_gsm = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_r_gsm_srvy_l2'',tint);', Param.SCall); % GSM [km] evaluates data for R 
        c_eval('event.R?_gsm = event.R?_gsm.resample(event.R1_gsm);',2:4);
        c_eval('event.Rres?  = event.R?_gsm.resample(event.B1_gsm);');
    
    % Curlometer technique (Method from Schwartz et al. 1998)
        % J, divB, jxB, divTshear, divPb (from Example_MMS_B_E_J.m)
        [event.j,event.raw.f.divB,event.raw.f.B,event.raw.f.jxB,event.raw.f.divTshear,event.raw.f.divPb] = c_4_j('event.Rres?','event.B?_gsm');
        % uses first three rows for x y z. (confirmed looking at c_4_j.m >
        % c_4_grad.m)
    
        % extract curlometer datapoints where divovercurl is < 1
        event.raw.divovercurl = event.raw.f.divB;     event.raw.divovercurl.data = abs(event.raw.divovercurl.data)./event.j.abs.data; % error on j in percentage of value
        % Test if j error is low enough to use j data
        j_ok = mean(event.raw.divovercurl.data, 'omitnan') < 0.3;
    
        if 1 %j_ok == 1
        % check tetraeder formation
            % confirming mms is in tetrahedral formation
            [tetr.volTensor,tetr.R_Center,tetr.dR1,tetr.dR2,tetr.dR3,tetr.dR4,tetr.L,tetr.E,tetr.P] = c_4_r(event.R1_gsm.data(:,2:4),event.R2_gsm.data(:,2:4),event.R3_gsm.data(:,2:4),event.R4_gsm.data(:,2:4)); %time is same as mmsd.R?_gsm
            tetr_E = mean(tetr.E, 'omitnan'); tetr_P = mean(tetr.P, 'omitnan');
            
            tetr_ok = mean(tetr.E, 'omitnan') < 0.3 && mean(tetr.P, 'omitnan') < 0.3;
        end
        
        if 1 %j_ok == 1 && tetr_ok == 1
        % check amount of peaks with stddev and pkprom
            %event.j_mag.data = sqrt(event.j.data(:,1).^2 + event.j.data(:,2).^2 + event.j.data(:,3).^2);
            event.j_mag = TSeries(event.j.time, sqrt(event.j.data(:,1).^2 + event.j.data(:,2).^2 + event.j.data(:,3).^2));
            [event.avgj, event.medj, event.stdj] = get_stats(event.j_mag.data);
            % Find peaks in Jmag and create corresponding timeseries
            nslv = event.medj + 2*event.stdj; % event.stdj*2
            % nslv = 2*event.stdj; 
            pkprom = 0.07*1e-6;
            [event.pksj,event.locsj]=findpeaks(event.j_mag.data,'MinPeakHeight',nslv, 'MinPeakProminence',pkprom); %'MinPeakProminence',0.07*1e-6
            pkprom_ok = ~isempty(event.locsj); 
        end
        if 1
            % Test if mean cone angle is quasi perpendicular (OMNI)
            tmp.B = irf_get_data_omni(tint, 'b,bx,byGSM,bzGSM', 'omni_min');% loads 1-min OMNI data flow speed velocity in gsm coordinates (OMNI v data is in gse) for selected tint
            omni.time = EpochUnix(tmp.B(:,1:1)); % convert downloaded time to unix data
            omni.Bxyz = TSeries(omni.time, [tmp.B(:,3:5)], 'to', 1); % GSM [nT]
            cone_a = mean(acosd(omni.Bxyz.data(:,1)./sqrt(omni.Bxyz.data(:,1).^2+omni.Bxyz.data(:,2).^2+omni.Bxyz.data(:,3).^2)), 'omitnan');
            cone_a_ok = cone_a > 60 && cone_a < 120;
        end
        ok_arr(i,:) = [j_ok, tetr_ok, pkprom_ok, cone_a_ok];
    end
end