function is_burst_available = burstmode_availability(sc, tint)
    % Checks if burst mode data is available for MMS spacecraft within a time interval
    %
    % Parameters:
    %   sc - spacecraft number ('1', '2', '3', or '4')
    %   tint - time interval as EpochTT object, e.g., tint = irf.tint('2015-10-16T10:00:00Z/2015-10-16T10:05:00Z');
    %
    % Returns:
    %   is_burst_available - logical value indicating if burst mode data is available (true or false)
    
    % Define data type and data source (adjust these parameters as needed)
    datatype = 'burst';
    varname = 'mms?_dfg_srvy_b';  % Magnetic field data in burst mode (modify if you are interested in other data types)
    mms_id = num2str(sc);  % Convert spacecraft number to string if needed
    
    % Replace ? in varname with the spacecraft number
    varname = strrep(varname, '?', mms_id);

    % Check data availability using mms.db_list_files
    try
        % List available burst mode data files for specified time interval
        file_list = mms.db_list_files(varname, tint);
        
        % If file list is not empty, burst mode data is available
        if ~isempty(file_list)
            is_burst_available = true;
            fprintf('Burst mode data is available for MMS%s in the specified interval.\n', sc);
        else
            is_burst_available = false;
            fprintf('No burst mode data available for MMS%s in the specified interval.\n', sc);
        end
    catch
        % Handle error if mms.db_list_files fails
        is_burst_available = false;
        fprintf('Error accessing burst mode data for MMS%s.\n', sc);
    end
end