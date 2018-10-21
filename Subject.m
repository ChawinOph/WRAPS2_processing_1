classdef Subject < handle
    %Subject: Class of a human subject
    %   Detailed explanation goes here
    properties (Constant)
        vicon_folder_name = 'Vicon data'; % the name in the directory that leads to the
        freq_fplate       = 1000;
        freq_marker       = 100;
    end
    
    properties (Access = public)
        sbj_folder_name = '';            % string of a folder name that contains all subject data, e.g., 'w001_2018-09-18'
        sbj_index = 0;                   % index number of subject
        sbj_id = '';                     % string of a subject id, e.g., w001.
        sbj_name = '';                   % string of a subject name, e.g., 'Chawin'
        weight = 72;                     % kg
        height = 173;                    % cm
        trial_data = struct()            % stucture of recoded raw data from the vicon
        sbj_anthro_models = {}           % 
        sbj_WRAPS2_models = {}           % 
    end
    
    methods
        %% Constructor
        function this = Subject(sbj_index, sbj_folder_names, sbj_ids, sbj_names)
            %Subject Construct an instance of this class
            %   Detailed explanation goes here
            this.sbj_folder_name = sbj_folder_names{sbj_index};
            this.sbj_index = sbj_index;
            this.sbj_id = sbj_ids{sbj_index};
            this.sbj_name = sbj_names{sbj_index};
        end
        
        %% Member functions
        
        %% Import data
        function importMarkerData_csv(this, trial_no, trial_file_names, sorted_segment_names, sorted_marker_names)
            % ImportMarkerData
            %   foldername: folder that stores the text files imported from VICON
            %   textfilename: name of a text file (without .csv)
            %
            %   This script filters data with a 4th order butterworth low-pass filter, and
            %   then sorts by by marker or time.
            
            % store trial name
            this.trial_data(trial_no).marker_trial_name =  trial_file_names{this.sbj_index};
            
            % Setup file directory
            d = dir([this.sbj_folder_name,'\', this.vicon_folder_name, '\', trial_file_names{this.sbj_index}, '.csv']);
            filename = [d.folder '\' d.name];
            
            % Store raw matrix of marker positions
            M_raw = csvread(filename,5,2);
            
            % Apply filter
            M_trim = rmmissing(M_raw);
            
            % Create the low-pass filter (4th order Butterworth)
            Fs = 100; %this is the sampling frequency (frame rate)
            Fc = 5; %this is the cutoff frequency for the low-pass filter.
            Wn = (Fc*2)/Fs;
            [b,a] = butter(4, Wn);
            
            % Apply filter
            M_filt = filtfilt(b, a, M_trim);
            
            % Import headername for sorting
            M_headers = importdata(filename,',',3);
            
            % split header name
            raw_headernames = strsplit(char(M_headers{3,1}),',');
            
            % get rid of empty cells
            raw_headernames = raw_headernames(~cellfun('isempty',raw_headernames));
            
            % allocate each segment data to structure
            for seg_no = 1:length(sorted_marker_names)
                [var, var_indice] = this.extractMarkers(sorted_marker_names{seg_no}, raw_headernames, M_filt);
                marker_pos = this.sortMarker(var, length(var_indice));
                this.trial_data(trial_no).marker_data(seg_no).segment_names = sorted_segment_names{seg_no};
                this.trial_data(trial_no).marker_data(seg_no).marker_names = sorted_marker_names{seg_no};
                this.trial_data(trial_no).marker_data(seg_no).marker_pos =  marker_pos;
            end

            
        end
        
        function importForcePlateData_csv(this, trial_no, trial_file_names, sorted_forceplate_names, sorted_forceplate_var_names)
            % ImportForcePlateData
            %   foldername: folder that stores the text files imported from VICON
            %   textfilename: name of a text file (without .txt)
            %
            %   This script filters data with a 4th order butterworth low-pass filter, and
            %   then sorts by by marker or time.
            % Setup file directory
            d = dir([this.sbj_folder_name,'\', this.vicon_folder_name, '\', trial_file_names{this.sbj_index}, '.csv']);
            filename = [d.folder '\' d.name];
            
            % Store raw matrix of marker positions
            F = csvread(filename,5,2);
            
            % Trim any rows containing NaN and filter data
            F_trim = rmmissing(F);
            
            % Make the filter (4th order Butterworth)
            Fs = 1000; %this is the sampling frequency (frame rate)
            Fc = 8; %this is the cutoff frequency for the low-pass filter.
            Wn = (Fc*2)/Fs;
            [b,a] = butter(4, Wn);
            
            % Apply filter
            F_filt = filtfilt(b, a, F_trim);
            
            % Get a list of all marker names
            % Import headername for sorting
            F_headers = importdata(filename,',',3);
            
            % split header name
            headernames = strsplit(char(F_headers{3,1}),',');
            
            % get rid of empty cells
            headernames = headernames(~cellfun('isempty',headernames));
            
            this.trial_data(trial_no).fplate_trial_name =  trial_file_names{this.sbj_index};
            
            % store data base on the force plate and variable names
            for plate_no = 1:length(sorted_forceplate_names)
                this.trial_data(trial_no).fplate_data(plate_no).fplate_name = sorted_forceplate_names(plate_no);
                this.trial_data(trial_no).fplate_data(plate_no).fplate_var_names = sorted_forceplate_var_names;
                for var_no = 1:length(sorted_forceplate_var_names)
                    % Assign a variable name
                    var_indice = [];                   
                    for k = 1:size(headernames, 2)
                        if strcmp(headernames{k}, [sorted_forceplate_names{plate_no}, ' - ', sorted_forceplate_var_names{var_no}])
                            var_indice = [var_indice, k]; %#ok<AGROW>
                        end
                    end
                    fplate_var = [];
                    for i = 1:length(var_indice)
                        fplate_var = [fplate_var, F_filt(:, 3*(var_indice(i) - 1) + 1 : 3*(var_indice(i) - 1) + 3)]; %#ok<AGROW>
                    end
                    this.trial_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no) = fplate_var;                    
                end
            end
        end
        
        function [var, var_indice] = extractMarkers(this, marker_names, raw_headernames, M)
            % select only the positions of markers in the list
            var_indice = [];
            for i = 1:length(marker_names)
                for j = 1:size(raw_headernames, 2)
                    if strcmp(raw_headernames{j}, [this.sbj_name, ':', marker_names{i}])
                        var_indice = [var_indice, j]; %#ok<AGROW>
                    end
                end
            end
            var = [];
            for i = 1:length(var_indice)
                var = [var, M(:, 3*(var_indice(i)-1) + 1 : 3*(var_indice(i)-1) + 3)]; %#ok<AGROW>
            end
        end
        
        function markerdata = sortMarker(~, data, n_markers)
            [rows, ~] = size(data);
            markerdata = zeros(rows,3, n_markers);
            % Sort by marker
            for n = 0: n_markers - 1
                markerdata(:,:, n+1) = data(:, 3*n + 1 : 3*n + 3);
            end
        end
        
%         function importMarkerCAD(this, marker_position_CAD)
%             
%         end
                
        %% Calculate Transformations from marker and CAD data
        
        %% Visualization
        function plotCoPvsTime(this, trial_no, plate_name)
            var_name = 'CoP';
            plate_no = find(strcmp([this.trial_data(trial_no).fplate_data.fplate_name], plate_name));
            var_no = find(strcmp(this.trial_data(trial_no).fplate_data(plate_no).fplate_var_names, var_name));
            var = this.trial_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no); %#ok<FNDSB>
            t = (0:1:length(var)-1)/this.freq_fplate;
            plot(t, var);
            xlim([0 max(t)]);
            title([plate_name, ': ', var_name])
        end
        
        function plotFvsTime(this, trial_no, plate_name)
            var_name = 'Force';
            plate_no = find(strcmp([this.trial_data(trial_no).fplate_data.fplate_name], plate_name));
            var_no = find(strcmp(this.trial_data(trial_no).fplate_data(plate_no).fplate_var_names, var_name));
            var = this.trial_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no); %#ok<FNDSB>
            t = (0:1:length(var)-1)/this.freq_fplate;
            plot(t, var);
            xlim([0 max(t)]);
            title([plate_name, ': ', var_name])
        end
        
        function plotTrajCoP(this, trial_no, plate_name)
            var_name = 'CoP';
            plate_no = find(strcmp([this.trial_data(trial_no).fplate_data.fplate_name], plate_name));
            var_no = find(strcmp(this.trial_data(trial_no).fplate_data(plate_no).fplate_var_names, var_name));
            var = this.trial_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no); %#ok<FNDSB>
%             plot(var(:,1)-var(1,1), var(:,2)-var(1,2));
%             plot(var(:,1), var(:,2));
            scat = scatter3(var(:,1), var(:,2), var(:, 3));
            scat.Marker = '.';
            grid on; grid minor;
            axis equal
            title([plate_name, ': CoP Trajectory'])
        end
        
    end
end

