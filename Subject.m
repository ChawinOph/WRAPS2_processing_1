classdef Subject < handle
    %Subject: Class of a human subject
    %   Detailed explanation goes here
    properties (Constant)
        vicon_folder_name = 'Vicon data'; % the name in the directory that leads to the
        freq_fplate       = 1000;
        freq_marker       = 100;
    end
    
    properties (Access = public)
        sbj_folder_name = '';               % string of a folder name that contains all subject data, e.g., 'w001_2018-09-18'
        sbj_index = 0;                      % index number of subject
        sbj_id = '';                        % string of a subject id, e.g., w001.
        sbj_name = '';                      % string of a subject name, e.g., 'Chawin'
        sbj_anthro_measurement              % struc of subject measurement data
        sbj_marker_cluster_pos = struct();  % struc of marker clusters on subject 
        raw_data = struct()               % stucture of recoded raw/processed data from the vicon
        sbj_anthro                          % object of human model containing transforms of all body segments
        sbj_WRAPS2 = WRAPS_2()              % object of WRAPS on the subject containing transforms of rings from CAD
    end
    
    methods
        %% Constructor
        function this = Subject(sbj_index, sbj_folder_names, sbj_ids, sbj_names, marker_cluster_pos)
            %Subject Construct an instance of this class
            %   Detailed explanation goes here
            this.sbj_folder_name = sbj_folder_names{sbj_index};
            this.sbj_index = sbj_index;
            this.sbj_id = sbj_ids{sbj_index};
            this.sbj_name = sbj_names{sbj_index};
            this.sbj_anthro_measurement.sbj_id = this.sbj_id;
            this.sbj_anthro_measurement.sbj_name = this.sbj_name;
            this.sbj_anthro_measurement.weight_kg = 72; % kg
            this.sbj_anthro_measurement.height_cm = 173; % cm
            this.sbj_marker_cluster_pos = this.importMarkerClusterPos(marker_cluster_pos);
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
            this.raw_data(trial_no).marker_trial_name =  trial_file_names{this.sbj_index};
            
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
                this.raw_data(trial_no).marker_data(seg_no).segment_names = sorted_segment_names{seg_no};
                this.raw_data(trial_no).marker_data(seg_no).marker_names = sorted_marker_names{seg_no};
                this.raw_data(trial_no).marker_data(seg_no).marker_pos =  marker_pos;
            end
            
            disp('Imported raw marker data')
            calcTransformation(this, trial_no, 'Pelvis')
            calcTransformation(this, trial_no, 'Thorax')
            
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
            
            this.raw_data(trial_no).fplate_trial_name =  trial_file_names{this.sbj_index};
            
            % store data base on the force plate and variable names
            for plate_no = 1:length(sorted_forceplate_names)
                this.raw_data(trial_no).fplate_data(plate_no).fplate_name = sorted_forceplate_names(plate_no);
                this.raw_data(trial_no).fplate_data(plate_no).fplate_var_names = sorted_forceplate_var_names;
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
                    this.raw_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no) = fplate_var;                    
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
        
        function sbj_marker_cluster_pos = importMarkerClusterPos(this, marker_cluster_pos)
            n = strcmp({marker_cluster_pos.sbj_id}, this.sbj_id);
            sbj_marker_cluster_pos = marker_cluster_pos(n).cluster;
        end
                
        %% Calculate Transformations of marker clusters in the specified trial   
        
        function calcTransformation(this, trial_no, segment_name)
            % match the cluster/segment names
            segment_no = strcmp({this.raw_data(trial_no).marker_data.segment_names}, segment_name);
            cluster_no = strcmp({this.sbj_marker_cluster_pos.cluster_names}, segment_name);
            used_marker_indcs = zeros(length(this.sbj_marker_cluster_pos(cluster_no).marker_names), 1);
            for i = 1: length(used_marker_indcs)
                % find the indcs of all recorded markers begin used in the
                % cluster
                used_marker_indcs(i) = find(strcmp(this.raw_data(trial_no).marker_data(segment_no).marker_names, ...
                    this.sbj_marker_cluster_pos(cluster_no).marker_names{i}));
            end
            
            % store the 3d matrix of marker pos from the trial in the same order as
            % the static cluster pos
            vicon_pos = this.raw_data(trial_no).marker_data(segment_no).marker_pos(:, :, used_marker_indcs);
            cluster_pos = this.sbj_marker_cluster_pos(cluster_no).marker_pos;        
            
            % Use Least Square Rigid Body Motion by SVD (http://www.igl.ethz.ch/projects/ARAP/svd_rot.pdf)
            cluster_pos_centroid = mean(cluster_pos)'; %
            X = cluster_pos' - cluster_pos_centroid;
            T_v2s = zeros(4, 4, length(vicon_pos)); % tranform from v to a segment
            
            for i = 1:length(vicon_pos)
                % get the set of marker pos at the current time step
                curr_vicon_pos = reshape(vicon_pos(i,:,:), size(vicon_pos, 2), []);
                
                % find the column vector centroid of the vicon pos
                curr_vicon_pos_centroid = mean(curr_vicon_pos, 2);
                Y = curr_vicon_pos - curr_vicon_pos_centroid;
                
                % Calculate the 3 x 3 covariance matrix
                S = X*Y';
                
                % Comput the SVD: S = U*Sigma*V'
                [U, Sigma, V] = svd(S);
                M = eye(size(Sigma, 1)); M(end,end) = det(V*U'); % correct reflection
                R = V*M*U'; % orthogonal rotational matrix
                
                % construct transformation from the vicon origin frame
                T_v2s(:, :, i) = eye(4);
                T_v2s(1:3, 1:3, i) = R;
                T_v2s(1:3, 4, i) =  curr_vicon_pos_centroid - R*cluster_pos_centroid;
            end
            
            this.raw_data(trial_no).marker_data(segment_no).used_marker_names = ...
                this.sbj_marker_cluster_pos(cluster_no).marker_names;
            this.raw_data(trial_no).marker_data(segment_no).used_marker_pos = vicon_pos;
            this.raw_data(trial_no).marker_data(segment_no).transforms_vicon2seg= T_v2s;
            disp(['Updated ', segment_name, ' segment transformations in trial no. ', num2str(trial_no)])
        end
       
        %% Visualization
        
        function plotCoPvsTime(this, trial_no, plate_name)
            var_name = 'CoP';
            plate_no = find(strcmp([this.raw_data(trial_no).fplate_data.fplate_name], plate_name));
            var_no = find(strcmp(this.raw_data(trial_no).fplate_data(plate_no).fplate_var_names, var_name));
            var = this.raw_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no); %#ok<FNDSB>
            t = (0:1:length(var)-1)/this.freq_fplate;
            plot(t, var);
            xlim([0 max(t)]);
            title([plate_name, ': ', var_name])
        end
        
        function plotFvsTime(this, trial_no, plate_name)
            var_name = 'Force';
            plate_no = find(strcmp([this.raw_data(trial_no).fplate_data.fplate_name], plate_name));
            var_no = find(strcmp(this.raw_data(trial_no).fplate_data(plate_no).fplate_var_names, var_name));
            var = this.raw_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no); %#ok<FNDSB>
            t = (0:1:length(var)-1)/this.freq_fplate;
            plot(t, var);
            xlim([0 max(t)]);
            title([plate_name, ': ', var_name])
        end
        
        function plotTrajCoP(this, trial_no, plate_name)
            var_name = 'CoP';
            plate_no = find(strcmp([this.raw_data(trial_no).fplate_data.fplate_name], plate_name));
            var_no = find(strcmp(this.raw_data(trial_no).fplate_data(plate_no).fplate_var_names, var_name));
            var = this.raw_data(trial_no).fplate_data(plate_no).fplate_var(:, :, var_no); %#ok<FNDSB>
%             plot(var(:,1)-var(1,1), var(:,2)-var(1,2));
%             plot(var(:,1), var(:,2));
            scat = scatter3(var(:,1), var(:,2), var(:, 3));
            scat.Marker = '.';
            grid on; grid minor;
            axis equal
            title([plate_name, ': CoP Trajectory'])
        end
                
        function vizTrial(this, trial_no)
            % vizTrial: visualization of each trial
            figure;
            this.plotCoordinatesTransform(eye(4), 250); hold on;
            this.plotTrajCoP(trial_no, 'Seat Plate');
            this.plotTrajCoP(trial_no, 'Foot Plate');
            
            % show initial pos of marker clusters of both rings
            pelvis_segment_no = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Pelvis');
            thorax_segment_no = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Thorax');
            
            init_pelvis_marker_pos = reshape(this.raw_data(trial_no).marker_data(pelvis_segment_no).used_marker_pos(1, :, :), 3, []);
            init_thorax_marker_pos = reshape(this.raw_data(trial_no).marker_data(thorax_segment_no).used_marker_pos(1, :, :), 3, []);
            scatter3(init_pelvis_marker_pos(1, :), init_pelvis_marker_pos(2, :), init_pelvis_marker_pos(3, :))
%             plot3(init_pelvis_marker_pos(1, :), init_pelvis_marker_pos(2, :), init_pelvis_marker_pos(3, :))
            scatter3(init_thorax_marker_pos(1, :), init_thorax_marker_pos(2, :), init_thorax_marker_pos(3, :))
            
%             extreme_time_step = 1000;
%             extreme_pelvis_marker_pos = reshape(this.raw_data(trial_no).marker_data(pelvis_segment_no).used_marker_pos(extreme_time_step, :, :), 3, []);
%             extreme_thorax_marker_pos = reshape(this.raw_data(trial_no).marker_data(thorax_segment_no).used_marker_pos(extreme_time_step, :, :), 3, []);
%             scatter3(extreme_pelvis_marker_pos(1, :), extreme_pelvis_marker_pos(2, :), extreme_pelvis_marker_pos(3, :))
%             scatter3(extreme_thorax_marker_pos(1, :), extreme_thorax_marker_pos(2, :), extreme_thorax_marker_pos(3, :))
            
            title(this.raw_data(trial_no).marker_trial_name)
            for time_step = 1: 100 : length(this.raw_data(trial_no).marker_data(1).transforms_vicon2seg)
                T_v2pelvis = this.raw_data(trial_no).marker_data(pelvis_segment_no).transforms_vicon2seg(:,:,time_step);
                T_v2thorax = this.raw_data(trial_no).marker_data(thorax_segment_no).transforms_vicon2seg(:,:,time_step);
                this.plotCoordinatesTransform(T_v2pelvis, 100);
                this.plotCoordinatesTransform(T_v2thorax, 100);
            end
        end      
        
        % visual elements
        
        function plotCoordinatesTransform(~, T, scale)
            % PlotCoordinate: Plot coordinates in the XYZ-RGB sequence
            origin = T(1:3,4);
            unit_vecs = T(1:3, 1:3);
            axis_x = unit_vecs(:,1);
            axis_y = unit_vecs(:,2);
            axis_z = unit_vecs(:,3);
            x = origin(1); y = origin(2); z = origin(3);
            quiver3(x,y,z,axis_x(1),axis_x(2),axis_x(3),scale,'color','r'); hold on
            quiver3(x,y,z,axis_y(1),axis_y(2),axis_y(3),scale,'color','g')
            quiver3(x,y,z,axis_z(1),axis_z(2),axis_z(3),scale,'color','b')
        end
        
    end
end

