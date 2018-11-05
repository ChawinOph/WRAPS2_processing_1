classdef Subject < handle
    %Subject: Class of a human subject
    %   Detailed explanation goes here
    properties (Constant)
        vicon_folder_name = 'Vicon data'; % the name in the directory that leads to the
        freq_fplate       = 1000;
        freq_marker       = 100;
        g                 = [0, 0, -9.81]'; % relative to vicon frame
    end
    
    properties (Access = public)
        sbj_folder_name = '';                   % string of a folder name that contains all subject data, e.g., 'w001_2018-09-18'
        sbj_index = 0;                          % index number of subject
        sbj_id = '';                            % string of a subject id, e.g., w001.
        sbj_name = '';                          % string of a subject name, e.g., 'Chawin'
        sbj_gender = ''                         % string of gender ('male' or 'female')
        sbj_anthro_measurement                  % struct of subject measurement data
        sbj_anthro_table                        % stcuture of data from de Leva table
        sbj_marker_cluster_pos = struct();      % struct of marker clusters on subject
        sbj_landmark_config = struct();         % struct of landmark clusters on subject
        raw_data = struct()                     % stucture of recoded raw data from the vicon
        sbj_anthro = struct()                   % object of human model containing transforms of all body segments
        sbj_WRAPS2 = struct()                   % object of WRAPS on the subject containing transforms of rings from CAD
    end
    
    methods
        %% Constructor
        function this = Subject(sbj_index, sbj_folder_names, sbj_ids, sbj_names, sbj_gender, ...
                marker_cluster_pos, anthropomet_table, sbj_measurment_table)
            % Subject Construct an instance of this class
            %   Detailed explanation goes here
            this.sbj_folder_name = sbj_folder_names{sbj_index};
            this.sbj_index = sbj_index;
            this.sbj_id = sbj_ids{sbj_index};
            this.sbj_name = sbj_names{sbj_index};
            this.sbj_gender = sbj_gender{sbj_index};
            this.sbj_anthro_measurement.sbj_id = this.sbj_id;
            this.sbj_anthro_measurement.sbj_name = this.sbj_name;
            this.sbj_anthro_measurement.sbj_gender = this.sbj_gender;
            this.sbj_anthro_measurement.weight_kg = 72; % kg
            this.sbj_anthro_measurement.height_cm = 173; % cm
            this.importMarkerClusterPos(marker_cluster_pos);
            this.importAnthropometTable(anthropomet_table);
            this.importMeasurementTable(sbj_measurment_table);
        end
        
        %% Member functions
        
        
        %% Import data
        function importMarkerData_csv(this, trial_no, trial_file_names, sorted_segment_names, sorted_marker_names, removed_marker_names)
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
                if nargin > 5
                    % go through each set of marker names on each segment
                    % and remove the unwanted marker name
                    for i = 1:length(removed_marker_names)
                        sorted_marker_names{seg_no} = sorted_marker_names{seg_no}(~strcmp(sorted_marker_names{seg_no}, removed_marker_names{i}));
                    end
                end
                [var, var_indice] = this.extractMarkers(sorted_marker_names{seg_no}, raw_headernames, M_filt);
                marker_pos = this.sortMarker(var, length(var_indice));
                this.raw_data(trial_no).marker_data(seg_no).segment_names = sorted_segment_names{seg_no};
                this.raw_data(trial_no).marker_data(seg_no).marker_names = sorted_marker_names{seg_no};
                this.raw_data(trial_no).marker_data(seg_no).marker_pos =  marker_pos;
            end
            
            disp(['Imported raw marker data from trial no. ', num2str(trial_no), ' (', this.raw_data(trial_no).marker_trial_name, ')'])
            
            % calculate transformation of braces from raw vicon data
            this.sbj_WRAPS2(trial_no).trial_name = this.raw_data(trial_no).marker_trial_name;
            this.sbj_WRAPS2(trial_no).trial_transform_data = this.sbj_marker_cluster_pos;
            this.calcTransformation(trial_no, 'Pelvis', 'Pelvis Brace');
            this.calcTransformation(trial_no, 'Thorax', 'Thorax Brace');
            this.calcTransformationBetweenClusters(trial_no, 'Pelvis Brace', 'Thorax Brace')
            this.calcLandmarkPos(trial_no);
            this.calcSegmentTrans(trial_no);
            this.calcSegmentInertia(trial_no)
            disp(' ')
            
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
        
        function importMarkerClusterPos(this, marker_cluster_pos)
            n = strcmp({marker_cluster_pos.sbj_id}, this.sbj_id);
            this.sbj_marker_cluster_pos = marker_cluster_pos(n).WRAPS_cluster;
            this.sbj_landmark_config = marker_cluster_pos(n).landmark_config;
        end
        
        function importAnthropometTable(this, anthropomet_table)
            n = strcmp({anthropomet_table.gender}, this.sbj_gender);
            this.sbj_anthro_table = anthropomet_table(n).inertia_table;
        end
        
        function importMeasurementTable(this, sbj_measurment_table)
            n = strcmp({sbj_measurment_table.sbj_id}, this.sbj_id);
            this.sbj_anthro_measurement.measurement_table = sbj_measurment_table(n).measurement_table;
        end
        %% Calculate Transformations of marker clusters in a specified trial
        
        function calcTransformation(this, trial_no, vicon_segment_name, cluster_name)
            % match the cluster/segment names
            segment_no = find(strcmp({this.raw_data(trial_no).marker_data.segment_names}, vicon_segment_name));
            cluster_no = find(strcmp({this.sbj_marker_cluster_pos.cluster_name}, cluster_name));
            used_marker_indcs = [];
            
            for i = 1: length(this.sbj_marker_cluster_pos(cluster_no).marker_names)
                % find the indcs of all recorded markers used in cluster
                found_indx = find(strcmp(this.raw_data(trial_no).marker_data(segment_no).marker_names, ...
                    this.sbj_marker_cluster_pos(cluster_no).marker_names{i}));
                if found_indx ~= 0
                    used_marker_indcs = [used_marker_indcs, find(strcmp(this.raw_data(trial_no).marker_data(segment_no).marker_names, ...
                        this.sbj_marker_cluster_pos(cluster_no).marker_names{i}))]; %#ok<AGROW>
                else
                    % if the name is not found, remove the name in the stored trial transform data in WRAPS2 obj.
                    % and also the respectiive column of marker_static_pos
                    original_marker_names = this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).marker_names;
                    original_static_marker_pos = this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).marker_static_pos;
                    new_marker_names = original_marker_names(~strcmp(original_marker_names, original_marker_names{i}));
                    new_static_marker_pos = original_static_marker_pos(~strcmp(original_marker_names, original_marker_names{i}), :);
                    this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).marker_names = new_marker_names;
                    this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).marker_static_pos = new_static_marker_pos;
                end
            end
            
            % store the 3d matrix of marker pos from the trial in the same order as
            % the static cluster pos
            vicon_pos = this.raw_data(trial_no).marker_data(segment_no).marker_pos(:, :, used_marker_indcs);
            cluster_pos = this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).marker_static_pos;
            
            % Use Least Square Rigid Body Motion by SVD
            % (http://www.igl.ethz.ch/projects/ARAP/svd_rot.pdf)
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
            
            % store in the WRAPS2 instance
            this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).marker_pos = vicon_pos;
            this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).transforms_vicon2seg = T_v2s;
            disp(['Updated ', cluster_name, ' transformations in trial no. ', num2str(trial_no)])
            
        end
        
        % Calculate transformations between segment
        function calcTransformationBetweenClusters(this, trial_no, cluster_name_base, cluster_name_end)
            cluster_base_no = strcmp({this.sbj_WRAPS2(trial_no).trial_transform_data.cluster_name}, cluster_name_base);
            cluster_end_no = strcmp({this.sbj_WRAPS2(trial_no).trial_transform_data.cluster_name}, cluster_name_end);
            T_v2base = this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_base_no).transforms_vicon2seg;
            T_v2end = this.sbj_WRAPS2(trial_no).trial_transform_data(cluster_end_no).transforms_vicon2seg;  
            this.sbj_WRAPS2(trial_no).T_base2end =  this.multiplyTransforms(this.invTransformsMat(T_v2base), T_v2end);
            
             disp(['Updated ', cluster_name_base, ' to ', cluster_name_end, ' transformations in trial no. ', num2str(trial_no)])
        end
%         
        %% Calculate landmark positions
        
        function calcLandmarkPos(this, trial_no)
            this.sbj_anthro(trial_no).trial_name = this.raw_data(trial_no).marker_trial_name;
            this.sbj_anthro(trial_no).torso_landmark_names = {'UMBLC', 'R_ASIS', 'L_ASIS', 'R_PSIS', 'L_PSIS', 'R_HIP', 'L_HIP', 'SN', 'C7', 'XP', 'T8', 'VT', 'OC'};
            this.sbj_anthro(trial_no).landmark_pos = zeros(size(this.raw_data(trial_no).marker_data(1).marker_pos,1), 3, 12);
            this.sbj_anthro(trial_no).body_segment_transform = struct();
            this.sbj_anthro(trial_no).body_segment_transform(1).segment_name = 'Pelvis';
            this.sbj_anthro(trial_no).body_segment_transform(2).segment_name = 'Lumbar';
            this.sbj_anthro(trial_no).body_segment_transform(3).segment_name = 'Thorax';
            this.sbj_anthro(trial_no).body_segment_transform(4).segment_name = 'Head';
            % other body part without markers
            this.sbj_anthro(trial_no).body_segment_transform(5).segment_name = 'UpperArm_R';
            this.sbj_anthro(trial_no).body_segment_transform(6).segment_name = 'UpperArm_L';
            this.sbj_anthro(trial_no).body_segment_transform(7).segment_name = 'Forearm_R';
            this.sbj_anthro(trial_no).body_segment_transform(8).segment_name = 'Forearm_L';
            this.sbj_anthro(trial_no).body_segment_transform(9).segment_name = 'Hand_R';
            this.sbj_anthro(trial_no).body_segment_transform(10).segment_name = 'Hand_L';
            this.sbj_anthro(trial_no).body_segment_transform(11).segment_name = 'Thigh_R';
            this.sbj_anthro(trial_no).body_segment_transform(12).segment_name = 'Thigh_L';
            this.sbj_anthro(trial_no).body_segment_transform(13).segment_name = 'Shank_R';
            this.sbj_anthro(trial_no).body_segment_transform(14).segment_name = 'Shank_L';
            this.sbj_anthro(trial_no).body_segment_transform(15).segment_name = 'Foot_R';
            this.sbj_anthro(trial_no).body_segment_transform(16).segment_name = 'Foot_L';
            
            % get the pelvis landmark position
            pelvis_cluster_indx = strcmp({this.sbj_WRAPS2(trial_no).trial_transform_data.cluster_name}, 'Pelvis Brace');
            thorax_cluster_indx = strcmp({this.sbj_WRAPS2(trial_no).trial_transform_data.cluster_name}, 'Thorax Brace');
            T_v2pbrace = this.sbj_WRAPS2(trial_no).trial_transform_data(pelvis_cluster_indx).transforms_vicon2seg;
            T_v2tbrace = this.sbj_WRAPS2(trial_no).trial_transform_data(thorax_cluster_indx).transforms_vicon2seg;
            
            % get the five pelvis landmark based on the CAD
            for i = [1:5, 11]
                landmark_static_pos_indx = find(strcmp(this.sbj_landmark_config.WRAPS_landmark_names, this.sbj_anthro(trial_no).torso_landmark_names{i}));
                if i <= 5
                    T_pbrace2landmark = this.T_translate(this.sbj_landmark_config.WRAPS_landmark_static_pos(landmark_static_pos_indx, :));
                    for j = 1:size(T_v2pbrace, 3)
                        T = T_v2pbrace(:,:,j)*T_pbrace2landmark;
                        this.sbj_anthro(trial_no).landmark_pos(j,:,i) = T(1:3,4);
                    end
                else
                    T_tbrace2landmark = this.T_translate(this.sbj_landmark_config.WRAPS_landmark_static_pos(landmark_static_pos_indx, :));
                    for j = 1:size(T_v2tbrace, 3)
                        T = T_v2tbrace(:,:,j)*T_tbrace2landmark;
                        this.sbj_anthro(trial_no).landmark_pos(j,:,i) = T(1:3,4);
                    end
                end
            end
            
            % get raw marker positions
            %             % not using hip markers right now
            %             pelvis_raw_indx = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Pelvis');
            %             P_Hip_R_raw_indx = strcmp(this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_names, 'P_Hip_R');
            %             P_Hip_L_raw_indx = strcmp(this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_names, 'P_Hip_L');
            
            thorax_raw_indx = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Thorax');
            sternal_notch_raw_indx = strcmp(this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_names, 'T_SternalNotch');
            C7_raw_indx = strcmp(this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_names, 'T_C7');
            xyphoid_raw_indx = strcmp(this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_names, 'T_Xyphoid');
            
            head_raw_indx = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Head');
            vertex_raw_indx = strcmp(this.raw_data(trial_no).marker_data(head_raw_indx).marker_names, 'H_Vertex');
            forehead_raw_indx = strcmp(this.raw_data(trial_no).marker_data(head_raw_indx).marker_names, 'H_Forehead');
            
            %             % not using hip markers right now
            %             this.sbj_anthro(trial_no).landmark_pos(:,:,6) = this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_pos(:, :, P_Hip_R_raw_indx);
            %             this.sbj_anthro(trial_no).landmark_pos(:,:,7) = this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_pos(:, :, P_Hip_L_raw_indx);
            
            this.sbj_anthro(trial_no).landmark_pos(:,:,8) = this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_pos(:, :, sternal_notch_raw_indx);
            this.sbj_anthro(trial_no).landmark_pos(:,:,9) = this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_pos(:, :, C7_raw_indx);
            this.sbj_anthro(trial_no).landmark_pos(:,:,10) = this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_pos(:, :, xyphoid_raw_indx);
            this.sbj_anthro(trial_no).landmark_pos(:,:,12) = this.raw_data(trial_no).marker_data(head_raw_indx).marker_pos(:, :, vertex_raw_indx);
            this.sbj_anthro(trial_no).landmark_pos(:,:,13) = this.raw_data(trial_no).marker_data(head_raw_indx).marker_pos(:, :, forehead_raw_indx);
            
            disp(['Updated anthropometric landmark positions in trial no. ', num2str(trial_no)])
        end
        
        %% Calculate the segment transformation from obtianed landmark positions
        
        function calcSegmentTrans(this, trial_no)
            %% pelvis center frame using psis and asis points
            indx_r_asis = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'R_ASIS');
            indx_l_asis = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'L_ASIS');
            indx_r_psis = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'R_PSIS');
            indx_l_psis = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'L_PSIS');
            
            pos_r_asis = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_r_asis);
            pos_l_asis = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_l_asis);
            pos_r_psis = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_r_psis);
            pos_l_psis = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_l_psis);
            
            pos_m_psis = (pos_r_psis + pos_l_psis)/2;
            pos_m_asis = (pos_r_asis + pos_l_asis)/2;
            
            % i - unit vector parallel to line connection between <r_asis, l_asis>
            i_hat_pelvis = this.unitVecMat(pos_r_asis - pos_l_asis, 2);
            % k - unit vector perpendicular to the plane of <r_asis, l_asis, m_psis>
            k_hat_pelvis = this.unitVecMat(cross((pos_r_asis - pos_m_psis),(pos_l_asis - pos_m_psis), 2), 2);
            % j - unit vector from the cross product j = k x i
            j_hat_pelvis = cross(k_hat_pelvis, i_hat_pelvis, 2);
            
            origin_pelvis_distal = (pos_m_psis + pos_m_asis)/2; % this is technically distal pelvis frame since the proximal is align with the vicon frame
            
            pelvis_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Pelvis'));
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_proxim = ...
                repmat(eye(4), 1, 1, length(origin_pelvis_distal));
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal = ...
                this.createTransforms(origin_pelvis_distal, i_hat_pelvis, j_hat_pelvis, k_hat_pelvis);
            
            %% find the estimate hip position from asis and psis pos
            % R. B. Davis, S. Õunpuu, D. Tyburski, and J. R. Gage,
            % “A gait analysis data collection and reduction technique,”
            % Hum. Mov. Sci., vol. 10, no. 5, pp. 575–587, Oct. 1991.
            
            % creat a the frame based on the reference
            T_v2masis = this.createTransforms(pos_m_asis, i_hat_pelvis, j_hat_pelvis, k_hat_pelvis);
            
            dasis_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'dasis');
            asis2hip_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'asis2hip_horiz');
            thigh_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'thigh');
            shank_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'shank');
            
            dasis = this.sbj_anthro_measurement.measurement_table(dasis_indx).length_mm/1000; % meters
            asis2hip = this.sbj_anthro_measurement.measurement_table(asis2hip_indx).length_mm/1000; % meters
            L_leg = this.sbj_anthro_measurement.measurement_table(thigh_indx).length_mm/1000 + ...
                this.sbj_anthro_measurement.measurement_table(shank_indx).length_mm/1000;
            
            theta = 28.4*pi/180; %rad
            beta = 18*pi/180; % rad
            r_marker = 0.014/2; % m
            
            C = 0.115*L_leg - 0.0153;
            
            trans_vec_R = [ abs(C*sin(theta) - dasis/2),...
                - abs(-(asis2hip + r_marker)*cos(beta) + C*cos(theta)*sin(beta)),...
                - abs(-(asis2hip + r_marker)*sin(beta) - C*cos(theta)*cos(beta))];
            
            trans_vec_L = [ -abs((C*sin(theta) - dasis/2)),...
                - abs(-(asis2hip + r_marker)*cos(beta) + C*cos(theta)*sin(beta)),...
                - abs(-(asis2hip + r_marker)*sin(beta) - C*cos(theta)*cos(beta))];
            
            % convert back to mm and construct the vector matrices
            T_masis2hip_r = this.createTransformsTranslation(1000*repmat(trans_vec_R, length(pos_m_asis), 1));
            T_masis2hip_l = this.createTransformsTranslation(1000*repmat(trans_vec_L, length(pos_m_asis), 1));
            
            T_v2hip_r = this.multiplyTransforms(T_v2masis, T_masis2hip_r);
            T_v2hip_l = this.multiplyTransforms(T_v2masis, T_masis2hip_l);
            
            % store as the landmark pos
            pos_r_hip_joint = reshape(T_v2hip_r(1:3, 4, :), 3, [])';
            pos_l_hip_joint = reshape(T_v2hip_l(1:3, 4, :), 3, [])';
            this.sbj_anthro(trial_no).landmark_pos(:,:,6) = pos_r_hip_joint;
            this.sbj_anthro(trial_no).landmark_pos(:,:,7) = pos_l_hip_joint;
            
            % store the hip joint frames
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal_aux_r = T_v2hip_r;
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal_aux_l = T_v2hip_l;
            
            %% lumber using 'origin_pelvis', 'XP', 'T8'
            indx_xp = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'XP');
            indx_t8 = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'T8');
            
            pos_xp = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_xp);
            pos_t8 = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_t8);
            
            origin_lumbar_proxim = origin_pelvis_distal;
            origin_lumbar_distal = (pos_xp + pos_t8)/2;
            
            % k - unit vector parallel to line connection between <origin_lumbar_proxim, origin_lumbar_distal>
            k_hat_lumbar = this.unitVecMat(origin_lumbar_distal - origin_lumbar_proxim, 2);
            % i - unit vector perpendicular to the plane of <origin_lumbar_proxim, l_asis, m_psis>
            i_hat_lumbar = this.unitVecMat(cross((pos_xp - origin_lumbar_proxim),(pos_t8 - origin_lumbar_proxim), 2), 2);
            % j - unit vector from the cross product j = k x i
            j_hat_lumbar = cross(k_hat_lumbar, i_hat_lumbar, 2);
            
            lumber_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Lumbar'));
            this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_proxim = ...
                this.createTransforms(origin_lumbar_proxim, i_hat_lumbar, j_hat_lumbar, k_hat_lumbar);
            this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_distal = ...
                this.createTransforms(origin_lumbar_distal, i_hat_lumbar, j_hat_lumbar, k_hat_lumbar);
            
            %% thorax using 'SN', 'C7', 'XP', 'T8'
            indx_sn = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'SN');
            indx_c7 = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'C7');
            
            pos_sn = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_sn);
            pos_c7 = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_c7);
            
            origin_thorax_proxim = origin_lumbar_distal;
            origin_thorax_distal = (pos_sn + pos_c7)/2;
            
            % k - unit vector parallel to line connection between <origin_thorax_proxim, origin_thorax_distal>
            k_hat_thorax = this.unitVecMat(origin_thorax_distal - origin_thorax_proxim, 2);
            % i - unit vector perpendicular to the plane of <origin_lumbar_proxim, l_asis, m_psis>
            i_hat_thorax = this.unitVecMat(cross((pos_sn - origin_thorax_proxim),(pos_c7 - origin_thorax_proxim), 2), 2);
            % j - unit vector from the cross product j = k x i
            j_hat_thorax = cross(k_hat_thorax, i_hat_thorax, 2);
            
            thorax_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thorax'));
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_proxim = ...
                this.createTransforms(origin_thorax_proxim, i_hat_thorax, j_hat_thorax, k_hat_thorax);
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal = ...
                this.createTransforms(origin_thorax_distal, i_hat_thorax, j_hat_thorax, k_hat_thorax);
            
            %% head using 'SN', 'C7', 'VT'
            indx_vt = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'VT');
            indx_oc = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'OC');
            
            pos_vt = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_vt);
            pos_fh = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_oc);
            
            origin_head_proxim = origin_thorax_distal;
            origin_head_distal = pos_vt;
            
            % k - unit vector parallel to line connection between <origin_head_proxim, origin_head_distal>
            k_hat_head = this.unitVecMat(origin_head_distal - origin_head_proxim, 2);
            % i - unit vector perpendicular to the plane of <origin_lumbar_proxim, l_asis, m_psis>
            i_hat_head = this.unitVecMat(cross((pos_fh - origin_head_proxim),(origin_head_distal - origin_head_proxim), 2), 2);
            % j - unit vector from the cross product j = k x i
            j_hat_head = cross(k_hat_head, i_hat_head, 2);
            
            head_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Head'));
            this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_proxim = ...
                this.createTransforms(origin_head_proxim, i_hat_head, j_hat_head, k_hat_head);
            this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_distal = ...
                this.createTransforms(origin_head_distal, i_hat_head, j_hat_head, k_hat_head);
            
            %% upper Arms transformation using the thorax coordinate frame
            upper_arm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_R'));
            upper_arm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_L'));
            
            shoulder_width_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'shoulder_width');
            z_sn2shouder_center_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'sternal2shoulder_center');
            upper_arm_length_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'upper_arm');
            
            shoulder_width = this.sbj_anthro_measurement.measurement_table(shoulder_width_indx).length_mm;
            z_sn2shouder_center = this.sbj_anthro_measurement.measurement_table(z_sn2shouder_center_indx).length_mm;
            l_upper_arm = this.sbj_anthro_measurement.measurement_table(upper_arm_length_indx).length_mm;
            
            % get the calculated thorax distal frame
            T_v2thorax_distal = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal;
            
            trans_thorax2shdr_R_vec_mat = repmat([shoulder_width/2, 0, -z_sn2shouder_center], length(T_v2thorax_distal), 1);
            trans_thorax2shdr_L_vec_mat = repmat([-shoulder_width/2, 0, -z_sn2shouder_center], length(T_v2thorax_distal), 1);
            
            T_thorax_distal2thorax_shdr_joint_R = this.createTransformsTranslation(trans_thorax2shdr_R_vec_mat);
            T_thorax_distal2thorax_shdr_joint_L = this.createTransformsTranslation(trans_thorax2shdr_L_vec_mat);
            
            T_v2thorax_shdr_R = this.multiplyTransforms(T_v2thorax_distal, T_thorax_distal2thorax_shdr_joint_R);
            T_v2thorax_shdr_L = this.multiplyTransforms(T_v2thorax_distal, T_thorax_distal2thorax_shdr_joint_L);
            
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal_aux_r = T_v2thorax_shdr_R;
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal_aux_l = T_v2thorax_shdr_L;
            
            % specify the shoulder joint rotation
            trial_length = length(T_v2thorax_shdr_R);
            T_thorax_shdr_R2uarm_R_proxim  = ...
                this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), zeros(trial_length, 1)); % no rotation
            T_thorax_shdr_L2uarm_L_proxim  = ...
                this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), zeros(trial_length, 1)); % try rotation around x by +pi/2
            
            T_v2uarm_R_proxim = this.multiplyTransforms(T_v2thorax_shdr_R, T_thorax_shdr_R2uarm_R_proxim);
            T_v2uarm_L_proxim = this.multiplyTransforms(T_v2thorax_shdr_L, T_thorax_shdr_L2uarm_L_proxim);
            
            T_uarm_R_proxim2uarm_R_distal = this.createTransformsTranslation(repmat([0, 0, -l_upper_arm], trial_length, 1));
            T_uarm_L_proxim2uarm_L_distal = this.createTransformsTranslation(repmat([0, 0, -l_upper_arm], trial_length, 1));
            
            T_v2uarm_R_distal = this.multiplyTransforms(T_v2uarm_R_proxim, T_uarm_R_proxim2uarm_R_distal);
            T_v2uarm_L_distal = this.multiplyTransforms(T_v2uarm_L_proxim, T_uarm_L_proxim2uarm_L_distal);
            
            % store proximal and distal transforms into the structure
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).T_v2seg_proxim = T_v2uarm_R_proxim ;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).T_v2seg_distal = T_v2uarm_R_distal;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).T_v2seg_proxim = T_v2uarm_L_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).T_v2seg_distal = T_v2uarm_L_distal;
            
            %% forearm
            forearm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_R'));
            forearm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_L'));
            
            forearm_length_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'forearm');
            l_forearm = this.sbj_anthro_measurement.measurement_table(forearm_length_indx).length_mm;
            
            % elbow rotation (restricted in x axis)
            T_uarm_R_distal2farm_R_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), zeros(trial_length, 1)); % no rotation
            T_uarm_L_distal2farm_L_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), zeros(trial_length, 1)); % no rotation
            
            T_v2farm_R_proxim = this.multiplyTransforms(T_v2uarm_R_distal, T_uarm_R_distal2farm_R_proxim);
            T_v2farm_L_proxim = this.multiplyTransforms(T_v2uarm_L_distal, T_uarm_L_distal2farm_L_proxim);
            
            T_farm_R_proxim2farm_R_distal = this.createTransformsTranslation(repmat([0, 0, -l_forearm], trial_length, 1));
            T_farm_L_proxim2farm_L_distal = this.createTransformsTranslation(repmat([0, 0, -l_forearm], trial_length, 1));
            
            T_v2farm_R_distal = this.multiplyTransforms(T_v2farm_R_proxim, T_farm_R_proxim2farm_R_distal);
            T_v2farm_L_distal = this.multiplyTransforms(T_v2farm_L_proxim, T_farm_L_proxim2farm_L_distal);
            
            % store proximal and distal transforms into the structure
            this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2seg_proxim = T_v2farm_R_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2seg_distal = T_v2farm_R_distal;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2seg_proxim = T_v2farm_L_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2seg_distal = T_v2farm_L_distal;
            
            %% hand
            hand_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_R'));
            hand_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_L'));
            
            hand_length_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'hand');
            l_hand = this.sbj_anthro_measurement.measurement_table(hand_length_indx).length_mm;
            
            % wrist rotation (restricted in y axis)
            T_hand_R_distal2hand_R_proxim = this.createRotationAxisAngle(repmat([0 1 0], trial_length, 1), zeros(trial_length, 1));
            T_hand_L_distal2hand_L_proxim = this.createRotationAxisAngle(repmat([0 1 0], trial_length, 1), zeros(trial_length, 1));
            
            T_v2hand_R_proxim = this.multiplyTransforms(T_v2farm_R_distal, T_hand_R_distal2hand_R_proxim);
            T_v2hand_L_proxim = this.multiplyTransforms(T_v2farm_L_distal, T_hand_L_distal2hand_L_proxim);
            
            T_hand_R_proxim2hand_R_distal = this.createTransformsTranslation(repmat([0, 0, -l_hand], trial_length, 1));
            T_hand_L_proxim2hand_L_distal = this.createTransformsTranslation(repmat([0, 0, -l_hand], trial_length, 1));
            
            T_v2hand_R_distal = this.multiplyTransforms(T_v2hand_R_proxim, T_hand_R_proxim2hand_R_distal);
            T_v2hand_L_distal = this.multiplyTransforms(T_v2hand_L_proxim, T_hand_L_proxim2hand_L_distal);
            
            % store proximal and distal transforms into the structure
            this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2seg_proxim = T_v2hand_R_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2seg_distal = T_v2hand_R_distal;
            this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2seg_proxim = T_v2hand_L_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2seg_distal = T_v2hand_L_distal;
            
            %% thigh
            thigh_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_R'));
            thigh_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_L'));
            
            thigh_length_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'thigh');
            l_thigh = this.sbj_anthro_measurement.measurement_table(thigh_length_indx).length_mm;
            
            % hip rotation (sphrical joint but limit in x rotation for now)
            T_hip_r2thigh_r_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), pi/2*ones(trial_length, 1));
            T_hip_l2thigh_l_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), pi/2*ones(trial_length, 1));
            
            T_v2thigh_r_proxim = this.multiplyTransforms(T_v2hip_r, T_hip_r2thigh_r_proxim);
            T_v2thigh_l_proxim = this.multiplyTransforms(T_v2hip_l, T_hip_l2thigh_l_proxim);
            
            T_thigh_r_proxim2thigh_r_distal = this.createTransformsTranslation(repmat([0, 0, -l_thigh], trial_length, 1));
            T_thigh_l_proxim2thigh_l_proxim = this.createTransformsTranslation(repmat([0, 0, -l_thigh], trial_length, 1));
            
            T_v2thigh_r_distal = this.multiplyTransforms(T_v2thigh_r_proxim, T_thigh_r_proxim2thigh_r_distal);
            T_v2thigh_l_distal = this.multiplyTransforms(T_v2thigh_l_proxim, T_thigh_l_proxim2thigh_l_proxim);
            
            % store proximal and distal transforms into the structure
            this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_proxim = T_v2thigh_r_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_distal = T_v2thigh_r_distal;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_proxim = T_v2thigh_l_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_distal = T_v2thigh_l_distal;
            
            
            %% shank
            shank_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_R'));
            shank_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_L'));
            
            shank_length_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'shank');
            l_shank = this.sbj_anthro_measurement.measurement_table(shank_length_indx).length_mm;
            
            % hip rotation (sphrical joint but limit in x rotation for now)
            T_thigh_r_distal2shank_r_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), -pi/2*ones(trial_length, 1));
            T_thigh_l_distal2shank_l_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), -pi/2*ones(trial_length, 1));
            
            T_v2shank_r_proxim = this.multiplyTransforms(T_v2thigh_r_distal, T_thigh_r_distal2shank_r_proxim);
            T_v2shank_l_proxim = this.multiplyTransforms(T_v2thigh_l_distal, T_thigh_l_distal2shank_l_proxim);
            
            T_shank_r_proxim2shank_r_distal = this.createTransformsTranslation(repmat([0, 0, -l_shank], trial_length, 1));
            T_shank_l_proxim2shank_l_proxim = this.createTransformsTranslation(repmat([0, 0, -l_shank], trial_length, 1));
            
            T_v2shank_r_distal = this.multiplyTransforms(T_v2shank_r_proxim, T_shank_r_proxim2shank_r_distal);
            T_v2shank_l_distal = this.multiplyTransforms(T_v2shank_l_proxim, T_shank_l_proxim2shank_l_proxim);
            
            % store proximal and distal transforms into the structure
            this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_proxim = T_v2shank_r_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_distal = T_v2shank_r_distal;
            this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_proxim = T_v2shank_l_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_distal = T_v2shank_l_distal;
            
            
            
            %% ankle
            foot_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_R'));
            foot_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_L'));
            
            foot_length_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'foot');
            heel2ankle_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'heel2ankle');
            ankle2ground_indx = strcmp({this.sbj_anthro_measurement.measurement_table.name}, 'ankle2ground');
            
            l_foot = this.sbj_anthro_measurement.measurement_table(foot_length_indx).length_mm;
            l_heel2ankle = this.sbj_anthro_measurement.measurement_table(heel2ankle_indx).length_mm;
            l_ankle2ground = this.sbj_anthro_measurement.measurement_table(ankle2ground_indx).length_mm;
            
            % ankle rotation (sphrical joint but limit in x rotation for now)
            T_shank_r_distal2foot_r_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), zeros(trial_length, 1));
            T_shank_l_distal2foot_l_proxim = this.createRotationAxisAngle(repmat([1 0 0], trial_length, 1), zeros(trial_length, 1));
            
            T_v2foot_r_proxim = this.multiplyTransforms(T_v2shank_r_distal, T_shank_r_distal2foot_r_proxim);
            T_v2foot_l_proxim = this.multiplyTransforms(T_v2shank_l_distal, T_shank_l_distal2foot_l_proxim);
            
            T_foot_r_proxim2foot_r_distal = this.createTransformsTranslation(repmat([0, l_foot - l_heel2ankle, -l_ankle2ground], trial_length, 1));
            T_foot_l_proxim2foot_l_distal = this.createTransformsTranslation(repmat([0, l_foot - l_heel2ankle, -l_ankle2ground], trial_length, 1));
            T_foot_r_proxim2heel_r = this.createTransformsTranslation(repmat([0, - l_heel2ankle, -l_ankle2ground], trial_length, 1));
            T_foot_l_proxim2heel_l = this.createTransformsTranslation(repmat([0, - l_heel2ankle, -l_ankle2ground], trial_length, 1));
            
            T_v2foot_r_distal = this.multiplyTransforms(T_v2foot_r_proxim, T_foot_r_proxim2foot_r_distal);
            T_v2foot_l_distal = this.multiplyTransforms(T_v2foot_l_proxim, T_foot_l_proxim2foot_l_distal);
            T_v2heel_r = this.multiplyTransforms(T_v2foot_r_proxim, T_foot_r_proxim2heel_r);
            T_v2heel_l = this.multiplyTransforms(T_v2foot_l_proxim, T_foot_l_proxim2heel_l);
            
            % store proximal and distal transforms into the structure
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_proxim = T_v2foot_r_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal = T_v2foot_r_distal;
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal_aux_r =T_v2heel_r;
            
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_proxim = T_v2foot_l_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal = T_v2foot_l_distal;
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal_aux_l = T_v2heel_l;
            
            %% Correct the lower limb configuration using inverse kinematics
            % store the first proximal foot frame transformation and use
            % it as the starting point of the inverse kinematics. This
            % idea is based on the fixed ankle assumption
            T_v2foot_r_fix = repmat(T_v2foot_r_proxim(:, :, 1), 1, 1, trial_length);
            T_v2foot_l_fix = repmat(T_v2foot_l_proxim(:, :, 1), 1, 1, trial_length);
            
            [T_hip_r2thigh_r_prox_new, T_thigh_r_dist2shank_r_prox_new, T_shank_r_dist2foot_r_prox_new] = ...
                this.inverseKinLowerLimbs(T_v2hip_r, T_v2foot_r_fix, l_thigh, l_shank);
            [T_hip_l2thigh_l_prox_new, T_thigh_l_dist2shank_l_prox_new, T_shank_L_dist2foot_l_prox_new] = ...
                this.inverseKinLowerLimbs(T_v2hip_l, T_v2foot_l_fix, l_thigh, l_shank);
            
            % redo the forward calculation
            T_v2thigh_r_proxim = this.multiplyTransforms(T_v2hip_r, T_hip_r2thigh_r_prox_new);
            T_v2thigh_l_proxim = this.multiplyTransforms(T_v2hip_l, T_hip_l2thigh_l_prox_new);
            T_v2thigh_r_distal = this.multiplyTransforms(T_v2thigh_r_proxim, T_thigh_r_proxim2thigh_r_distal);
            T_v2thigh_l_distal = this.multiplyTransforms(T_v2thigh_l_proxim, T_thigh_l_proxim2thigh_l_proxim);
            
            T_v2shank_r_proxim = this.multiplyTransforms(T_v2thigh_r_distal, T_thigh_r_dist2shank_r_prox_new);
            T_v2shank_l_proxim = this.multiplyTransforms(T_v2thigh_l_distal, T_thigh_l_dist2shank_l_prox_new);
            T_v2shank_r_distal = this.multiplyTransforms(T_v2shank_r_proxim, T_shank_r_proxim2shank_r_distal);
            T_v2shank_l_distal = this.multiplyTransforms(T_v2shank_l_proxim, T_shank_l_proxim2shank_l_proxim);
            
            T_v2foot_r_proxim = this.multiplyTransforms(T_v2shank_r_distal, T_shank_r_dist2foot_r_prox_new);
            T_v2foot_l_proxim = this.multiplyTransforms(T_v2shank_l_distal, T_shank_L_dist2foot_l_prox_new);
            T_v2foot_r_distal = this.multiplyTransforms(T_v2foot_r_proxim, T_foot_r_proxim2foot_r_distal);
            T_v2foot_l_distal = this.multiplyTransforms(T_v2foot_l_proxim, T_foot_l_proxim2foot_l_distal);
            T_v2heel_r = this.multiplyTransforms(T_v2foot_r_proxim, T_foot_r_proxim2heel_r);
            T_v2heel_l = this.multiplyTransforms(T_v2foot_l_proxim, T_foot_l_proxim2heel_l);
            
            % store proximal and distal transforms into the structure
            % thigh
            this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_proxim = T_v2thigh_r_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_distal = T_v2thigh_r_distal;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_proxim = T_v2thigh_l_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_distal = T_v2thigh_l_distal;
            
            % shank
            this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_proxim = T_v2shank_r_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_distal = T_v2shank_r_distal;
            this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_proxim = T_v2shank_l_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_distal = T_v2shank_l_distal;
            
            % ankle
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_proxim = T_v2foot_r_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal = T_v2foot_r_distal;
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal_aux_r =T_v2heel_r;
            
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_proxim = T_v2foot_l_proxim;
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal = T_v2foot_l_distal;
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal_aux_l = T_v2heel_l;
            
            disp(['Updated anthropometric segment transformations in trial no. ', num2str(trial_no)])
        end
        
        %% Inverse kinematic for lower limb positions
        function [T_hip2thigh_prox, T_thigh_dist2shank_prox, T_shank_dis2foot_prox] = inverseKinLowerLimbs(this, T_v2hip, T_v2foot_fix, l_thigh, l_shank)
            % inverseKinLowerLimbs: calculate reverse joint angles of the
            %   from the position of the moving hip wrt proximal ankle frame
            %   and calculate the transformation matrices from the hip to the
            %   ankle joint in the forward manner
            % the kinematic chains from ankle to hip is Uyx->Rx->S(hip)
            
            pos_hip_wrt_foot_fix = this.calcPosInNewFrame(T_v2foot_fix, reshape(T_v2hip(1:3, 4, :), 3, [])');
            
            L_thigh = l_thigh*ones(length(pos_hip_wrt_foot_fix), 1);
            L_shank = l_shank*ones(length(pos_hip_wrt_foot_fix), 1);
            %
            D = ((sum(pos_hip_wrt_foot_fix.^2, 2) - L_thigh.^2 - L_shank.^2))./(2*L_thigh.*L_shank);
            q_invs_knee = atan2(sqrt(1 - D.^2), D);
            q_invs_ankle_x = abs(atan2(pos_hip_wrt_foot_fix(:, 2), sqrt(pos_hip_wrt_foot_fix(:, 1).^2 + pos_hip_wrt_foot_fix(:, 3).^2)))...
                - atan2(L_thigh.*sin(q_invs_knee), L_shank + L_thigh.*cos(q_invs_knee)); % assuming that the second term is the small angle (elbow down)
            q_invs_ankle_y = atan2(pos_hip_wrt_foot_fix(:, 1), pos_hip_wrt_foot_fix(:, 3));
            
            T_invs_ankle_y = this.createRotationAxisAngle(repmat([0 1 0], length(pos_hip_wrt_foot_fix), 1), q_invs_ankle_y);
            T_invs_ankle_x = this.createRotationAxisAngle(repmat([1 0 0], length(pos_hip_wrt_foot_fix), 1), q_invs_ankle_x);
            
            T_foot_prox2shank_dist = this.multiplyTransforms(T_invs_ankle_y, T_invs_ankle_x);
            T_shank_dist2shank_prox = this.createTransformsTranslation(repmat([0, 0, l_shank], length(pos_hip_wrt_foot_fix), 1));
            T_shank_prox2thigh_dist = this.createRotationAxisAngle(repmat([1 0 0], length(pos_hip_wrt_foot_fix), 1), q_invs_knee);
            T_thigh_dist2thigh_prox =  this.createTransformsTranslation(repmat([0, 0, l_thigh], length(pos_hip_wrt_foot_fix), 1));
            
            T_v2shank_dist = this.multiplyTransforms(T_v2foot_fix, T_foot_prox2shank_dist);
            T_v2shank_prox = this.multiplyTransforms(T_v2shank_dist, T_shank_dist2shank_prox);
            T_v2thigh_dist = this.multiplyTransforms(T_v2shank_prox, T_shank_prox2thigh_dist);
            T_v2thigh_prox = this.multiplyTransforms(T_v2thigh_dist, T_thigh_dist2thigh_prox);
            
            % return solutions
            T_shank_dis2foot_prox = this.invTransformsMat(T_foot_prox2shank_dist); % forward ankle joint
            T_thigh_dist2shank_prox = this.invTransformsMat(T_shank_prox2thigh_dist); % forward ankle joint
            T_hip2thigh_prox = this.multiplyTransforms(this.invTransformsMat(T_v2hip), T_v2thigh_prox);
        end
        
        %% Calculate segment intertia
        function calcSegmentInertia(this, trial_no)
            % find landmarks projected on the local segment frame
            anthro_pelvis_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Pelvis');
            anthro_lumbar_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Lumbar');
            anthro_thorax_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Thorax');
            anthro_head_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Head');
            anthro_uarm_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'UpperArm');
            anthro_farm_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Forearm');
            anthro_hand_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Hand');
            anthro_thigh_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Thigh');
            anthro_shank_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Shank');
            anthro_foot_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Foot');
            
            cm_pos_ratio_pelvis_umblc2hip = this.sbj_anthro_table(anthro_pelvis_table_indx).cm_pos_ratio  ; % male  0.6115
            cm_pos_ratio_lumbar_xp2umblc = this.sbj_anthro_table(anthro_lumbar_table_indx).cm_pos_ratio  ; % 0.4502
            cm_pos_ratio_thorax_sn2xp = this.sbj_anthro_table(anthro_thorax_table_indx).cm_pos_ratio  ; % 0.2999
            cm_pos_ratio_head_vt2c7 = this.sbj_anthro_table(anthro_head_table_indx).cm_pos_ratio  ; % 0.5002
            cm_pos_ratio_uarm_prox2dist = this.sbj_anthro_table(anthro_uarm_table_indx).cm_pos_ratio  ; % 0.5772
            cm_pos_ratio_farm_prox2dist = this.sbj_anthro_table(anthro_farm_table_indx).cm_pos_ratio  ; % 0.4574
            cm_pos_ratio_hand_prox2dist = this.sbj_anthro_table(anthro_hand_table_indx).cm_pos_ratio  ; % 0.3624
            cm_pos_ratio_thigh_prox2dist = this.sbj_anthro_table(anthro_thigh_table_indx).cm_pos_ratio  ; % 0.4095
            cm_pos_ratio_shank_prox2dist = this.sbj_anthro_table(anthro_shank_table_indx).cm_pos_ratio  ; % 0.4459
            cm_pos_ratio_foot_heel2toe = this.sbj_anthro_table(anthro_foot_table_indx).cm_pos_ratio  ; % 0.4415
            
            mass_ratio_pelvis = this.sbj_anthro_table(anthro_pelvis_table_indx).mass_ratio; % male 0.1117
            mass_ratio_lumbar = this.sbj_anthro_table(anthro_lumbar_table_indx).mass_ratio; % 0.1633
            mass_ratio_thorax = this.sbj_anthro_table(anthro_thorax_table_indx).mass_ratio; % 0.1596
            mass_ratio_head = this.sbj_anthro_table(anthro_head_table_indx).mass_ratio; % 0.0694
            mass_ratio_uarm = this.sbj_anthro_table(anthro_uarm_table_indx).mass_ratio; % 0.0271
            mass_ratio_farm = this.sbj_anthro_table(anthro_farm_table_indx).mass_ratio; % 0.0162
            mass_ratio_hand = this.sbj_anthro_table(anthro_hand_table_indx).mass_ratio; % 0.0061
            mass_ratio_thigh = this.sbj_anthro_table(anthro_thigh_table_indx).mass_ratio; % 0.1416
            mass_ratio_shank = this.sbj_anthro_table(anthro_shank_table_indx).mass_ratio; % 0.0433
            mass_ratio_foot = this.sbj_anthro_table(anthro_foot_table_indx).mass_ratio; % 0.0137
            
            sbj_mass = this.sbj_anthro_measurement.weight_kg;
            
            indx_r_hip = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'R_HIP');
            indx_l_hip = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'L_HIP');
            indx_umblc = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'UMBLC');
            indx_xp = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'XP');
            indx_c7 = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'C7');
            indx_sn = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'SN');
            indx_vt = strcmp(this.sbj_anthro(trial_no).torso_landmark_names, 'VT');
            
            pos_umblc = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_umblc);
            pos_xp = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_xp);
            pos_r_hip = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_r_hip);
            pos_l_hip = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_l_hip);
            pos_c7 = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_c7);
            pos_sn = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_sn);
            pos_vt = this.sbj_anthro(trial_no).landmark_pos(:, :, indx_vt);
            pos_m_hip = (pos_r_hip+pos_l_hip)/2;
            
            pelvis_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Pelvis'));
            lumber_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Lumbar'));
            thorax_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thorax'));
            head_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Head'));
            upper_arm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_R'));
            upper_arm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_L'));
            forearm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_R'));
            forearm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_L'));
            hand_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_R'));
            hand_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_L'));
            thigh_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_R'));
            thigh_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_L'));
            shank_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_R'));
            shank_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_L'));
            foot_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_R'));
            foot_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_L'));
            
            T_v2pelv_dist = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal;
            T_v2lumb_dist = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_distal;
            T_v2thor_dist = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal;
            T_v2head_dist = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_distal;
            
            [T_v2cm_pelvis, pos_cm_pelvis] = this.calcCMTransforms(T_v2pelv_dist, pos_umblc, pos_m_hip, cm_pos_ratio_pelvis_umblc2hip);
            [T_v2cm_lumbar, pos_cm_lumbar] = this.calcCMTransforms(T_v2lumb_dist, pos_xp, pos_umblc, cm_pos_ratio_lumbar_xp2umblc);
            [T_v2cm_thorax, pos_cm_thorax] = this.calcCMTransforms(T_v2thor_dist, pos_sn, pos_xp, cm_pos_ratio_thorax_sn2xp);
            [T_v2cm_head, pos_cm_head] = this.calcCMTransforms(T_v2head_dist, pos_vt, pos_c7, cm_pos_ratio_head_vt2c7);
            
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2cm_seg = T_v2cm_pelvis;
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).pos_cm_seg = pos_cm_pelvis;
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).seg_mass = sbj_mass*mass_ratio_pelvis;
            
            this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2cm_seg = T_v2cm_lumbar;
            this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).pos_cm_seg = pos_cm_lumbar;
            this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).seg_mass = sbj_mass*mass_ratio_lumbar;
            
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2cm_seg = T_v2cm_thorax;
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).pos_cm_seg = pos_cm_thorax;
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).seg_mass =  sbj_mass*mass_ratio_thorax;
            
            this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2cm_seg = T_v2cm_head;
            this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).pos_cm_seg = pos_cm_head;
            this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).seg_mass =  sbj_mass*mass_ratio_head;
            
            %% calculate limb CoM transformation
            T_v2uarm_R_proxim = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).T_v2seg_proxim;
            T_v2uarm_L_proxim = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).T_v2seg_proxim;
            T_v2uarm_R_distal = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).T_v2seg_distal;
            T_v2uarm_L_distal = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).T_v2seg_distal;
            
            T_v2farm_R_proxim = this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2seg_proxim;
            T_v2farm_L_proxim = this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2seg_proxim;
            T_v2farm_R_distal = this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2seg_distal;
            T_v2farm_L_distal = this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2seg_distal;
            
            T_v2hand_R_proxim = this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2seg_proxim;
            T_v2hand_L_proxim = this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2seg_proxim;
            T_v2hand_R_distal = this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2seg_distal;
            T_v2hand_L_distal = this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2seg_distal;
            
            
            % proximal joints are superior frames this case
            % upper limbs
            [T_v2cm_uarm_R, pos_cm_uarm_R] = this.calcLimbCMTranforms(T_v2uarm_R_proxim, T_v2uarm_R_distal, cm_pos_ratio_uarm_prox2dist);
            [T_v2cm_uarm_L, pos_cm_uarm_L] = this.calcLimbCMTranforms(T_v2uarm_L_proxim, T_v2uarm_L_distal, cm_pos_ratio_uarm_prox2dist);
            [T_v2cm_farm_R, pos_cm_farm_R] = this.calcLimbCMTranforms(T_v2farm_R_proxim, T_v2farm_R_distal, cm_pos_ratio_farm_prox2dist);
            [T_v2cm_farm_L, pos_cm_farm_L] = this.calcLimbCMTranforms(T_v2farm_L_proxim, T_v2farm_L_distal, cm_pos_ratio_farm_prox2dist);
            [T_v2cm_hand_R, pos_cm_hand_R] = this.calcLimbCMTranforms(T_v2hand_R_proxim, T_v2hand_R_distal, cm_pos_ratio_hand_prox2dist);
            [T_v2cm_hand_L, pos_cm_hand_L] = this.calcLimbCMTranforms(T_v2hand_L_proxim, T_v2hand_L_distal, cm_pos_ratio_hand_prox2dist);
            
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).T_v2cm_seg = T_v2cm_uarm_R;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).pos_cm_seg = pos_cm_uarm_R;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).seg_mass = sbj_mass*mass_ratio_uarm;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).T_v2cm_seg = T_v2cm_uarm_L;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).pos_cm_seg = pos_cm_uarm_L;
            this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).seg_mass = sbj_mass*mass_ratio_uarm;
            
            this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2cm_seg = T_v2cm_farm_R;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).pos_cm_seg = pos_cm_farm_R;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).seg_mass = sbj_mass*mass_ratio_farm;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2cm_seg = T_v2cm_farm_L;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).pos_cm_seg = pos_cm_farm_L;
            this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).seg_mass = sbj_mass*mass_ratio_farm;
            
            this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2cm_seg = T_v2cm_hand_R;
            this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).pos_cm_seg = pos_cm_hand_R;
            this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).seg_mass = sbj_mass*mass_ratio_hand;
            this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2cm_seg = T_v2cm_hand_L;
            this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).pos_cm_seg = pos_cm_hand_L;
            this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).seg_mass = sbj_mass*mass_ratio_hand;
            
            % get transformation lower limbs
            T_v2thigh_r_proxim = this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_proxim;
            T_v2thigh_r_distal = this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_distal;
            T_v2thigh_l_proxim = this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_proxim;
            T_v2thigh_l_distal = this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_distal;
            
            T_v2shank_r_proxim = this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_proxim;
            T_v2shank_r_distal = this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_distal;
            T_v2shank_l_proxim = this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_proxim;
            T_v2shank_l_distal = this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_distal;
            
            T_v2foot_r_proxim = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_proxim;
            T_v2foot_r_distal = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal;
            T_v2heel_r = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal_aux_r;
            
            T_v2foot_l_proxim = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_proxim;
            T_v2foot_l_distal = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal;
            T_v2heel_l = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal_aux_l;
                     
            % CoM of thighs and shanks
            [T_v2cm_thigh_R, pos_cm_thigh_R] = this.calcLimbCMTranforms(T_v2thigh_r_proxim, T_v2thigh_r_distal, cm_pos_ratio_thigh_prox2dist);
            [T_v2cm_thigh_L, pos_cm_thigh_L] = this.calcLimbCMTranforms(T_v2thigh_l_proxim, T_v2thigh_l_distal, cm_pos_ratio_thigh_prox2dist);
            [T_v2cm_shank_R, pos_cm_shank_R] = this.calcLimbCMTranforms(T_v2shank_r_proxim, T_v2shank_r_distal, cm_pos_ratio_shank_prox2dist);
            [T_v2cm_shank_L, pos_cm_shank_L] = this.calcLimbCMTranforms(T_v2shank_l_proxim, T_v2shank_l_distal, cm_pos_ratio_shank_prox2dist);            
            
            [T_v2cm_foot_R, pos_cm_foot_R] = this.calcFootCMTranforms(T_v2foot_r_proxim, T_v2foot_r_distal, T_v2heel_r, cm_pos_ratio_foot_heel2toe);
            [T_v2cm_foot_L, pos_cm_foot_L] = this.calcFootCMTranforms(T_v2foot_l_proxim, T_v2foot_l_distal, T_v2heel_l, cm_pos_ratio_foot_heel2toe);
            

            this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2cm_seg = T_v2cm_thigh_R;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).pos_cm_seg =  pos_cm_thigh_R;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).seg_mass =  sbj_mass*mass_ratio_thigh;
            
            this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2cm_seg = T_v2cm_thigh_L;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).pos_cm_seg = pos_cm_thigh_L;
            this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).seg_mass =  sbj_mass*mass_ratio_thigh;
            
            this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2cm_seg = T_v2cm_shank_R;
            this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).pos_cm_seg = pos_cm_shank_R;
            this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).seg_mass =  sbj_mass*mass_ratio_shank;
            
            this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2cm_seg = T_v2cm_shank_L;
            this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).pos_cm_seg = pos_cm_shank_L;
            this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).seg_mass =  sbj_mass*mass_ratio_shank;
            
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2cm_seg = T_v2cm_foot_R;
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).pos_cm_seg = pos_cm_foot_R;
            this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).seg_mass =  sbj_mass*mass_ratio_foot;
            
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2cm_seg = T_v2cm_foot_L;
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).pos_cm_seg = pos_cm_foot_L;
            this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).seg_mass =  sbj_mass*mass_ratio_foot;
            
            %% CoM of total body 
            seg_masses = [this.sbj_anthro(trial_no).body_segment_transform.seg_mass];
            weighted_sum = zeros(length(T_v2cm_foot_L), 3);
            for i = 1:length(seg_masses)
                weighted_sum = weighted_sum + seg_masses(i)*this.sbj_anthro(trial_no).body_segment_transform(i).pos_cm_seg;               
            end
            this.sbj_anthro(trial_no).total_cm_pos = weighted_sum./sum(seg_masses);
            
            disp(['Updated anthropometric segmental inertia in trial no. ', num2str(trial_no)])
        end
        
        %% Transformation Functions
        function T = T_translate(~, vec)
            T = eye(4);
            T(1:3,4) = vec;
        end
        
        function S = skew(~, v)
            % vec: 3 x 1 vector column
            S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        end
        
        function unit_V = unitVecMat(~, vec_mat, vec_dim)
            % UnitVecmat: Calculate the unit vector from a Matrix of
            % n vectors (if n x 3: use vec_dim = 2; 3 x n: use vec_dim = 1)
            unit_V = vec_mat./vecnorm(vec_mat, 2, vec_dim); % norm of each row
        end
        
        function T_3dmat = createTransforms(~, origin_mat, i_hat_mat, j_hat_mat, k_hat_mat)
            T_3dmat = repmat(eye(4), 1, 1, length(origin_mat));
            T_3dmat(1:3, 1, :) = reshape(i_hat_mat', 3, 1, []);
            T_3dmat(1:3, 2, :) = reshape(j_hat_mat', 3, 1, []);
            T_3dmat(1:3, 3, :) = reshape(k_hat_mat', 3, 1, []);
            T_3dmat(1:3, 4, :) = reshape(origin_mat', 3, 1, []);
        end
        
        function T_3dmat = createRotationAxisAngle(this, axis_mat, angle)
            % axis_mat: n x 3 double array
            % angle: 1 x n double array (rad)
            T_3dmat = repmat(eye(4), 1, 1, length(axis_mat));
            w = this.unitVecMat(axis_mat, 2); % normalize the axis vectors
            for i = 1:length(axis_mat)
                c = cos(angle(i)); s = sin(angle(i)); v = 1 - c;
                % I3*c + skew(w)*s + wT*w*v
                T_3dmat(1:3, 1:3, i) = eye(3)*c + this.skew(w(i, :))*s + w(i, :)'*w(i, :)*v;
            end
        end
        
        function T_3dmat = createTransformsTranslation(~, vec_mat)
            % vec_mat: n x 3 double array
            T_3dmat = repmat(eye(4), 1, 1, length(vec_mat));
            T_3dmat(1:3, 4, :) = reshape(vec_mat', 3, 1, []);
        end
        
        function invT_3dmat = invTransformsMat(~, T_3dmat)
            invT_3dmat = repmat(eye(4), 1, 1, length(T_3dmat));
            for i = 1:length(T_3dmat)
                R = T_3dmat(1:3, 1:3, i);
                t = T_3dmat(1:3, 4, i);
                invT_3dmat(:, :, i) = [R', -R'*t; [0 0 0 1]];
            end
        end
        
        function T_mul = multiplyTransforms(~, TL, TR)
            T_mul = repmat(eye(4), 1, 1, length(TL));
            for i = 1 : length(TL)
                T_mul(:, :, i) = TL(:, :, i)*TR(:, :, i);
            end
        end
        
        function pos_wrt_new_frame = calcPosInNewFrame(this, T_v2new_frame, pos_wrt_vicon)
            T_pos_wrt_vicon = this.createTransformsTranslation(pos_wrt_vicon);
            T_pos_wrt_new_frame = this.multiplyTransforms(this.invTransformsMat(T_v2new_frame), T_pos_wrt_vicon);
            pos_wrt_new_frame = reshape(T_pos_wrt_new_frame(1:3, 4, :), 3, [])';
        end
        
        function z_pos_wrt_new_frame = calcZPosInNewFrame(this, T_v2new_frame, pos_wrt_vicon)
            T_pos_wrt_vicon = this.createTransformsTranslation(pos_wrt_vicon);
            T_pos_wrt_new_frame = this.multiplyTransforms(this.invTransformsMat(T_v2new_frame), T_pos_wrt_vicon);
            z_pos_wrt_new_frame = reshape(T_pos_wrt_new_frame(1:3, 4, :), 3, [])';
            z_pos_wrt_new_frame(:, 1:2) = 0;
        end
        
        function [T_v2cm_seg, pos_cm_seg] = calcCMTransforms(this, T_v2sup_frame, pos_sup_landmark, pos_inf_landmark, ratio_sup2inf)
            z_sup_wrt_sup_frame = this.calcZPosInNewFrame(T_v2sup_frame, pos_sup_landmark);
            z_inf_wrt_sup_frame = this.calcZPosInNewFrame(T_v2sup_frame, pos_inf_landmark);
            T_sup_frame2z_sup = this.createTransformsTranslation(z_sup_wrt_sup_frame);
            T_v2z_sup = this.multiplyTransforms(T_v2sup_frame, T_sup_frame2z_sup);
            z_sup2cm_seg = ratio_sup2inf*(z_inf_wrt_sup_frame - z_sup_wrt_sup_frame);
            T_z_sup2cm_seg = this.createTransformsTranslation(z_sup2cm_seg);
            T_v2cm_seg = this.multiplyTransforms(T_v2z_sup , T_z_sup2cm_seg);
            % get only the position vectors
            pos_cm_seg = reshape(T_v2cm_seg(1:3, 4, :), 3, [])';
        end
        
        function [T_v2cm_seg, pos_cm_seg] = calcLimbCMTranforms(this, T_v2sup_frame, T_v2inf_frame, ratio_sup2inf)
            % get origins of both ends
            sup_joint_pos = reshape(T_v2sup_frame(1:3, 4, :), 3, [])';
            inf_joint_pos = reshape(T_v2inf_frame(1:3, 4, :), 3, [])';
            
            z_sup_wrt_sup_frame = this.calcZPosInNewFrame(T_v2sup_frame, sup_joint_pos);
            z_inf_wrt_sup_frame = this.calcZPosInNewFrame(T_v2sup_frame, inf_joint_pos);
            T_sup_frame2z_sup = this.createTransformsTranslation(z_sup_wrt_sup_frame);
            T_v2z_sup = this.multiplyTransforms(T_v2sup_frame, T_sup_frame2z_sup); % may be a redundant step
            z_sup2cm_seg = ratio_sup2inf*(z_inf_wrt_sup_frame - z_sup_wrt_sup_frame);
            T_z_sup2cm_seg = this.createTransformsTranslation(z_sup2cm_seg);
            T_v2cm_seg = this.multiplyTransforms(T_v2z_sup , T_z_sup2cm_seg);
            % get only the position vectors
            pos_cm_seg = reshape(T_v2cm_seg(1:3, 4, :), 3, [])';
        end
        
         function [T_v2cm_foot, pos_cm_foot] = calcFootCMTranforms(this, T_v2foot_proxim, T_v2foot_distal, T_v2heel, cm_ratio_heel2toe)
            % get origins of both ends
            heel_pos = reshape(T_v2heel(1:3, 4, :), 3, [])';
            toe_pos = reshape(T_v2foot_distal(1:3, 4, :), 3, [])';
            ankle_pos = reshape(T_v2foot_proxim(1:3, 4, :), 3, [])';    
            
            heel_pos_wrt_heel_frame = this.calcPosInNewFrame(T_v2heel, heel_pos); % should be zeros
            toe_pos_wrt_heel_frame = this.calcPosInNewFrame(T_v2heel, toe_pos); % should be along y axis only
            ankle_pos_wrt_heel_frame = this.calcPosInNewFrame(T_v2heel, ankle_pos ); % should be along y and z axes 
            ankle_pos_wrt_heel_frame(:, [1,2]) = 0; % use only the z component
            
            heel2cm_foot = cm_ratio_heel2toe*(toe_pos_wrt_heel_frame - heel_pos_wrt_heel_frame) + ...
                0.5*ankle_pos_wrt_heel_frame(); % go up half of the foot height
            T_heel2cm_foot = this.createTransformsTranslation(heel2cm_foot);
            
            T_v2cm_foot = this.multiplyTransforms(T_v2heel , T_heel2cm_foot);
            % get only the position vectors
            pos_cm_foot = reshape(T_v2cm_foot(1:3, 4, :), 3, [])';
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
            scat = scatter3(var(:,1), var(:,2), var(:, 3));
            scat.Marker = '.';
            grid on; grid minor;
            axis equal
            title([plate_name, ': CoP Trajectory'])
        end
        
        function plotSegmentTrans(this, trial_no, viz_time_step)
            coordinate_scale = 50;
            % pelvis
            pelvis_segment_indx = strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Pelvis');
            %             Tp1 = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Tp2 = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            %             this.plotCoordinateTransform(Tp1, 150)
            Thip_r = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal_aux_r(:,:, viz_time_step);
            Thip_l = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal_aux_l(:,:, viz_time_step);
            this.plotCoordinateTransform(Tp2, coordinate_scale)
            this.plotCoordinateTransform(Thip_r, coordinate_scale)
            this.plotCoordinateTransform(Thip_l, coordinate_scale)
            % plot line from pelvis center to hips
            pelv2hip_r_line = [Tp2(1:3, 4),Thip_r(1:3, 4)];
            pelv2hip_l_line = [Tp2(1:3, 4),Thip_l(1:3, 4)];
            plot3(pelv2hip_r_line(1, :), pelv2hip_r_line(2, :), pelv2hip_r_line(3, :), 'k:');
            plot3(pelv2hip_l_line(1, :), pelv2hip_l_line(2, :), pelv2hip_l_line(3, :), 'k:');
            
            % lumbar
            lumber_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Lumbar'));
            Tl1 = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Tl2 = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(Tl1, coordinate_scale)
            this.plotCoordinateTransform(Tl2, coordinate_scale)
            % plot segment line
            lumber_line = [Tl1(1:3, 4),Tl2(1:3, 4)];
            plot3(lumber_line(1, :), lumber_line(2, :), lumber_line(3, :), 'k:');
            
            % thorax
            thorax_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thorax'));
            Tt1 = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Tt2 = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            Tshldr_r = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal_aux_r(:,:, viz_time_step);
            Tshldr_l = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal_aux_l(:,:, viz_time_step);
            this.plotCoordinateTransform(Tt1, coordinate_scale)
            this.plotCoordinateTransform(Tt2, coordinate_scale)
            this.plotCoordinateTransform(Tshldr_r, coordinate_scale)
            this.plotCoordinateTransform(Tshldr_l, coordinate_scale)
            % plot segment line
            thorax_line = [Tt1(1:3, 4),Tt2(1:3, 4)];
            thorax2shldr_r_line = [Tt2(1:3, 4),Tshldr_r(1:3, 4)];
            thorax2shldr_l_line = [Tt2(1:3, 4),Tshldr_l(1:3, 4)];
            plot3(thorax_line(1, :), thorax_line(2, :), thorax_line(3, :), 'k:');
            plot3(thorax2shldr_r_line(1, :), thorax2shldr_r_line(2, :), thorax2shldr_r_line(3, :), 'k:');
            plot3(thorax2shldr_l_line(1, :), thorax2shldr_l_line(2, :), thorax2shldr_l_line(3, :), 'k:');
            
            % upperarm
            uarm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_R'));
            uarm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_L'));
            T_uarm_r_prox = this.sbj_anthro(trial_no).body_segment_transform(uarm_r_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_uarm_l_prox = this.sbj_anthro(trial_no).body_segment_transform(uarm_l_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_uarm_r_dist = this.sbj_anthro(trial_no).body_segment_transform(uarm_r_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            T_uarm_l_dist = this.sbj_anthro(trial_no).body_segment_transform(uarm_l_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(T_uarm_r_prox, coordinate_scale)
            this.plotCoordinateTransform(T_uarm_l_prox, coordinate_scale)
            this.plotCoordinateTransform(T_uarm_r_dist, coordinate_scale)
            this.plotCoordinateTransform(T_uarm_l_dist, coordinate_scale)
            % plot segment lines
            uram_r_line = [T_uarm_r_prox(1:3, 4), T_uarm_r_dist(1:3, 4)];
            uram_l_line = [T_uarm_l_prox(1:3, 4), T_uarm_l_dist(1:3, 4)];
            plot3(uram_r_line(1, :), uram_r_line(2, :), uram_r_line(3, :), 'k:');
            plot3(uram_l_line(1, :), uram_l_line(2, :), uram_l_line(3, :), 'k:');
            
            % forearm
            forearm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_R'));
            forearm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_L'));
            T_farm_r_prox = this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_farm_l_prox = this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_farm_r_dist = this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            T_farm_l_dist = this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(T_farm_r_prox, coordinate_scale)
            this.plotCoordinateTransform(T_farm_l_prox, coordinate_scale)
            this.plotCoordinateTransform(T_farm_r_dist, coordinate_scale)
            this.plotCoordinateTransform(T_farm_l_dist, coordinate_scale)
            % plot segment lines
            fram_r_line = [T_farm_r_prox(1:3, 4), T_farm_r_dist(1:3, 4)];
            fram_l_line = [T_farm_l_prox(1:3, 4), T_farm_l_dist(1:3, 4)];
            plot3(fram_r_line(1, :), fram_r_line(2, :), fram_r_line(3, :), 'k:');
            plot3(fram_l_line(1, :), fram_l_line(2, :), fram_l_line(3, :), 'k:');
            
            % hand
            hand_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_R'));
            hand_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_L'));
            T_hand_r_prox = this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_hand_l_prox = this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_hand_r_dist = this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            T_hand_l_dist = this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(T_hand_r_prox, coordinate_scale)
            this.plotCoordinateTransform(T_hand_l_prox, coordinate_scale)
            this.plotCoordinateTransform(T_hand_r_dist, coordinate_scale)
            this.plotCoordinateTransform(T_hand_l_dist, coordinate_scale)
            % plot segment lines
            hand_r_line = [T_hand_r_prox(1:3, 4), T_hand_r_dist(1:3, 4)];
            hand_l_line = [T_hand_l_prox(1:3, 4), T_hand_l_dist(1:3, 4)];
            plot3(hand_r_line(1, :), hand_r_line(2, :), hand_r_line(3, :), 'k:');
            plot3(hand_l_line(1, :), hand_l_line(2, :), hand_l_line(3, :), 'k:');
            
            % head
            head_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Head'));
            Th1 = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Th2 = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(Th1, coordinate_scale)
            this.plotCoordinateTransform(Th2, coordinate_scale)
            % plot segment line
            head_line = [Th1(1:3, 4),Th2(1:3, 4)];
            plot3(head_line(1, :), head_line(2, :), head_line(3, :), 'k:');
            
            % thigh
            thigh_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_R'));
            thigh_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_L'));
            T_thigh_r_prox = this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_thigh_r_dist = this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            T_thigh_l_prox = this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_thigh_l_dist = this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(T_thigh_r_prox, coordinate_scale)
            this.plotCoordinateTransform(T_thigh_r_dist, coordinate_scale)
            this.plotCoordinateTransform(T_thigh_l_prox, coordinate_scale)
            this.plotCoordinateTransform(T_thigh_l_dist, coordinate_scale)
            % plot segment lines
            thigh_r_line = [T_thigh_r_prox(1:3, 4), T_thigh_r_dist(1:3, 4)];
            thigh_l_line = [T_thigh_l_prox(1:3, 4), T_thigh_l_dist(1:3, 4)];
            plot3(thigh_r_line(1, :), thigh_r_line(2, :), thigh_r_line(3, :), 'k:');
            plot3(thigh_l_line(1, :), thigh_l_line(2, :), thigh_l_line(3, :), 'k:');
            
            % shank
            shank_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_R'));
            shank_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_L'));
            T_shank_r_prox = this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_shank_r_dist = this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            T_shank_l_prox = this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_shank_l_dist = this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(T_shank_r_prox, coordinate_scale)
            this.plotCoordinateTransform(T_shank_r_dist, coordinate_scale)
            this.plotCoordinateTransform(T_shank_l_prox, coordinate_scale)
            this.plotCoordinateTransform(T_shank_l_dist, coordinate_scale)
            % plot segment lines
            shank_r_line = [T_shank_r_prox(1:3, 4), T_shank_r_dist(1:3, 4)];
            shank_l_line = [T_shank_l_prox(1:3, 4), T_shank_l_dist(1:3, 4)];
            plot3(shank_r_line(1, :), shank_r_line(2, :), shank_r_line(3, :), 'k:');
            plot3(shank_l_line(1, :), shank_l_line(2, :), shank_l_line(3, :), 'k:');
            
            % foot
            foot_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_R'));
            foot_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_L'));
            T_foot_r_prox = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_foot_r_dist = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            T_heel_r = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2seg_distal_aux_r(:,:, viz_time_step);
            
            T_foot_l_prox = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            T_foot_l_dist = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            T_heel_l = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2seg_distal_aux_l(:,:, viz_time_step);
            
            this.plotCoordinateTransform(T_foot_r_prox, coordinate_scale)
            this.plotCoordinateTransform(T_foot_r_dist, coordinate_scale)
            this.plotCoordinateTransform(T_heel_r, coordinate_scale)
            
            this.plotCoordinateTransform(T_foot_l_prox, coordinate_scale)
            this.plotCoordinateTransform(T_foot_l_dist, coordinate_scale)
            this.plotCoordinateTransform(T_heel_l, coordinate_scale)
            
            % plot segment lines
            ankle_r_line = [T_foot_r_prox(1:3, 4), T_foot_r_dist(1:3, 4), T_heel_r(1:3, 4), T_foot_r_prox(1:3, 4)];
            ankle_l_line = [T_foot_l_prox(1:3, 4), T_foot_l_dist(1:3, 4), T_heel_l(1:3, 4), T_foot_l_prox(1:3, 4)];
            
            plot3(ankle_r_line(1, :), ankle_r_line(2, :), ankle_r_line(3, :), 'k:');
            plot3(ankle_l_line(1, :), ankle_l_line(2, :), ankle_l_line(3, :), 'k:');
            
        end
        
        function plotSegmentCM(this, trial_no, viz_time_step)
            coordinate_scale_CM = 25;
            marker_full_size = 2000;
            
            anthro_pelvis_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Pelvis');
            anthro_lumbar_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Lumbar');
            anthro_thorax_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Thorax');
            anthro_head_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Head');
            anthro_uarm_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'UpperArm');
            anthro_farm_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Forearm');
            anthro_hand_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Hand');
            anthro_thigh_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Thigh');
            anthro_shank_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Shank');
            anthro_foot_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Foot');
            
            pelvis_segment_indx = strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Pelvis');
            lumber_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Lumbar'));
            thorax_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thorax'));
            head_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Head'));
            upper_arm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_R'));
            upper_arm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'UpperArm_L'));
            forearm_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_R'));
            forearm_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Forearm_L'));
            hand_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_R'));
            hand_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Hand_L'));
            thigh_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_R'));
            thigh_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thigh_L'));
            shank_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_R'));
            shank_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Shank_L'));
            foot_r_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_R'));
            foot_l_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Foot_L'));
            
            mass_ratio_pelvis = this.sbj_anthro_table(anthro_pelvis_table_indx).mass_ratio; % male 0.1117
            mass_ratio_lumbar = this.sbj_anthro_table(anthro_lumbar_table_indx).mass_ratio; % 0.1633
            mass_ratio_thorax = this.sbj_anthro_table(anthro_thorax_table_indx).mass_ratio; % 0.1596
            mass_ratio_head = this.sbj_anthro_table(anthro_head_table_indx).mass_ratio; % 0.0694
            mass_ratio_uarm = this.sbj_anthro_table(anthro_uarm_table_indx).mass_ratio; % 0.0271
            mass_ratio_farm = this.sbj_anthro_table(anthro_farm_table_indx).mass_ratio; % 0.0162
            mass_ratio_hand = this.sbj_anthro_table(anthro_hand_table_indx).mass_ratio; % 0.0061
            mass_ratio_thigh = this.sbj_anthro_table(anthro_thigh_table_indx).mass_ratio; % 0.1416
            mass_ratio_shank = this.sbj_anthro_table(anthro_shank_table_indx).mass_ratio; % 0.0433
            mass_ratio_foot = this.sbj_anthro_table(anthro_foot_table_indx).mass_ratio; % 0.0137
            
            % pelvis CM
            T_v2cm_pelvis = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_pelvis = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_pelvis, coordinate_scale_CM);
            scatter3(pos_cm_pelvis(1), pos_cm_pelvis(2), pos_cm_pelvis(3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_pelvis*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % lumbar CM
            T_v2cm_lumbar = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_lumbar = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_lumbar, coordinate_scale_CM);
            scatter3(pos_cm_lumbar(1), pos_cm_lumbar(2), pos_cm_lumbar(3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_lumbar*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % thorax CM
            T_v2cm_thorax = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_thorax = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_thorax, coordinate_scale_CM);
            scatter3(pos_cm_thorax(1), pos_cm_thorax(2), pos_cm_thorax(3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_thorax*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % head CM
            T_v2cm_head = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_head = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_head, coordinate_scale_CM);
            scatter3(pos_cm_head(1), pos_cm_head(2), pos_cm_head(3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_head*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % upper arm CoM
            T_v2cm_uarm_R = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_uarm_R = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_r_segment_indx).pos_cm_seg(viz_time_step, :);
            T_v2cm_uarm_L = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_uarm_L = this.sbj_anthro(trial_no).body_segment_transform(upper_arm_l_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_uarm_L, coordinate_scale_CM);
            this.plotCoordinateTransform(T_v2cm_uarm_R, coordinate_scale_CM);
            pos_cm_uarm = [pos_cm_uarm_R; pos_cm_uarm_L];
            scatter3(pos_cm_uarm(:, 1), pos_cm_uarm(:, 2), pos_cm_uarm(:, 3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_uarm*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % forearm CoM
            T_v2cm_farm_R = this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_farm_R = this.sbj_anthro(trial_no).body_segment_transform(forearm_r_segment_indx).pos_cm_seg(viz_time_step, :);
            T_v2cm_farm_L = this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_farm_L = this.sbj_anthro(trial_no).body_segment_transform(forearm_l_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_farm_L, coordinate_scale_CM);
            this.plotCoordinateTransform(T_v2cm_farm_R, coordinate_scale_CM);
            pos_cm_farm = [pos_cm_farm_R; pos_cm_farm_L];
            scatter3(pos_cm_farm(:, 1), pos_cm_farm(:, 2), pos_cm_farm(:, 3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_farm*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % hand CoM
            T_v2cm_hand_R = this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_hand_R = this.sbj_anthro(trial_no).body_segment_transform(hand_r_segment_indx).pos_cm_seg(viz_time_step, :);
            T_v2cm_hand_L = this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_hand_L = this.sbj_anthro(trial_no).body_segment_transform(hand_l_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_hand_L, coordinate_scale_CM);
            this.plotCoordinateTransform(T_v2cm_hand_R, coordinate_scale_CM);
            pos_cm_hand = [pos_cm_hand_R; pos_cm_hand_L];
            scatter3(pos_cm_hand(:, 1), pos_cm_hand(:, 2), pos_cm_hand(:, 3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_hand*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % thigh CoM
            T_v2cm_thigh_R = this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            T_v2cm_thigh_L = this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_thigh_R = this.sbj_anthro(trial_no).body_segment_transform(thigh_r_segment_indx).pos_cm_seg(viz_time_step, :);
            pos_cm_thigh_L = this.sbj_anthro(trial_no).body_segment_transform(thigh_l_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_thigh_R, coordinate_scale_CM);
            this.plotCoordinateTransform(T_v2cm_thigh_L, coordinate_scale_CM);
            pos_cm_thigh = [pos_cm_thigh_R; pos_cm_thigh_L];
            scatter3(pos_cm_thigh(:, 1), pos_cm_thigh(:, 2), pos_cm_thigh(:, 3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_thigh*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % shank CoM
            T_v2cm_shank_R = this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            T_v2cm_shank_L = this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_shank_R = this.sbj_anthro(trial_no).body_segment_transform(shank_r_segment_indx).pos_cm_seg(viz_time_step, :);
            pos_cm_shank_L = this.sbj_anthro(trial_no).body_segment_transform(shank_l_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_shank_L, coordinate_scale_CM);
            this.plotCoordinateTransform(T_v2cm_shank_R, coordinate_scale_CM);
            pos_cm_shank = [pos_cm_shank_R; pos_cm_shank_L];
            scatter3(pos_cm_shank(:, 1), pos_cm_shank(:, 2), pos_cm_shank(:, 3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_shank*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % Foot CoM
            T_v2cm_foot_R = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            T_v2cm_foot_L = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_foot_R = this.sbj_anthro(trial_no).body_segment_transform(foot_r_segment_indx).pos_cm_seg(viz_time_step, :);
            pos_cm_foot_L = this.sbj_anthro(trial_no).body_segment_transform(foot_l_segment_indx).pos_cm_seg(viz_time_step, :);
            this.plotCoordinateTransform(T_v2cm_foot_L, coordinate_scale_CM);
            this.plotCoordinateTransform(T_v2cm_foot_R, coordinate_scale_CM);
            pos_cm_foot = [pos_cm_foot_R; pos_cm_foot_L];
            scatter3(pos_cm_foot(:, 1), pos_cm_foot(:, 2), pos_cm_foot(:, 3), 'MarkerFaceColor', 'r',...
                'SizeData', mass_ratio_foot*marker_full_size, 'MarkerEdgeColor', 'k');
            
            % total CoM
            pos_cm_total = this.sbj_anthro(trial_no).total_cm_pos(viz_time_step, :);
            scatter3(pos_cm_total(:, 1), pos_cm_total(:, 2), pos_cm_total(:, 3), 'MarkerFaceColor', 'b',...
                'SizeData', 0.25*marker_full_size, 'MarkerEdgeColor', 'k');
            
        end
        
        function vizTrial(this, trial_no, viz_time_step)
            % vizTrial: visualization of each trial
            figure('pos',[100 100 900 600])
            this.plotCoordinateTransform(eye(4), 250); hold on; % vicon origin
            this.plotTrajCoP(trial_no, 'Seat Plate');
            this.plotTrajCoP(trial_no, 'Foot Plate');
            
            % show initial pos of marker clusters of both rings
            pelvis_cluster_no = strcmp({this.sbj_WRAPS2(trial_no).trial_transform_data.cluster_name}, 'Pelvis Brace');
            thorax_cluster_no = strcmp({this.sbj_WRAPS2(trial_no).trial_transform_data.cluster_name}, 'Thorax Brace');
            
            init_pelvis_marker_pos = reshape(this.sbj_WRAPS2(trial_no).trial_transform_data(pelvis_cluster_no).marker_pos(viz_time_step, :, :), 3, []);
            init_thorax_marker_pos = reshape(this.sbj_WRAPS2(trial_no).trial_transform_data(thorax_cluster_no).marker_pos(viz_time_step, :, :), 3, []);
            scatter3(init_pelvis_marker_pos(1, :), init_pelvis_marker_pos(2, :), init_pelvis_marker_pos(3, :))
            scatter3(init_thorax_marker_pos(1, :), init_thorax_marker_pos(2, :), init_thorax_marker_pos(3, :))
            
            title(this.raw_data(trial_no).marker_trial_name)
            %             for time_step = 1: 100 : length(this.sbj_WRAPS2(trial_no).trial_transform_data(pelvis_cluster_no).transforms_vicon2seg)
            T_v2pelvis = this.sbj_WRAPS2(trial_no).trial_transform_data(pelvis_cluster_no).transforms_vicon2seg(:,:,viz_time_step);
            T_v2thorax = this.sbj_WRAPS2(trial_no).trial_transform_data(thorax_cluster_no).transforms_vicon2seg(:,:,viz_time_step);
            this.plotCoordinateTransform(T_v2pelvis, 75);
            this.plotCoordinateTransform(T_v2thorax, 75);
            %             end
            
            % landmark pos
            landmark_vec = reshape(this.sbj_anthro(trial_no).landmark_pos(viz_time_step,:,:), 3, []);
            scatter3(landmark_vec(1,:), landmark_vec(2,:), landmark_vec(3,:), '*');
            
            % segment transforms
            this.plotSegmentTrans(trial_no, viz_time_step)
            this.plotSegmentCM(trial_no, viz_time_step)
            
        end
        
        % visual elements
        function plotCoordinateTransform(~, T, scale)
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

