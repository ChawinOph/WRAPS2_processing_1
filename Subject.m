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
        function this = Subject(sbj_index, sbj_folder_names, sbj_ids, sbj_names, sbj_gender, marker_cluster_pos, anthropomet_table)
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
        
        %% Calculate landmark positions
        
        function calcLandmarkPos(this, trial_no)
            % pelvis landmark
            this.sbj_anthro(trial_no).trial_name = this.raw_data(trial_no).marker_trial_name;
            this.sbj_anthro(trial_no).torso_landmark_names = {'UMBLC', 'R_ASIS', 'L_ASIS', 'R_PSIS', 'L_PSIS', 'R_HIP', 'L_HIP', 'SN', 'C7', 'XP', 'T8', 'VT', 'OC'};
            this.sbj_anthro(trial_no).landmark_pos = zeros(size(this.raw_data(trial_no).marker_data(1).marker_pos,1), 3, 12);
            this.sbj_anthro(trial_no).body_segment_transform = struct();
            this.sbj_anthro(trial_no).body_segment_transform(1).segment_name = 'Pelvis';
            this.sbj_anthro(trial_no).body_segment_transform(2).segment_name = 'Lumbar';
            this.sbj_anthro(trial_no).body_segment_transform(3).segment_name = 'Thorax';
            this.sbj_anthro(trial_no).body_segment_transform(4).segment_name = 'Head';
            
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
            pelvis_raw_indx = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Pelvis');
            P_Hip_R_raw_indx = strcmp(this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_names, 'P_Hip_R');
            P_Hip_L_raw_indx = strcmp(this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_names, 'P_Hip_L');
            
            thorax_raw_indx = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Thorax');
            sternal_notch_raw_indx = strcmp(this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_names, 'T_SternalNotch');
            C7_raw_indx = strcmp(this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_names, 'T_C7');
            xyphoid_raw_indx = strcmp(this.raw_data(trial_no).marker_data(thorax_raw_indx).marker_names, 'T_Xyphoid');
            
            head_raw_indx = strcmp({this.raw_data(trial_no).marker_data.segment_names}, 'Head');
            vertex_raw_indx = strcmp(this.raw_data(trial_no).marker_data(head_raw_indx).marker_names, 'H_Vertex');
            forehead_raw_indx = strcmp(this.raw_data(trial_no).marker_data(head_raw_indx).marker_names, 'H_Forehead');
            
            this.sbj_anthro(trial_no).landmark_pos(:,:,6) = this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_pos(:, :, P_Hip_R_raw_indx);
            this.sbj_anthro(trial_no).landmark_pos(:,:,7) = this.raw_data(trial_no).marker_data(pelvis_raw_indx).marker_pos(:, :, P_Hip_L_raw_indx);
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
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_proxim = repmat(eye(4), 1, 1, length(origin_pelvis_distal));
            this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal = this.createTransforms(origin_pelvis_distal, i_hat_pelvis, j_hat_pelvis, k_hat_pelvis);
            
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
            this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_proxim = this.createTransforms(origin_lumbar_proxim, i_hat_lumbar, j_hat_lumbar, k_hat_lumbar);
            this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_distal = this.createTransforms(origin_lumbar_distal, i_hat_lumbar, j_hat_lumbar, k_hat_lumbar);
                        
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
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_proxim = this.createTransforms(origin_thorax_proxim, i_hat_thorax, j_hat_thorax, k_hat_thorax);
            this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal = this.createTransforms(origin_thorax_distal, i_hat_thorax, j_hat_thorax, k_hat_thorax);
            
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
            this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_proxim = this.createTransforms(origin_head_proxim, i_hat_head, j_hat_head, k_hat_head);
            this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_distal = this.createTransforms(origin_head_distal, i_hat_head, j_hat_head, k_hat_head);
                      
            disp(['Updated anthropometric segment transformations in trial no. ', num2str(trial_no)])
            
        end
        
        function calcSegmentInertia(this, trial_no)
            % find landmarks projected on the local segment frame
            anthro_pelvis_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Pelvis');
            anthro_lumbar_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Lumbar');
            anthro_thorax_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Thorax');
            anthro_head_table_indx = strcmp({this.sbj_anthro_table.segment_name}, 'Head');
                             
            cm_pos_ratio_pelvis_umblc2hip = this.sbj_anthro_table(anthro_pelvis_table_indx).cm_pos_ratio  ; % male  0.6115
            cm_pos_ratio_lumbar_xp2umblc = this.sbj_anthro_table(anthro_lumbar_table_indx).cm_pos_ratio  ; % 0.4502
            cm_pos_ratio_thorax_sn2xp = this.sbj_anthro_table(anthro_thorax_table_indx).cm_pos_ratio  ; % 0.2999
            cm_pos_ratio_thorax_vt2c7 = this.sbj_anthro_table(anthro_head_table_indx).cm_pos_ratio  ; % 0.5002
            
            mass_ratio_pelvis = this.sbj_anthro_table(anthro_pelvis_table_indx).mass_ratio; % male 0.1117
            mass_ratio_lumbar = this.sbj_anthro_table(anthro_lumbar_table_indx).mass_ratio; % 0.1633
            mass_ratio_thorax = this.sbj_anthro_table(anthro_thorax_table_indx).mass_ratio; % 0.1596
            mass_ratio_head = this.sbj_anthro_table(anthro_head_table_indx).mass_ratio; % 0.0694
            
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
            
            T_v2pelv_dist = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal;
            T_v2lumb_dist = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_distal;
            T_v2thor_dist = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal;
            T_v2head_dist = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_distal;
            
            [T_v2cm_pelvis, pos_cm_pelvis] = this.calCMTransforms(T_v2pelv_dist, pos_umblc, pos_m_hip, cm_pos_ratio_pelvis_umblc2hip);
            [T_v2cm_lumbar, pos_cm_lumbar] = this.calCMTransforms(T_v2lumb_dist, pos_xp, pos_umblc, cm_pos_ratio_lumbar_xp2umblc);
            [T_v2cm_thorax, pos_cm_thorax] = this.calCMTransforms(T_v2thor_dist, pos_sn, pos_xp, cm_pos_ratio_thorax_sn2xp);
            [T_v2cm_head, pos_cm_head] = this.calCMTransforms(T_v2head_dist, pos_vt, pos_c7, cm_pos_ratio_thorax_vt2c7);
              
            
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
            
            disp(['Updated anthropometric segmental inertia in trial no. ', num2str(trial_no)])
        end
        
        
        %% Transformation Functions
        function T = T_translate(~, vec)
            T = eye(4);
            T(1:3,4) = vec;
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
        
        function T_3dmat = createTransformsTranslation(~, vec_mat)
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
        
         function [T_v2cm_seg, pos_cm_seg] = calCMTransforms(this, T_v2distal_frame, pos_sup_landmark, pos_inf_landmark, ratio_sup2inf)
             z_sup_wrt_distal_frame = this.calcZPosInNewFrame(T_v2distal_frame, pos_sup_landmark);
             z_inf_wrt_distal_frame = this.calcZPosInNewFrame(T_v2distal_frame, pos_inf_landmark);
             T_dist_frame2z_sup = this.createTransformsTranslation(z_sup_wrt_distal_frame);
             T_v2z_sup = this.multiplyTransforms(T_v2distal_frame, T_dist_frame2z_sup);            
             z_sup2cm_seg = ratio_sup2inf*(z_inf_wrt_distal_frame - z_sup_wrt_distal_frame);             
             T_z_sup2cm_seg = this.createTransformsTranslation(z_sup2cm_seg);
             T_v2cm_seg = this.multiplyTransforms(T_v2z_sup , T_z_sup2cm_seg);
             % get only the position vectors
             pos_cm_seg = reshape(T_v2cm_seg(1:3, 4, :), 3, [])';
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
            % pelvis
            pelvis_segment_indx = strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Pelvis');
%             Tp1 = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Tp2 = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2seg_distal(:,:, viz_time_step);
%             this.plotCoordinateTransform(Tp1, 150)
            this.plotCoordinateTransform(Tp2, 150)
            
            % lumbar
            lumber_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Lumbar'));
            Tl1 = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Tl2 = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(Tl1, 150)
            this.plotCoordinateTransform(Tl2, 150)
            % plot segment line
            lumber_line = [Tl1(1:3,4),Tl2(1:3,4)];
            plot3(lumber_line(1, :), lumber_line(2, :), lumber_line(3, :), 'k:');
            
            % thorax
            thorax_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thorax'));
            Tt1 = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Tt2 = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(Tt1, 150)
            this.plotCoordinateTransform(Tt2, 150)
             % plot segment line
            thorax_line = [Tt1(1:3,4),Tt2(1:3,4)];
            plot3(thorax_line(1, :), thorax_line(2, :), thorax_line(3, :), 'k:');
            
            % head
            head_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Head'));
            Th1 = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_proxim(:,:, viz_time_step);
            Th2 = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2seg_distal(:,:, viz_time_step);
            this.plotCoordinateTransform(Th1, 150)
            this.plotCoordinateTransform(Th2, 150)
            % plot segment line
            head_line = [Th1(1:3,4),Th2(1:3,4)];
            plot3(head_line(1, :), head_line(2, :), head_line(3, :), 'k:');
            
        end
        
        function plotSegmentCM(this, trial_no, viz_time_step)
            pelvis_segment_indx = strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Pelvis');
            lumber_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Lumbar'));
            thorax_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Thorax'));
            head_segment_indx = find(strcmp({this.sbj_anthro(trial_no).body_segment_transform.segment_name}, 'Head'));
            
            % pelvis CM
            T_v2cm_pelvis = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_pelvis = this.sbj_anthro(trial_no).body_segment_transform(pelvis_segment_indx).pos_cm_seg(viz_time_step, :);  
            this.plotCoordinateTransform(T_v2cm_pelvis, 50);
            scatter3(pos_cm_pelvis(1), pos_cm_pelvis(2), pos_cm_pelvis(3), 'MarkerFaceColor', 'r',...
                'SizeData', 50, 'MarkerEdgeColor', 'k');   
            
             % lumbar CM
            T_v2cm_lumbar = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_lumbar = this.sbj_anthro(trial_no).body_segment_transform(lumber_segment_indx).pos_cm_seg(viz_time_step, :);  
            this.plotCoordinateTransform(T_v2cm_lumbar, 50);
            scatter3(pos_cm_lumbar(1), pos_cm_lumbar(2), pos_cm_lumbar(3), 'MarkerFaceColor', 'r',...
                'SizeData', 50, 'MarkerEdgeColor', 'k');   
            
             % thorax CM
            T_v2cm_thorax = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_thorax = this.sbj_anthro(trial_no).body_segment_transform(thorax_segment_indx).pos_cm_seg(viz_time_step, :);  
            this.plotCoordinateTransform(T_v2cm_thorax, 50);
            scatter3(pos_cm_thorax(1), pos_cm_thorax(2), pos_cm_thorax(3), 'MarkerFaceColor', 'r',...
                'SizeData', 50, 'MarkerEdgeColor', 'k');   
            
             % head CM
            T_v2cm_head = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).T_v2cm_seg(:,:,viz_time_step);
            pos_cm_head = this.sbj_anthro(trial_no).body_segment_transform(head_segment_indx).pos_cm_seg(viz_time_step, :);  
            this.plotCoordinateTransform(T_v2cm_head, 50);
            scatter3(pos_cm_head(1), pos_cm_head(2), pos_cm_head(3), 'MarkerFaceColor', 'r',...
                'SizeData', 25, 'MarkerEdgeColor', 'k');            
        end
        
        function vizTrial(this, trial_no, viz_time_step)
            % vizTrial: visualization of each trial
            figure;
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
            %                 T_v2pelvis = this.sbj_WRAPS2(trial_no).trial_transform_data(pelvis_cluster_no).transforms_vicon2seg(:,:,time_step);
            %                 T_v2thorax = this.sbj_WRAPS2(trial_no).trial_transform_data(thorax_cluster_no).transforms_vicon2seg(:,:,time_step);
            %                 this.plotCoordinateTransform(T_v2pelvis, 100);
            %                 this.plotCoordinateTransform(T_v2thorax, 100);
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

