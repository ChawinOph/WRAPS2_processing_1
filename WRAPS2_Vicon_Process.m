% Process WRAPS 2.0 Data collection
clc; clear; close all;
tic
%% Setup folder directory
sbj_folder_names = {'w001_2018-09-18'};
sbj_idcs = {'w001'};
sbj_names = {'WRAPS2_Chawin'};
sbj_gender = {'male'};

load 'marker_cluster_pos.mat'; % can be edited in future like saving the cluster of (virtual) markers for static trials
load 'anthropomet_table'; % can be edited in future like saving the cluster of (virtual) markers for static trials
load 'sbj_measurement_table'; % can be edited in future like saving the cluster of (virtual) markers for static trials
%% Setup sorted marker names for each segments (May vary for different subjects)
pelvis_sorted_marker_names = {'P_Top_R', 'P_Front_L', 'P_Top_L', 'P_Back_L', 'P_Hip_R', 'P_ASIS_R', 'P_ASIS_L', 'P_Hip_L'};
thorax_sorted_marker_names = {'T_Front_R', 'T_Front_L', 'T_Side_L', 'T_SternalNotch', 'T_C7', 'T_T8_Top', 'T_T8_L', 'T_T8_R', 'T_Xyphoid'};
head_sorted_marker_names = {'H_Ear_R', 'H_Forehead', 'H_Ear_L', 'H_Occiput', 'H_Vertex'};
sorted_segment_names = {'Pelvis', 'Thorax', 'Head'};
sorted_marker_names = {pelvis_sorted_marker_names; thorax_sorted_marker_names; head_sorted_marker_names}; % segment no 1,2,3

%% Setup Force Plate names
forceplate_names = {'Seat Plate', 'Foot Plate'};
fplate_var_names = {'Force', 'CoP'};

%% Setup VICON markers trial file names by the order of subjects
%# Markers_static_baseline = {''};
markers_static_max_flexion_trials = {'Static Max Flexion03 Markers'}; % match with trial no. 1 (construct a 2D cell later)

markers_reach_max_knee_ipsi_trials = {'Max Reach Knee Ipsi Markers'}; % match with trial no. 2

markers_static_max_rbend_trials = {'Static Max R-Side Bend01 Markers'}; % match with trial no. 3

markers_reach_max_shoulder_center_trials = {'Max Reach Shoulder Center01 Markers'}; % trial no.4

markers_reach_max_shoulder_contra_trials = {'Max Reach Shoulder Contra01 Markers'}; % trial no. 5

markers_fast_flexion_extension_trials = {'Slow Flexion-Extension Markers'}; % trial no.6

markers_static_neutral_trials = {'Static Neutral02 Markers'}; % trial no.7

markers_static_max_extension_trials = {'Static Max Extension02 Markers'};

markers_static_max_lbend_trials = {'Static Max L-Side Bend01 Markers'};

markers_static_max_rtwist_trials = {'Static Max R-Twist Markers'};

markers_static_max_ltwsist_trials = {'Static Max L-Twist Markers'};

markers_reach_max_shoulder_ipsi_trials = {'Max Reach Shoulder Ipsi Markers'};

markers_reach_max_knee_contra_trials = {'Max Reach Knee Contra Markers'};

markers_reach_max_knee_center_trials = {'Max Reach Knee Center Markers'};

markers_slow_flexion_extension_trials = {'Slower Flexion-Extension Markers'};

markers_fast_lateral_trials = {'Slow Lateral Markers'};

markers_slow_lateral_trials = {'Slower Lateral Markers'};

markers_fast_twist_trials = {'Slow Twisting Markers'};

markers_slow_twist_trials = {'Slower Twisting Markers'};

markers_fast_rolling_trials = {'Slow Rolling Markers'};

markers_slow_rolling_trials = {'Slower Rolling01 Markers'};

%% Setup VICON force plate trial file names by the order of subjects
% forceplate_static_baseline = {'Static Max Flexion03 ForcePlate'};
force_plate_static_max_flexion_trials = {'Static Max Flexion03 ForcePlate'}; % trial no. 1

force_plate_reach_max_knee_ipsi_trials = {'Max Reach Knee Ipsi ForcePlate'}; % trial no. 2

force_plate_static_max_rbend_trials = {'Static Max R-Side Bend01 ForcePlate'}; % trial no. 3

force_plate_reach_max_shoulder_center_trials = {'Max Reach Shoulder Center01 ForcePlate'}; % trial no. 4

force_plate_reach_max_shoulder_contra_trials = {'Max Reach Shoulder Contra01 ForcePlate'}; % trial no. 5

force_plate_fast_flexion_extension_trials = {'Slow Flexion-Extension ForcePlate'}; % trial no.6

force_plate_static_neutral_trials = {'Static Neutral02 ForcePlate'}; % trial no.7

force_plate_static_max_extension_trials =  {'Static Max Extension02 ForcePlate'};

force_plate_static_max_lbend_trials = {'Static Max L-Side Bend01 ForcePlate'};

force_plate_static_max_rtwist_trials = {'Static Max R-Twist ForcePlate'};

force_plate_static_max_ltwsist_trials = {'Static Max L-Twist ForcePlate'};

force_plate_reach_max_shoulder_ipsi_trials = {'Max Reach Shoulder Ipsi ForcePlate'};

force_plate_reach_max_knee_contra_trials = {'Max Reach Knee Contra ForcePlate'};

force_plate_reach_max_knee_center_trials = {'Max Reach Knee Center ForcePlate'};

force_plate_slow_flexion_extension_trials = {'Slower Flexion-Extension ForcePlate'};

force_plate_fast_lateral_trials = {'Slow Lateral ForcePlate'};

force_plate_slow_lateral_trials = {'Slower Lateral ForcePlate'};

force_plate_fast_twist_trials = {'Slow Twisting ForcePlate'};

force_plate_slow_twist_trials = {'Slower Twisting ForcePlate'};

force_plate_fast_rolling_trials = {'Slow Rolling ForcePlate'};

force_plate_slow_rolling_trials = {'Slower Rolling01 ForcePlate'};

%% Import Vicon marker data
sbj_index = 1; trial_no = 1;
sbj1 = Subject(sbj_index, sbj_folder_names, sbj_idcs, sbj_names, sbj_gender, marker_cluster_pos, anthropomet_table, sbj_measurement_table);

%% import data
sbj1.importMarkerData_csv(trial_no, markers_static_max_flexion_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(trial_no, force_plate_static_max_flexion_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(2, markers_reach_max_knee_ipsi_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(2, force_plate_reach_max_knee_ipsi_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(3, markers_static_max_rbend_trials, sorted_segment_names, sorted_marker_names, {'P_Back_L'});
sbj1.importForcePlateData_csv(3, force_plate_static_max_rbend_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(4, markers_reach_max_shoulder_center_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(4, force_plate_reach_max_shoulder_center_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(5, markers_reach_max_shoulder_contra_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(5, force_plate_reach_max_shoulder_contra_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(6, markers_fast_flexion_extension_trials, sorted_segment_names, sorted_marker_names,  {'P_Hip_L'});
sbj1.importForcePlateData_csv(6, force_plate_fast_flexion_extension_trials, forceplate_names, fplate_var_names);

%% new files
sbj1.importMarkerData_csv(7, markers_static_neutral_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(7, force_plate_static_neutral_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(8, markers_static_max_extension_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(8, force_plate_static_max_extension_trials, forceplate_names, fplate_var_names);
 
sbj1.importMarkerData_csv(9, markers_static_max_lbend_trials, sorted_segment_names, sorted_marker_names,  {'P_Hip_L'});
sbj1.importForcePlateData_csv(9, force_plate_static_max_lbend_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(10, markers_static_max_rtwist_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(10, force_plate_static_max_rtwist_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(11, markers_static_max_ltwsist_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(11, force_plate_static_max_ltwsist_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(12, markers_reach_max_shoulder_ipsi_trials, sorted_segment_names, sorted_marker_names,  {'P_Top_R'});
sbj1.importForcePlateData_csv(12, force_plate_reach_max_shoulder_ipsi_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(13, markers_reach_max_knee_contra_trials, sorted_segment_names, sorted_marker_names,  {'P_Front_L'});
sbj1.importForcePlateData_csv(13, force_plate_reach_max_knee_contra_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(14, markers_reach_max_knee_center_trials, sorted_segment_names, sorted_marker_names,  {'P_Hip_L'});
sbj1.importForcePlateData_csv(14, force_plate_reach_max_knee_center_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(15, markers_slow_flexion_extension_trials, sorted_segment_names, sorted_marker_names, {'P_Hip_L'});
sbj1.importForcePlateData_csv(15, force_plate_slow_flexion_extension_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(16, markers_fast_lateral_trials, sorted_segment_names, sorted_marker_names, {'P_Hip_L'});
sbj1.importForcePlateData_csv(16, force_plate_fast_lateral_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(17, markers_slow_lateral_trials, sorted_segment_names, sorted_marker_names, {'P_Hip_L'});
sbj1.importForcePlateData_csv(17, force_plate_slow_lateral_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(18, markers_fast_twist_trials, sorted_segment_names, sorted_marker_names, {'P_Hip_L'});
sbj1.importForcePlateData_csv(18, force_plate_fast_twist_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(19, markers_slow_twist_trials, sorted_segment_names, sorted_marker_names, {'P_Hip_L'});
sbj1.importForcePlateData_csv(19, force_plate_slow_twist_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(20, markers_fast_rolling_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(20, force_plate_fast_rolling_trials, forceplate_names, fplate_var_names);

sbj1.importMarkerData_csv(21, markers_slow_rolling_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(21, force_plate_slow_rolling_trials, forceplate_names, fplate_var_names);

toc
%% Calculate CM of WRAPS using CAD and scale

%% visualize the data

for i = 1:21
    disp(['viztrial: ' num2str(i)])
    if i == 7
        sbj1.vizTrial(i, 1);
%         pause
        close all
    else
        sbj1.vizTrial(i, 1500);
%         pause
        close all
    end
end

sbj1.vizTrial(trial_no, 1000);
sbj1.plotISAdata(trial_no);
% 
% sbj1.vizTrial(2, 1000);
% sbj1.plotISAdata(2);
% 
sbj1.vizTrial(3, 1000);
sbj1.plotISAdata(3);
% 
% sbj1.vizTrial(4, 1000);
% sbj1.plotISAdata(4);
% 
% sbj1.vizTrial(5, 1000);
% sbj1.plotISAdata(5);

%% calc the force and moment
trial_no = 2;
total_cm_pos = sbj1.sbj_anthro(trial_no).total_cm_pos; % mm
mg = sbj1.g*sbj1.sbj_anthro_measurement.weight_kg; % N
mg_mat = repmat(mg', length(total_cm_pos), 1);
moment_mg = cross(total_cm_pos, mg_mat, 2);

% get the force plate data and down-sample it 
down_sample_ratio = sbj1.freq_fplate/sbj1.freq_marker;

seat_plate_indx = find(strcmp([sbj1.raw_data(trial_no).fplate_data.fplate_name],'Seat Plate'));
seat_cop_indx = find(strcmp([sbj1.raw_data(trial_no).fplate_data(seat_plate_indx).fplate_var_names],'CoP'));
seat_f_indx = find(strcmp([sbj1.raw_data(trial_no).fplate_data(seat_plate_indx).fplate_var_names],'Force'));

foot_plate_indx = find(strcmp([sbj1.raw_data(trial_no).fplate_data.fplate_name],'Foot Plate'));
foot_cop_indx = find(strcmp([sbj1.raw_data(trial_no).fplate_data(foot_plate_indx).fplate_var_names],'CoP'));
foot_f_indx = find(strcmp([sbj1.raw_data(trial_no).fplate_data(foot_plate_indx).fplate_var_names],'Force'));

CoP_seat = sbj1.raw_data(trial_no).fplate_data(seat_plate_indx).fplate_var(:,:,seat_cop_indx);
GRF_seat = sbj1.raw_data(trial_no).fplate_data(seat_plate_indx).fplate_var(:,:,seat_f_indx);

CoP_foot = sbj1.raw_data(trial_no).fplate_data(foot_plate_indx).fplate_var(:,:,foot_cop_indx);
GRF_foot = sbj1.raw_data(trial_no).fplate_data(foot_plate_indx).fplate_var(:,:,foot_f_indx);

moment_seat = -cross(CoP_seat, GRF_seat, 2);
moment_foot = -cross(CoP_foot, GRF_foot, 2);

moment_seat = moment_seat(1:down_sample_ratio:end, :);
moment_foot = moment_foot(1:down_sample_ratio:end, :);
GRF_seat = -GRF_seat(1:down_sample_ratio:end, :);
GRF_foot = -GRF_foot(1:down_sample_ratio:end, :);

T = sbj1.calcViconTime(length(moment_mg));
figure; 
plot(T, moment_mg, '-', T, moment_seat, ':', T, moment_foot);
title('moment separate')

figure; 
plot(T, (moment_mg + moment_seat + moment_foot)/1000);
title('moment')
figure;
plot(T, mg_mat + GRF_seat + GRF_foot);
title('force')

%              u = uicontrol('Style','slider','Position',[10 50 20 340],...
%                 'Min',1,'Max',16,'Value',1);


