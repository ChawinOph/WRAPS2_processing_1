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
sorted_segment_names = {'Pelvis','Thorax','Head'};
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

%% Setup VICON force plate trial file names by the order of subjects
% forceplate_static_baseline = {'Static Max Flexion03 ForcePlate'};
force_plate_static_max_flexion_trials = {'Static Max Flexion03 ForcePlate'}; % trial no. 1

force_plate_reach_max_knee_ipsi_trials = {'Max Reach Knee Ipsi ForcePlate'}; % trial no. 2

force_plate_static_max_rbend_trials = {'Static Max R-Side Bend01 ForcePlate'}; % trial no. 3

force_plate_reach_max_shoulder_center_trials = {'Max Reach Shoulder Center01 ForcePlate'}; % trial no. 4

force_plate_reach_max_shoulder_contra_trials = {'Max Reach Shoulder Contra01 ForcePlate'}; % trial no. 5

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

% Calculate transformation of the WRAPS2.0 from vicon
toc
sbj1.vizTrial(1, 1000)
sbj1.vizTrial(3, 1000)

%% Calculate CM of WRAPS using CAD and scale

%% Calculate the Velocities
T = sbj1.sbj_WRAPS2(trial_no).T_base2end;
time = 0: 1/sbj1.freq_marker: (length(T) - 1)/sbj1.freq_marker;
var = reshape(T(1:3, 4, :), 3, [])';
[time, dvar_fwd] = sbj1.calcFirstOrderDerivative(time, var, 'forward');
[time, dvar_bwd] = sbj1.calcFirstOrderDerivative(time, var, 'backward');
[time, dvar_cen] = sbj1.calcFirstOrderDerivative(time, var, 'center');
[time, dvar_cen5] = sbj1.calcFirstOrderDerivative(time, var, 'center_5point');

figure;
subplot(2,1,1);
plot(time, var); hold on
legend('x', 'y', 'z')
subplot(2,1, 2);
plot(time, dvar_cen, '--'); hold on
legend('v_x', 'v_y', 'v_z')
% plot(time, dvar_bwd, ':r'); 
% plot(time, dvar_cen, '--b'); 

%% smooth the marker pos on the brace before finding the velocities of each marker and IARs
% [1] A. Page, P. Candelas, and F. Belmar, “On the use of local fitting techniques 
% for the analysis of physical dynamic systems,” Eur. J. Phys., vol. 27, no. 2, p. 273, 2006.
% [2] A. Page, H. de Rosario, V. Mata, R. Porcar, J. Solaz, and M. J. Such, “Kinematics of the 
% trunk in sitting posture: An analysis based on the instantaneous axis of rotation,” Ergonomics, 
% vol. 52, no. 6, pp. 695–706, Jun. 2009.
pelvis_brace_marker_pos = sbj1.sbj_WRAPS2(trial_no).trial_transform_data(1).marker_pos;
thorax_brace_marker_pos = sbj1.sbj_WRAPS2(trial_no).trial_transform_data(2).marker_pos;
T = 0: 1/sbj1.freq_marker : (length(pelvis_brace_marker_pos) - 1)/sbj1.freq_marker;
var = pelvis_brace_marker_pos(:,3,4);
P = zeros(length(T), 2);
h = 0.05;

for i = 1:length(T)
    [p,S,mu] = localCubicRegression(T(i), T', var, h);
    P(i, :) = p(3:-1:2);
end

% diff method
[~, dvar_P] = sbj1.calcFirstOrderDerivative(T, var, 'center');

% close all;
figure;
% position y
plot(T, var); hold on
plot(T, P(:,1));
legend('raw', 'cubic filter');

figure;
% velocity y
plot(T, dvar_P); hold on
plot(T, P(:,2)); 
legend('diff raw', 'cubic filter');

figure;
T_diff = 10*ones(size(T)) - T;
W = 1./sqrt(2*pi)*exp(-(T_diff.^2/(2*h^2)));
hold on;
plot(T_diff, W)

