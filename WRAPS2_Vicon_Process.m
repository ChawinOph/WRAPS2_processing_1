%% Process WRAPS 2.0 Data collection 
clc; clear; close all;
tic
%% Setup folder directory
sbj_folder_names = {'w001_2018-09-18'};
sbj_idcs = {'w001'};
sbj_names = {'WRAPS2_Chawin'};

load 'marker_pos_CAD.mat';

%% Setup sorted marker names for each segments (May subject to cahnge for each subject)
pelvis_sorted_marker_names = {'P_Top_R', 'P_Front_L', 'P_Top_L', 'P_Back_L', 'P_Hip_R', 'P_ASIS_R', 'P_ASIS_L', 'P_Hip_L'};
thorax_sorted_marker_names = {'T_Front_R', 'T_Front_L', 'T_Side_L', 'T_SternalNotch', 'T_C7', 'T_T8_Top', 'T_T8_L', 'T_T8_R', 'T_Xyphoid'};
head_sorted_marker_names = {'H_Ear_R', 'H_Forehead', 'H_Ear_L', 'H_Occiput', 'H_Vertex'};
sorted_marker_names = {pelvis_sorted_marker_names; thorax_sorted_marker_names; head_sorted_marker_names}; % segment no 1,2,3

%% Setup Force Plate names
forceplate_names = {'Seat Plate', 'Foot Plate'};
fplate_var_names = {'Force', 'CoP'};
% Name list in the stucture output will be  Seat Plate - Force, Seat Plate - CoP, Foot Plate - ...

%% Setup VICON markers trial file names by the order of subjects
%# Markers_static_baseline = {''};
markers_static_max_flexion_trials = {'Static Max Flexion03 Markers'}; % trial no. 1

%% Setup VICON force plate trial file names by the order of subjects
% forceplate_static_baseline = {'Static Max Flexion03 ForcePlate'};
forcePlate_static_max_flexion_trials = {'Static Max Flexion03 ForcePlate'}; % trial no. 1

%% Import Vicon marker data
sbj_index = 1; trial_no = 1;
sbj1 = Subject(sbj_index, sbj_folder_names, sbj_idcs, sbj_names);

%% import data
sbj1.importMarkerData_csv(trial_no, markers_static_max_flexion_trials, sorted_marker_names);
sbj1.ImportForcePlateData_csv(trial_no, forcePlate_static_max_flexion_trials, forceplate_names, fplate_var_names);
toc

figure;
sbj1.plotFvsTime(trial_no, 'Seat Plate'); 
title('Seat Plate')

figure;
sbj1.plotFvsTime(trial_no, 'Foot Plate'); 
title('Foot Plate')

figure;
sbj1.plotTrajCoP(trial_no, 'Seat Plate'); 
title('Seat Plate')




