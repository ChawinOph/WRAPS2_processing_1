%% Process WRAPS 2.0 Data collection
clc; clear; close all;
tic
%% Setup folder directory
sbj_folder_names = {'w001_2018-09-18'};
sbj_idcs = {'w001'};
sbj_names = {'WRAPS2_Chawin'};

load 'marker_pos_CAD.mat';

%% Setup sorted marker names for each segments (May vary for different subjects)
pelvis_sorted_marker_names = {'P_Top_R', 'P_Front_L', 'P_Top_L', 'P_Back_L', 'P_Hip_R', 'P_ASIS_R', 'P_ASIS_L', 'P_Hip_L'};
thorax_sorted_marker_names = {'T_Front_R', 'T_Front_L', 'T_Side_L', 'T_SternalNotch', 'T_C7', 'T_T8_Top', 'T_T8_L', 'T_T8_R', 'T_Xyphoid'};
head_sorted_marker_names = {'H_Ear_R', 'H_Forehead', 'H_Ear_L', 'H_Occiput', 'H_Vertex'};
sorted_segment_names = {'Pelvis','Thorax','Head'};
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
sbj1.importMarkerData_csv(trial_no, markers_static_max_flexion_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(trial_no, forcePlate_static_max_flexion_trials, forceplate_names, fplate_var_names);
sbj1.importMarkerData_csv(2, markers_static_max_flexion_trials, sorted_segment_names, sorted_marker_names);
sbj1.importForcePlateData_csv(2, forcePlate_static_max_flexion_trials, forceplate_names, fplate_var_names);
toc

figure;
sbj1.plotCoPvsTime(trial_no, 'Foot Plate');

figure;
sbj1.plotFvsTime(trial_no, 'Seat Plate');

figure;
sbj1.plotTrajCoP(trial_no, 'Seat Plate');

figure;
sbj1.plotTrajCoP(trial_no, 'Foot Plate');

%% Calculate the segment transformations

% calculate the thorax ring transformation w.r.t the vicoin frame in a trial (4 x 4 x no. of
% time steps using least meac square problem
n = find(strcmp({marker_pos_CAD.sbjID}, sbj1.sbj_id));

used_marker_indcs = zeros(length(marker_pos_CAD(n).t_brace_marker_names), 1);
for i = 1: length(marker_pos_CAD(n).t_brace_marker_names)
    % match the thorax brace CAD name with the subject trial data
    m = find(strcmp({sbj1.trial_data(trial_no).marker_data.segment_names}, 'Thorax'));
    used_marker_indcs(i) = find(strcmp(sbj1.trial_data(trial_no).marker_data(m).marker_names, ...
        marker_pos_CAD(n).t_brace_marker_names{i})); 
end

% store the 3d matrix of marker pos from the trial in the same order as
% the cad data
vicon_pos = sbj1.trial_data(trial_no).marker_data(m).marker_pos(:, :, used_marker_indcs);
vicon_pos = cat(2, vicon_pos, ones(size(vicon_pos, 1), 1, size(vicon_pos, 3)));
% append ones as an extra column
CAD_pos = marker_pos_CAD(n).t_brace_marker_pos;

% construct the an A matrix for the linear equation
A = zeros(4*length(CAD_pos), 16);
for n = 1 : length(CAD_pos)
    vec = [CAD_pos(n, :), 1];
    A(1 + 4*(n - 1), 1 : 4) = vec;
    A(2 + 4*(n - 1), 5 : 8) = vec;
    A(3 + 4*(n - 1), 9 : 12) = vec;
    A(4 + 4*(n - 1), 13 : 16) = vec;
end

for i = 1:length(vicon_pos)
    b = reshape(vicon_pos(i,:,:), [], 1);
    T = reshape(A\b, 4, 4)';
    % calculate the correct R by SVD
    [U,S,V] = svd(T(1:3, 1:3));
    T(1:3, 1:3) = (V*U')';
end

figure
plotCoordinatesTransform(T, 100); hold on;
sbj1.plotTrajCoP(trial_no, 'Seat Plate');
sbj1.plotTrajCoP(trial_no, 'Foot Plate');

scatter3(reshape(vicon_pos(1,1,:), 1, []), reshape(vicon_pos(1,2,:), 1, []), reshape(vicon_pos(1,3,:), 1, []))

sbj1.trial_data(trial_no).marker_data










