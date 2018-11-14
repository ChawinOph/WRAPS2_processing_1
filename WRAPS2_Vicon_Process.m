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
legend('v_x,cent5pt', 'v_y,cent5pt', 'v_z,cent5pt')
% plot(time, dvar_bwd, ':r'); 
% plot(time, dvar_cen, '--b'); 

%% smooth the marker pos on the brace before finding the velocities of each marker and IARs
% [1] A. Page, P. Candelas, and F. Belmar, “On the use of local fitting
% techniques for the analysis of physical dynamic systems,” Eur. J. Phys.,
% vol. 27, no. 2, p. 273, 2006. [2] A. Page, H. de Rosario, V. Mata, R.
% Porcar, J. Solaz, and M. J. Such, “Kinematics of the trunk in sitting
% posture: An analysis based on the instantaneous axis of rotation,”
% Ergonomics, vol. 52, no. 6, pp. 695–706, Jun. 2009.
% pelvis_brace_marker_pos = sbj1.sbj_WRAPS2(trial_no).trial_transform_data(1).marker_pos;
% thorax_brace_marker_pos = sbj1.sbj_WRAPS2(trial_no).trial_transform_data(2).marker_pos;
% T = 0: 1/sbj1.freq_marker : (length(pelvis_brace_marker_pos) - 1)/sbj1.freq_marker;
% var = pelvis_brace_marker_pos(:,3,4);
% P = zeros(length(T), 2);
% h = 0.05;
% 
% for i = 1:length(T)
%     [p,S,mu] = localCubicRegression(T(i), T', var, h);
%     P(i, :) = p(3:-1:2);
% end
% 
% % diff method
% [~, dvar_P] = sbj1.calcFirstOrderDerivative(T, var, 'center');
% 
% % close all;
% figure;
% % position y
% plot(T, var); hold on
% plot(T, P(:,1));
% legend('raw', 'cubic filter');
% 
% [time, dvar_cen5] = sbj1.calcFirstOrderDerivative(time, P(:,1), 'center_5point');
% figure;
% % velocity y
% plot(T, dvar_P); hold on
% plot(T, dvar_cen5); 
% legend('diff raw', 'cubic filter')
% 
% 
% figure;
% T_diff = 10*ones(size(T)) - T;
% W = 1./sqrt(2*pi)*exp(-(T_diff.^2/(2*h^2)));
% hold on;
% plot(T_diff, W)

%% Try other filtering techniques
% use point from thorax rings
t_marker_pos_filt = sbj1.raw_data(trial_no).marker_data(2).marker_pos;
t_marker_pos_raw = sbj1.raw_data(trial_no).marker_data(2).marker_pos_raw;
T = 0: 1/sbj1.freq_marker : (length(t_marker_pos_filt) - 1)/sbj1.freq_marker;

marker_no = 1;

var_raw = t_marker_pos_raw(:,:,marker_no); % raw
var_filt = t_marker_pos_filt(:,:,marker_no); % butter 4th order lowpass 6 Hz

% take the first derivative
[T, dvar_filt_6hzbutter] = sbj1.calcFirstOrderDerivative(T, var_filt, 'center_5point');
[T, dvar_raw] = sbj1.calcFirstOrderDerivative(T, var_raw, 'center_5point');
[T, ddvar_filt_6hzbutter] = sbj1.calcSecondOrderDerivative(T, var_filt, 'center_5point');
[T, ddvar_raw] = sbj1.calcSecondOrderDerivative(T, var_raw, 'center_5point');

% [1] F. J. Alonso, J. M. Del Castillo, and P. Pintado, “An Automatic
% Filtering Procedure for Processing Biomechanical Kinematic Signals,”
% Springer, Berlin, Heidelberg, 2004, pp. 281–291. 
% [2] R. Aissaoui, S. Husse, H. Mecheri, G. Parent, and J. a. D. Guise,
% “Automatic filtering techniques for three-dimensional kinematics data
% using 3D motion capture system,” in 2006 IEEE International Symposium on
% Industrial Electronics, 2006, vol. 1, pp. 614–619. 
% (Golyandina et al., 2001, Chapter 6 has the suggestion for choosing window length).

% Step 1 Embedding
var_raw_1d = var_raw(:, 2); % use y
N = length(var_raw_1d); % signal length
L = round(N/60); % signal length (suggested round(N/60))
%  construct the Hankel matrix (for 1D data), the element in i+j = constant
%  are equal (somtimes it referred to as the trajectory matrix)
X = hankel(var_raw_1d(1:L), var_raw_1d(L: end)); % (size L x N - L + 1)

% Step 2 SVD (performs a singular value decomposition of matrix A, such that A = U*Sigma*V'.)
% S = X*X'; [eigvec, eig_diag] = eig(S);
[U,Sigma,V] = svd(X); 

% Step 3 Grouping
r = 10; % number of first elementary matrices used (r <= L)
X_est = U(:, 1:r)*Sigma(1:r, 1:r)*V(:, 1:r)';

% Step 4: Reconstruction (Diagonal Averaging)


%% plot all results
figure;
title_font_size = 15;
% plot v_x, v_y, v_z
raw_v_alpha = 0.25;
subplot(3,2,1)
p_raw_x = plot(T, dvar_raw(:, 1), 'k'); hold on; 
plot(T, dvar_filt_6hzbutter(:, 1), 'r');
title('$v_x$', 'Interpreter', 'latex', 'fontsize', title_font_size)
p_raw_x.Color(4) = 0.25;
xlim([0 T(end)]);
ylabel('$mm/s$', 'Interpreter', 'latex');
grid on; grid minor;

subplot(3,2,3)
p_raw_y = plot(T, dvar_raw(:, 2), 'k'); hold on; 
plot(T, dvar_filt_6hzbutter(:, 2), 'g');
title('$v_y$', 'Interpreter', 'latex', 'fontsize', title_font_size)
p_raw_y.Color(4) = 0.25;
xlim([0 T(end)]);
ylabel('$mm/s$', 'Interpreter', 'latex');
grid on;  grid minor;

subplot(3,2,5)
p_raw_z = plot(T, dvar_raw(:, 3), 'k'); hold on; 
plot(T, dvar_filt_6hzbutter(:, 3), 'b');
title('$v_z$', 'Interpreter', 'latex', 'fontsize', title_font_size)
p_raw_z.Color(4) = 0.25;
xlim([0 T(end)]);
xlabel('time (s)', 'Interpreter', 'latex')
ylabel('$mm/s$', 'Interpreter', 'latex');
grid on;  grid minor;

% plot a_x, a_y, a_z
raw_a_alpha = 0.1;
raw_filt_ylim_ratio = 20;

subplot(3,2,2)
p_raw_x = plot(T, ddvar_raw(:, 1), 'k'); hold on; 
plot(T, ddvar_filt_6hzbutter(:, 1), 'r');
title('$a_x$', 'Interpreter', 'latex', 'fontsize', title_font_size);
p_raw_x.Color(4) = raw_a_alpha;
xlim([0 T(end)]);
y_lim = ylim();
ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_ratio); % symmetric ylim
ylabel('$mm/s^2$', 'Interpreter', 'latex');
grid on;  grid minor;

subplot(3,2,4)
p_raw_y = plot(T, ddvar_raw(:, 2), 'k'); hold on; 
plot(T, ddvar_filt_6hzbutter(:, 2), 'g');
title('$a_y$', 'Interpreter', 'latex', 'fontsize', title_font_size);
p_raw_y.Color(4) = raw_a_alpha;
xlim([0 T(end)]);
y_lim = ylim();
ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_ratio);
ylabel('$mm/s^2$', 'Interpreter', 'latex')
grid on;  grid minor;

subplot(3,2,6)
p_raw_z = plot(T, ddvar_raw(:, 3), 'k'); hold on; 
plot(T, ddvar_filt_6hzbutter(:, 3), 'b');
title('$a_z$', 'Interpreter', 'latex', 'fontsize', title_font_size);
p_raw_z.Color(4) = raw_a_alpha;
xlim([0 T(end)]);
y_lim = ylim();
ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_ratio);
xlabel('time (s)', 'Interpreter', 'latex');
ylabel('$mm/s^2$', 'Interpreter', 'latex');
grid on;  grid minor;



