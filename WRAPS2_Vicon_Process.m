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

markers_fast_flexion_extension_trials = {'Slow Flexion-Extension Markers'}; % trial no.6

%% Setup VICON force plate trial file names by the order of subjects
% forceplate_static_baseline = {'Static Max Flexion03 ForcePlate'};
force_plate_static_max_flexion_trials = {'Static Max Flexion03 ForcePlate'}; % trial no. 1

force_plate_reach_max_knee_ipsi_trials = {'Max Reach Knee Ipsi ForcePlate'}; % trial no. 2

force_plate_static_max_rbend_trials = {'Static Max R-Side Bend01 ForcePlate'}; % trial no. 3

force_plate_reach_max_shoulder_center_trials = {'Max Reach Shoulder Center01 ForcePlate'}; % trial no. 4

force_plate_reach_max_shoulder_contra_trials = {'Max Reach Shoulder Contra01 ForcePlate'}; % trial no. 5

force_plate_fast_flexion_extension_trials = {'Slow Flexion-Extension ForcePlate'}; % trial no.6

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


% Calculate transformation of the WRAPS2.0 from vicon
toc
% sbj1.vizTrial(1, 1000)
% sbj1.vizTrial(6, 1500)

%% Calculate CM of WRAPS using CAD and scale


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
% trial_no = 1;
% t_marker_pos_filt = sbj1.raw_data(trial_no).marker_data(2).marker_pos;
% t_marker_pos_raw = sbj1.raw_data(trial_no).marker_data(2).marker_pos_raw;
% T = 0: 1/sbj1.freq_marker : (length(t_marker_pos_filt) - 1)/sbj1.freq_marker;
% 
% marker_no = 1;
% 
% var_raw = t_marker_pos_raw(:,:,marker_no); % raw
% var_filt_6hzbutter = t_marker_pos_filt(:,:,marker_no); % butter 4th order lowpass 6 Hz
% 
% % take the first derivative
% dvar_filt_6hzbutter = sbj1.calcFirstOrderDerivative(var_filt_6hzbutter, 'center_5point');
% dvar_raw = sbj1.calcFirstOrderDerivative(var_raw, 'center_5point');
% ddvar_filt_6hzbutter = sbj1.calcSecondOrderDerivative(var_filt_6hzbutter, 'center_5point');
% ddvar_raw = sbj1.calcSecondOrderDerivative(var_raw, 'center_5point');
% 
% % [1] F. J. Alonso, J. M. Del Castillo, and P. Pintado, “An Automatic
% % Filtering Procedure for Processing Biomechanical Kinematic Signals,”
% % Springer, Berlin, Heidelberg, 2004, pp. 281–291. 
% %
% % [2] R. Aissaoui, S.
% % Husse, H. Mecheri, G. Parent, and J. a. D. Guise, “Automatic filtering
% % techniques for three-dimensional kinematics data using 3D motion capture
% % system,” in 2006 IEEE International Symposium on Industrial Electronics,
% % 2006, vol. 1, pp. 614–619. (Golyandina et al., 2001, Chapter 6 has the
% % suggestion for choosing window length). 
% % [1] N. Golyandina, V. Nekrutkin, and A. A. Zhigljavsky, Analysis of time
% % series structure: SSA and related techniques. Boca Raton, Florida:
% % Chapman and Hall/CRC, 2001.
% % 
% % [pxx,w] = periodogram(var_raw(:, 1));
% % figure; plot(w,pxx);
% 
% var_filt_SSA = zeros(size(var_raw));
% for dof =1:3 % y component
%     var_raw_1d = var_raw(:, dof); % use y
%     
%     g_prev = var_raw_1d;
%     ddotg_prev = sbj1.calcSecondOrderDerivative(g_prev, 'center_5point');
%     rms_ddotg_prev = rms(ddotg_prev);
%     
%     % recursively apply the SSA to the time series until change of the rms is
%     % smaller than 1% of the previous acceleration rms
%     N_max_iter = 1;
%     tic
%     for n = 1: N_max_iter
%         % Step 1 Embedding
%         N = length(g_prev); % signal length
% %         L = round(N/60); % signal length (suggested round(N/60))
%         L = 30;
%         %  construct the Hankel matrix (for 1D data), the element in i+j = constant
%         %  are equal (somtimes it referred to as the trajectory matrix)
%         X = hankel(g_prev(1:L), g_prev(L: end)); % (size L x N - L + 1)
%         
%         % Step 2 SVD (performs a singular value decomposition of matrix A, such that A = U*Sigma*V'.)
%         % S = X*X'; [eigvec, eig_diag] = eig(S); % cross check
%         % diag(sqrt(sort(diag((eig_diag)), 'descend'))) % cross check
%         [U, Sigma, V] = svd(X);
%         
%         % Step 3 Grouping (Eigentriple grouping) grouping with eigenvalues
%         % of S = X*X' that contribute to 99.999% of the sum
% %         eig_sum_threshold =  99.999/100*trace((Sigma(1:size(Sigma, 1), 1:size(Sigma, 1))));
% %         for r = 1:size(Sigma, 1) % number of first elementary matrices used (4 < r <= L)
% %             disp(['r = ', num2str(r)])
% %             disp(100*trace((Sigma(1:r + 1, 1:r + 1)))/trace((Sigma(1:size(Sigma, 1), 1:size(Sigma, 1)))));
% %             if trace((Sigma(1:r + 1, 1:r + 1))) > eig_sum_threshold
% %                 break;
% %             end
% %         end
%         r = 1;
% 
%         % rank truncation for the approximation of X
%         Y = U(:, 1:r)*Sigma(1:r, 1:r)*V(:, 1:r)';
%         
%         % Step 4: Reconstruction (Diagonal Averaging)
%         g_curr = sbj1.calcDiagonalAverage_SSA(Y);
%         ddotg_curr = sbj1.calcSecondOrderDerivative(g_curr, 'center_5point');
%         rms_ddotg_curr = rms(ddotg_curr);
%         
%         % check the manitude chage of the rms
% %         if abs(rms_ddotg_curr - rms_ddotg_prev) < 0.01*rms_ddotg_prev
% %             break;
% %         else
% %             rms_ddotg_prev = rms_ddotg_curr;
% %             g_prev = g_curr;
% %         end
%     end
%     n
%     toc
%     var_filt_SSA(:, dof) = g_curr;
% end
% 
% dvar_filt_SSA = sbj1.calcFirstOrderDerivative(var_filt_SSA, 'center_5point');
% ddvar_filt_SSA = sbj1.calcSecondOrderDerivative(var_filt_SSA, 'center_5point');
% 
% % plot all results
% close all;
% figure;
% title_font_size = 15;
% 
% % plot displacement
% raw_filt_ylim_disp_ratio = 1;
% 
% subplot(3,3,1)
% p_raw_x = plot(T, var_raw(:, 1), 'k'); hold on; 
% plot(T, var_filt_6hzbutter(:, 1), 'b');
% plot(T, var_filt_SSA(:, 1), 'r');
% title('$x$', 'Interpreter', 'latex', 'fontsize', title_font_size)
% p_raw_x.Color(4) = 1;
% xlim([0 T(end)]);
% % y_lim = ylim();
% % ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_disp_ratio);
% ylabel('$mm$', 'Interpreter', 'latex');
% grid on; grid minor;
% 
% subplot(3,3,2)
% p_raw_x = plot(T, var_raw(:, 2), 'k'); hold on; 
% plot(T, var_filt_6hzbutter(:, 2), 'b');
% plot(T, var_filt_SSA(:, 2), 'r');
% title('$y$', 'Interpreter', 'latex', 'fontsize', title_font_size)
% p_raw_x.Color(4) = 1;
% xlim([0 T(end)]);
% % y_lim = ylim();
% % ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_disp_ratio);
% ylabel('$mm$', 'Interpreter', 'latex');
% grid on; grid minor;
% 
% subplot(3,3,3)
% p_raw_x = plot(T, var_raw(:, 3), 'k'); hold on; 
% plot(T, var_filt_6hzbutter(:, 3), 'b');
% plot(T, var_filt_SSA(:, 3), 'r');
% title('$z$', 'Interpreter', 'latex', 'fontsize', title_font_size)
% p_raw_x.Color(4) = 1;
% xlim([0 T(end)]);
% % y_lim = ylim();
% % ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_disp_ratio);
% ylabel('$mm$', 'Interpreter', 'latex');
% grid on; grid minor;
% 
% % plot v_x, v_y, v_z
% raw_filt_ylim_v_ratio = 1;
% raw_v_alpha = 0.25;
% 
% subplot(3,3,4)
% p_raw_x = plot(T, dvar_raw(:, 1), 'k'); hold on; 
% plot(T, dvar_filt_6hzbutter(:, 1), 'b');
% plot(T, dvar_filt_SSA(:, 1), 'r');
% title('$v_x$', 'Interpreter', 'latex', 'fontsize', title_font_size)
% p_raw_x.Color(4) = 0.25;
% xlim([0 T(end)]);
% y_lim = ylim();
% ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_v_ratio);
% ylabel('$mm/s$', 'Interpreter', 'latex');
% grid on; grid minor;
% 
% subplot(3,3,5)
% p_raw_y = plot(T, dvar_raw(:, 2), 'k'); hold on; 
% plot(T, dvar_filt_6hzbutter(:, 2), 'b');
% plot(T, dvar_filt_SSA(:, 2), 'r');
% title('$v_y$', 'Interpreter', 'latex', 'fontsize', title_font_size)
% p_raw_y.Color(4) = 0.25;
% xlim([0 T(end)]);
% y_lim = ylim();
% ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_v_ratio);
% ylabel('$mm/s$', 'Interpreter', 'latex');
% grid on;  grid minor;
% 
% subplot(3,3,6)
% p_raw_z = plot(T, dvar_raw(:, 3), 'k'); hold on; 
% plot(T, dvar_filt_6hzbutter(:, 3), 'b');
% plot(T, dvar_filt_SSA(:, 3), 'r');
% title('$v_z$', 'Interpreter', 'latex', 'fontsize', title_font_size)
% p_raw_z.Color(4) = 0.25;
% xlim([0 T(end)]);
% y_lim = ylim();
% ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_v_ratio);
% xlabel('time (s)', 'Interpreter', 'latex')
% ylabel('$mm/s$', 'Interpreter', 'latex');
% grid on;  grid minor;
% 
% % plot a_x, a_y, a_z
% raw_a_alpha = 0.1;
% raw_filt_ylim_a_ratio = 30;
% 
% subplot(3,3,7)
% p_raw_x = plot(T, ddvar_raw(:, 1), 'k'); hold on; 
% plot(T, ddvar_filt_6hzbutter(:, 1), 'b');
% plot(T, ddvar_filt_SSA(:, 1), 'r');
% title('$a_x$', 'Interpreter', 'latex', 'fontsize', title_font_size);
% p_raw_x.Color(4) = raw_a_alpha;
% xlim([0 T(end)]);
% y_lim = ylim();
% ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_a_ratio); % symmetric ylim
% ylabel('$mm/s^2$', 'Interpreter', 'latex');
% grid on;  grid minor;
% 
% subplot(3,3,8)
% p_raw_y = plot(T, ddvar_raw(:, 2), 'k'); hold on; 
% plot(T, ddvar_filt_6hzbutter(:, 2), 'b');
% plot(T, ddvar_filt_SSA(:, 2), 'r');
% title('$a_y$', 'Interpreter', 'latex', 'fontsize', title_font_size);
% p_raw_y.Color(4) = raw_a_alpha;
% xlim([0 T(end)]);
% y_lim = ylim();
% ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_a_ratio);
% ylabel('$mm/s^2$', 'Interpreter', 'latex')
% grid on;  grid minor;
% 
% subplot(3,3,9)
% p_raw_z = plot(T, ddvar_raw(:, 3), 'k'); hold on; 
% plot(T, ddvar_filt_6hzbutter(:, 3), 'b');
% plot(T, ddvar_filt_SSA(:, 3), 'r');
% title('$a_z$', 'Interpreter', 'latex', 'fontsize', title_font_size);
% p_raw_z.Color(4) = raw_a_alpha;
% xlim([0 T(end)]);
% y_lim = ylim();
% ylim([-max(abs(y_lim)), max(abs(y_lim))]/raw_filt_ylim_a_ratio);
% xlabel('time (s)', 'Interpreter', 'latex');
% ylabel('$mm/s^2$', 'Interpreter', 'latex');
% grid on;  grid minor;
% 
% leg = legend('Raw', 'Butterworth (6 Hz)', 'Recursive SSA');
% leg.Interpreter =  'latex';


%% construct the screw axis (will be inserted into the importMarkerData function)
trial_no= 3;

[v_g_pelv, omega_pelv, ISA_pelv, centroid_pelv, theta_pelv, T,  T_used_indcs_pelv] = sbj1.calcISA(trial_no, 'Pelvis Brace');
[v_g_thor, omega_thor, ISA_thor, centroid_thor, theta_thor, ~,  T_used_indcs_thor] = sbj1.calcISA(trial_no, 'Thorax Brace');


sbj1.vizTrial(trial_no, 1000);
scatter3(centroid_thor(:,1), centroid_thor(:,2), centroid_thor(:,3));  hold on;
scatter3(centroid_pelv(:,1), centroid_pelv(:,2), centroid_pelv(:,3));

n_down_sample_ratio = 100;
Tt = T_used_indcs_thor(1:round(length(T_used_indcs_thor)/n_down_sample_ratio):end);
quiver3(ISA_thor(Tt,1), ISA_thor(Tt,2), ISA_thor(Tt,3), ...
    omega_thor(Tt,1), omega_thor(Tt,2), omega_thor(Tt,3), 10); hold on;

quiver3(centroid_thor(Tt,1), centroid_thor(Tt,2), centroid_thor(Tt,3), v_g_thor(Tt,1), v_g_thor(Tt,2), v_g_thor(Tt,3), 10);

Tp = T_used_indcs_pelv(1:round(length(T_used_indcs_thor)/n_down_sample_ratio):end);
quiver3(ISA_pelv(Tp,1), ISA_pelv(Tp,2), ISA_pelv(Tp,3),...
    omega_pelv(Tp,1), omega_pelv(Tp,2), omega_pelv(Tp,3), 10);
quiver3(centroid_pelv(Tp,1), centroid_pelv(Tp,2), centroid_pelv(Tp,3),v_g_pelv(Tp,1), v_g_pelv(Tp,2), v_g_pelv(Tp,3), 10);
axis equal;

% calculate the relative screw axis between pelvis and thorax

% get the psis positions
indx_r_psis = strcmp(sbj1.sbj_anthro(trial_no).torso_landmark_names, 'R_PSIS');
indx_l_psis = strcmp(sbj1.sbj_anthro(trial_no).torso_landmark_names, 'L_PSIS');
pos_r_psis = sbj1.sbj_anthro(trial_no).landmark_pos(:, :, indx_r_psis);
pos_l_psis = sbj1.sbj_anthro(trial_no).landmark_pos(:, :, indx_l_psis);
pos_m_psis = (pos_r_psis + pos_l_psis)/2;

pelv_g2m_psis = pos_m_psis - centroid_pelv;
thor_g2m_psis = pos_m_psis - centroid_thor;

v_m_psis_on_pelv = v_g_pelv + cross(omega_pelv, pelv_g2m_psis, 2);
v_m_psis_on_thor = v_g_thor + cross(omega_thor, thor_g2m_psis, 2);
v_m_psis_r = v_m_psis_on_thor - v_m_psis_on_pelv;

omega_r = omega_thor - omega_pelv;
omega_r_norm = vecnorm(omega_r, 2, 2);

PH = (cross(omega_r, v_m_psis_r, 2))./omega_r_norm.^2;
ISA_r =  pos_m_psis + PH;

% calculate the finite angle of rotation
Tr_plot = 0: 1/sbj1.freq_marker : (length(pos_m_psis) - 1)/sbj1.freq_marker;
theta_r = cumtrapz(Tr_plot, omega_r);

T_used_indcs_r = find(omega_r_norm >= 0.25*max(omega_r_norm));
 
scatter3(pos_m_psis(:,1), pos_m_psis(:,2), pos_m_psis(:,3));
Tr = T_used_indcs_r(1:round(length(T_used_indcs_r)/n_down_sample_ratio):end);
quiver3(ISA_r(Tr,1), ISA_r(Tr,2), ISA_r(Tr,3),...
    omega_r(Tr,1), omega_r(Tr,2), omega_r(Tr,3), 10);
quiver3(pos_m_psis(Tr,1), pos_m_psis(Tr,2), pos_m_psis(Tr,3), v_m_psis_r(Tr,1), v_m_psis_r(Tr,2), v_m_psis_r(Tr,3), 10);
axis equal;

figure;
plot(T, omega_pelv); hold on;
plot(T, omega_thor, ':');  
plot(T, omega_r, '--'); 
legend('\omega_{pelv,x}','\omega_{pelv,y}','\omega_{pelv,z}', ...
    '\omega_{thor,x}','\omega_{thor,y}','\omega_{thor,z}', ...
    '\omega_{r,x}','\omega_{r,y}','\omega_{r,z}');

figure;
plot(T, theta_pelv); hold on;
plot(T, theta_thor, ':'); 
plot(T, theta_r, '--');
legend('\theta_{pelv,x}','\theta_{pelv,y}','\theta_{pelv,z}',...
    '\theta_{thor,x}','\theta_{thor,y}','\theta_{thor,z}',...
    '\theta_{r,x}','\theta_{r,y}','\theta_{r,z}');

figure;
plot(T, v_g_pelv); hold on;
plot(T, v_g_thor, ':'); 
plot(T, v_m_psis_r, '--');
legend('v_{pelv,x}','v_{pelv,y}','v_{pelv,z}', ...
       'v_{thor,x}','v_{thor,y}','v_{thor,z}', ...
       'v_{r,x}','v_{r,y}','v_{r,z}');
   
% find ISA_r with respect to the pelvis brace moving frame
cluster_no = strcmp({sbj1.sbj_WRAPS2(trial_no).trial_transform_data.cluster_name}, 'Pelvis Brace');
T_v2pelv_brace = sbj1.sbj_WRAPS2(trial_no).trial_transform_data(cluster_no).transforms_vicon2seg;
ISA_r_wrt_pelvis_brace = sbj1.calcPosInNewFrame(T_v2pelv_brace, ISA_r);

figure;
quiver3(ISA_r(Tr,1), ISA_r(Tr,2), ISA_r(Tr,3),...
    omega_r(Tr,1), omega_r(Tr,2), omega_r(Tr,3), 10); hold on;
quiver3(ISA_r_wrt_pelvis_brace(Tr,1), ISA_r_wrt_pelvis_brace(Tr,2), ISA_r_wrt_pelvis_brace(Tr,3),...
    omega_r(Tr,1), omega_r(Tr,2), omega_r(Tr,3), 10);
axis equal;






