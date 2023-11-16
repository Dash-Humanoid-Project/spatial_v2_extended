clear
disp('Unit Test for Tello Spanning Tree Model')

%% build spanning tree model

model = tello_spanning_tree_model()

%% set test joint

% Floating base joint for the first body (torso)
% * fb_q: [a, bi, cj, dk, x, y, z] of quaternion rotation (a+bi+cj+dk) and x,y,z position
% * fb_qd: [wx, wy, wz, vx, vy, vz]
% * fb_qdd: [wx_dot, wy_dot, wz_dot, ax, ay, az]
fb_q = [1 0 0 0 0 0 0]';
fb_qd = [0 0 0 0 0 0]';
fb_qdd = [0 0 0 0 0 0]';

% Spanning tree joint for Tello
% * each leg has (2 + 4 + 4) DoF
st_q = zeros(20,1);
st_qd = zeros(20,1);
st_qdd = zeros(20,1);

% Combine floating base joint and Tello spanning tree joints
q = [fb_q; zeros(20,1)];
q_dot = [fb_qd; zeros(20,1)];
q_ddot = [fb_qdd; zeros(20,1)];

%% Test inverse and forward dynamics

[tau, out] = ID(model, q, q_dot, q_ddot);   % Inverse dynamics
q_ddot_ABA = FDab(model, q, q_dot, tau);    % Forward dynamics

%% get G matrix

qr_hip_diff_1 = st_q(3:4);
ql_hip_diff_1 = st_q(5:6);
qr_ka_diff_1 = st_q(7:8);
ql_ka_diff_1 = st_q(9:10);
qr_hip_diff_2 = st_q(13:14);
ql_hip_diff_2 = st_q(15:16);
qr_ka_diff_2 = st_q(17:18);
ql_ka_diff_2 = st_q(19:20);

gear_ratio = 6;
hip_diff = derive_tello_hip_differential()
ka_diff = derive_tello_knee_ankle_differential()
G = zeros(20,10);
G(1:2,1) = [gear_ratio 1]';
G(3:6,2:3) = full(hip_diff.G(qr_hip_diff_1, ql_hip_diff_1));
G(7:10,4:5) = full(ka_diff.G(qr_ka_diff_1, ql_ka_diff_1));
G(11:12,6) = [gear_ratio 1]';
G(13:16,7:8) = full(hip_diff.G(qr_hip_diff_2, ql_hip_diff_2));
G(17:20,9:10) = full(ka_diff.G(qr_ka_diff_2, ql_ka_diff_2));

%% get implicit constraint for hip differentials and knee-ankle differentials

phi = zeros(4,2);
phi(1,:) = full(hip_diff.phi(qr_hip_diff_1, ql_hip_diff_1));
phi(2,:) = full(ka_diff.phi(qr_ka_diff_1, ql_ka_diff_1));
phi(3,:) = full(hip_diff.phi(qr_hip_diff_2, ql_hip_diff_1));
phi(4,:) = full(ka_diff.phi(qr_ka_diff_2, ql_ka_diff_2));
