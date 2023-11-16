clear
%% Test Tello Spanning Tree
disp('Unit Test for Tello Spanning Tree Model')

model = tello_spanning_tree_model()

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

% Test inverse and forward dynamics
[tau, out] = ID(model, q, q_dot, q_ddot);   % Inverse dynamics
q_ddot_ABA = FDab(model, q, q_dot, tau);    % Forward dynamics