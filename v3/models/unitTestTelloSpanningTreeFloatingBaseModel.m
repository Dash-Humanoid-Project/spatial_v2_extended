clear
disp('Unit Test for Tello Spanning Tree Model')

%% build spanning tree model

model = tello_spanning_tree_floating_base_model();

%% set test joint

% Spanning tree joint for Tello
% * torso has:
%   - 7 position coordinate [a, bi, cj, dk, x, y, z] of quaternion rotation
%      and x,y,z position
%   - 6 velocity coordinate [drx, dry, drz, dx, dy, dz]
%   - 6 acceleration coordinate [ddrx, ddry, ddrz, ddx, ddy, ddz]
% * each leg has (2 + 4 + 4) DoF
q_torso = [1 0 0 0 0 0 0]';
q = [q_torso; zeros(27,1)];
qd = zeros(26,1);
qdd = zeros(26,1);


%% [sanity check] test inverse and forward dynamics

[tau, out] = ID(model, q, qd, qdd);   % Inverse dynamics
qdd_ABA = FDab(model, q, qd, tau);    % Forward dynamics


%% set variable

q = sym('q',[27 1], 'real'); % spanning tree coordinate
y = sym('y',[17 1], 'real'); % minimal coordinate
syms g real;
global gear_ratio
gear_ratio = 6;

model_sym = tello_spanning_tree_floating_base_model(g);

%% Define y (minimal coordinate)
y(8) = q(9); % y8
y(9) = q(10) / gear_ratio; % y9
y(10) = q(11) / gear_ratio; % y10
y(11) = q(14) / gear_ratio; % y11
y(12) = q(15) / gear_ratio; % y12
y(13) = q(19); % y13
y(14) = q(20) / gear_ratio; % y14
y(15) = q(21) / gear_ratio; % y15
y(16) = q(24) / gear_ratio; % y16
y(17) = q(25) / gear_ratio; % y17

%% Get implicit position constraints (phi)

phi(1) = q(8) - gear_ratio * y(8);
[phi(2),phi(3),G_hd1] = derive_tello_hip_differential(y(9),y(10),q(12),q(13));
[phi(4),phi(5),G_kad1] = derive_tello_knee_ankle_differential(y(11),y(12),q(16),q(17));
phi(6) = q(18) - gear_ratio * y(13);
[phi(7),phi(8),G_hd2] = derive_tello_hip_differential(y(14),y(15),q(22),q(23));
[phi(9),phi(10),G_kad2] = derive_tello_knee_ankle_differential(y(16),y(17),q(26),q(27));

%% get G matrix

G = sym(zeros(26,16));
G(1:6,1:6) = diag(ones(6,1));
G(7:8,7) = [gear_ratio 1]';
G(9:12,8:9) = G_hd1;
G(13:16,10:11) = G_kad1;
G(17:18,12) = [gear_ratio 1]';
G(19:22,13:14) = G_hd2;
G(23:26,15:16) = G_kad2;

%% functions

% * Spanning tree coordinate (q) of hip differential and knee-ankle differentials
% can be divided into q = [qr1; qr2; ql1; ql2] where:
% - qr is the "independent" coordinate of the spanning tree coordinate
% - ql is the "dependent" coordinate of the spanning tree coordinate
% - qr is the pre-gearbox rotor angle
% - ql is the joint angle.
% For example, in the first hip differential, q = [q10, q11, q12, q13];
% * Minimal coordinate (y) of hip differential and knee-ankle differential
% refers to the post-gearbox rotor output angle

function [phi1,phi2,G] = derive_tello_hip_differential(y1,y2,ql1,ql2)
  phi1 = ...
        (57*sin(y1))/2500 - (49*cos(ql1))/5000 - (399*sin(ql1))/20000 - ...
        (8*cos(y1)*cos(ql2))/625 - (57*cos(ql1)*sin(ql2))/2500 - ...
        (7*sin(y1)*sin(ql1))/625 + (7*sin(ql1)*sin(ql2))/625 - ...
        (8*cos(ql1)*sin(y1)*sin(ql2))/625 + 3021/160000;
  
  phi2 = ...
        (57*sin(y2))/2500 - (49*cos(ql1))/5000 + (399*sin(ql1))/20000 - ...
        (8*cos(y2)*cos(ql2))/625 - (57*cos(ql1)*sin(ql2))/2500 + ...
        (7*sin(y2)*sin(ql1))/625 - (7*sin(ql1)*sin(ql2))/625 - ...
        (8*cos(ql1)*sin(y2)*sin(ql2))/625 + 3021/160000;

  global gear_ratio
  phi = [phi1; phi2];
  qr1 = gear_ratio * y1;
  qr2 = gear_ratio * y2;
  qr = [qr1; qr2];
  ql = [ql1; ql2];
  Ki = jacobian(phi,qr);
  Kd = jacobian(phi,ql);
  J_dqr_2_dql = -(Kd\Ki);
  J_dy_2_dql = gear_ratio * J_dqr_2_dql;
  G = [gear_ratio 0;...
       0 gear_ratio;
       J_dy_2_dql];
end

function [phi1,phi2,G] = derive_tello_knee_ankle_differential(y1,y2,ql1,ql2)
  phi1 = ...
        (21*cos(y1/2 - y2/2 + (1979*pi)/4500))/6250 - ...
        (13*cos(y1/2 - y2/2 + (493*pi)/1500))/625 - (273*cos(pi/9))/12500 - ...
        (7*sin(y1/2 - y2/2 + ql2 + (231*pi)/500))/2500 + ...
        (91*sin(ql2 + (2*pi)/15))/5000 - (147*sin(ql2 + pi/45))/50000 + ...
        163349/6250000;
  
  phi2 = ql1 - y2/2 - y1/2;

  global gear_ratio
  phi = [phi1; phi2];
  qr1 = gear_ratio * y1;
  qr2 = gear_ratio * y2;
  qr = [qr1; qr2];
  ql = [ql1; ql2];
  Ki = jacobian(phi,qr);
  Kd = jacobian(phi,ql);
  J_dqr_2_dql = -(Kd\Ki);
  J_dy_2_dql = gear_ratio * J_dqr_2_dql;
  G = [gear_ratio 0;...
       0 gear_ratio;
       J_dy_2_dql];
end