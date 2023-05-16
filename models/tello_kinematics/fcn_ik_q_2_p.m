function [p] = fcn_ik_q_2_p(q)
% inverse kinematics: joint angles (q) --> actuator angles (p)
% q, p are 5x1 vectors of angles

p = zeros(5,1);

% first actuation group 
% hip yaw (q1) - hip yaw actuator (p1)
p(1) = q(1);


% second actuation group - hip differential
% {hip roll (q2), hip pitch (q3)} - {actuator2 (p2), actuator3 (p3)}
a1 = 57/2500 - (8*cos(q(2))*sin(q(3)))/625 - (7*sin(q(2)))/625;
b1 = -(8*cos(q(3)))/625;
c1 = (49*cos(q(2)))/5000 + (399*sin(q(2)))/20000 + ...
    (57*cos(q(2))*sin(q(3)))/2500 - (7*sin(q(2))*sin(q(3)))/625 - 3021/160000;
d1 = sqrt(a1^2 + b1^2 - c1^2);
a2 = (7*sin(q(2)))/625 - (8*cos(q(2))*sin(q(3)))/625 + 57/2500;
b2 = -(8*cos(q(3)))/625;
c2 = (49*cos(q(2)))/5000 - (399*sin(q(2)))/20000 + ...
    (57*cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*sin(q(3)))/625 - 3021/160000;
d2 = sqrt(a2^2 + b2^2 - c2^2);

p(2) = atan2(-b1*d1 + a1*c1, b1*c1+a1*d1);
p(3) = atan2(-b2*d2 + a2*c2, b2*c2+a2*d2);


% third actuation group - knee-ankle differential
% {knee pitch (q4), ankle pitch (q5)} - {actuator4 (p4), actuator5 (p5)}
alpha_0 = (180-59.16)*pi/180;  % alpha_0
beta_0 = (180-66)*pi/180;   % beta_0
beta_ = q(5) + beta_0;
a3 = - (21*sin(pi/9))/6250 - (7*sin(beta_))/2500;
b3 = 13/625 - (7*cos(beta_))/2500 - (21*cos(pi/9))/6250;
c3 = (273*cos(pi/9))/12500 + (91*cos(beta_))/5000 - ...
    (147*cos(pi/9)*cos(beta_))/50000 - ...
    (147*sin(pi/9)*sin(beta_))/50000 - 163349/6250000;
d3 = sqrt(a3^2 + b3^2 - c3^2);
alpha_ = atan2( b3*d3 + a3*c3, b3*c3-a3*d3);

p(4) = q(4) - (alpha_ - alpha_0);
p(5) = q(4) + (alpha_ - alpha_0);

end








% p = zeros(5,1);
% 
%   p(1,1)=q(1);
%   p(2,1)=atan2((8*cos(q(3))*(((7*sin(q(2)))/625 + (8*cos(q(2))*sin(q(3)))/625 - 57/2500)^2 + (64*...
%          cos(q(3))^2)/390625 - ((49*cos(q(2)))/5000 + (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 - (7*sin(q(2))*...
%          sin(q(3)))/625 - 3021/160000)^2)^(1/2))/625 - ((7*sin(q(2)))/625 + (8*cos(q(2))*sin(q(3)))/625 - 57/2500)*((49*cos(q(2)))/5000 + (399*...
%          sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 - (7*sin(q(2))*sin(q(3)))/625 - 3021/160000), - (8*cos(q(3))*((49*...
%          cos(q(2)))/5000 + (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 - (7*sin(q(2))*...
%          sin(q(3)))/625 - 3021/160000))/625 - ((7*sin(q(2)))/625 + (8*cos(q(2))*sin(q(3)))/625 - 57/2500)*(((7*sin(q(2)))/625 + (8*cos(q(2))*...
%          sin(q(3)))/625 - 57/2500)^2 + (64*cos(q(3))^2)/390625 - ((49*cos(q(2)))/5000 + (399*sin(q(2)))/20000 + (57*cos(q(2))*...
%          sin(q(3)))/2500 - (7*sin(q(2))*sin(q(3)))/625 - 3021/160000)^2)^(1/2));
%   p(3,1)=atan2(((7*sin(q(2)))/625 - (8*cos(q(2))*sin(q(3)))/625 + 57/2500)*((49*...
%          cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*sin(q(3)))/625 - 3021/160000) + (8*...
%          cos(q(3))*(((7*sin(q(2)))/625 - (8*cos(q(2))*sin(q(3)))/625 + 57/2500)^2 + (64*cos(q(3))^2)/390625 - ((49*...
%          cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*...
%          sin(q(3)))/625 - 3021/160000)^2)^(1/2))/625, ((7*sin(q(2)))/625 - (8*cos(q(2))*sin(q(3)))/625 + 57/2500)*(((7*sin(q(2)))/625 - (8*cos(q(2))*...
%          sin(q(3)))/625 + 57/2500)^2 + (64*cos(q(3))^2)/390625 - ((49*cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*...
%          cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*sin(q(3)))/625 - 3021/160000)^2)^(1/2) - (8*cos(q(3))*((49*...
%          cos(q(2)))/5000 - (399*sin(q(2)))/20000 + (57*cos(q(2))*sin(q(3)))/2500 + (7*sin(q(2))*sin(q(3)))/625 - 3021/160000))/625);
%   p(4,1)=q(4) + (1007*pi)/1500 - atan2(((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*...
%          ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
%           (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000) - (((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
%          cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
%           (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
%           (21*sin(pi/9))/6250)^2)^(1/2)*((7*cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625), ((7*...
%          cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625)*((147*cos((19*pi)/30 + q(5))*...
%          cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 + (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 +...
%           163349/6250000) + ((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*(((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
%          cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
%           (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
%           (21*sin(pi/9))/6250)^2)^(1/2));
%   p(5,1)=q(4) - (1007*pi)/1500 + atan2(((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*...
%          ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
%           (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000) - (((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
%          cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
%           (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
%           (21*sin(pi/9))/6250)^2)^(1/2)*((7*cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625), ((7*...
%          cos((19*pi)/30 + q(5)))/2500 + (21*cos(pi/9))/6250 - 13/625)*((147*cos((19*pi)/30 + q(5))*...
%          cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 + (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 +...
%           163349/6250000) + ((7*sin((19*pi)/30 + q(5)))/2500 + (21*sin(pi/9))/6250)*(((7*cos((19*pi)/30 + q(5)))/2500 + (21*...
%          cos(pi/9))/6250 - 13/625)^2 - ((147*cos((19*pi)/30 + q(5))*cos(pi/9))/50000 - (273*cos(pi/9))/12500 - (91*cos((19*pi)/30 + q(5)))/5000 +...
%           (147*sin((19*pi)/30 + q(5))*sin(pi/9))/50000 + 163349/6250000)^2 + ((7*sin((19*pi)/30 + q(5)))/2500 +...
%           (21*sin(pi/9))/6250)^2)^(1/2));

 