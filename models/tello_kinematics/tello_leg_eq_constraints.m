function h = tello_leg_eq_constraints
% System of constraints that defines relation between two spaces of angles
% h_i(p, q) = 0, i=1,2,3,4,5
% p: vector of actuator angles 5x1 (after reduction)
% q: vector of joint angles 5x1
% NOTE: h2, h3, h4, h5 may require additional term - constant angle offset.


syms p q [5 1] real
h = sym(zeros(5,1));

h(1) = p1-q1;

h(2) = ...
    (57*sin(p2))/2500 - (49*cos(q2))/5000 - (399*sin(q2))/20000 - ...
    (8*cos(p2)*cos(q3))/625 - (57*cos(q2)*sin(q3))/2500 - ...
    (7*sin(p2)*sin(q2))/625 + (7*sin(q2)*sin(q3))/625 - ...
    (8*cos(q2)*sin(p2)*sin(q3))/625 + 3021/160000;

h(3) = ...
    (57*sin(p3))/2500 - (49*cos(q2))/5000 + (399*sin(q2))/20000 - ...
    (8*cos(p3)*cos(q3))/625 - (57*cos(q2)*sin(q3))/2500 + ...
    (7*sin(p3)*sin(q2))/625 - (7*sin(q2)*sin(q3))/625 - ...
    (8*cos(q2)*sin(p3)*sin(q3))/625 + 3021/160000;

h(4) = ...
    (21*cos(p4/2 - p5/2 + (1979*pi)/4500))/6250 - ...
    (13*cos(p4/2 - p5/2 + (493*pi)/1500))/625 - (273*cos(pi/9))/12500 - ...
    (7*sin(p4/2 - p5/2 + q5 + (231*pi)/500))/2500 + ...
    (91*sin(q5 + (2*pi)/15))/5000 - (147*sin(q5 + pi/45))/50000 + ...
    163349/6250000;

h(5) = q4 - p5/2 - p4/2;

 


