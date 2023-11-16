function model = tello_spanning_tree_model(gravity)

import casadi.*

if nargin == 1
    model.gravity = [0; 0; gravity];
end

model = RBD_model();

model.NB = 21;
model.parent = [0, ...
                1, 1, 3, 3, 3, 6, 7, 7, 7, 10, ...
                1, 1, 13, 13, 13, 16, 17, 17 ,17, 20];

% import tello inertial parameters from solidworks-generated data files
inertial = fcn_import_params_solidworks('tello_params/inertial_*.txt');

% 1st body: torso
model.joint{1} = floatingBaseJoint();
model.Xtree{1} = eye(6);
model.I{1} = mcI(inertial.mass{1},...
                 inertial.CoM{1},...
                 inertial.Ic{1});

rotor_mass = 0.07;
rotor_CoM = [0 0 0]';
rotor_Ixx = 2.5984e-05;
rotor_Iyy = rotor_Ixx;
rotor_Izz = 5.1512e-05;
I_rotor = [rotor_Ixx rotor_Iyy rotor_Izz]' .* eye(3);
rotorInertia = mcI(rotor_mass, rotor_CoM, I_rotor);

% joint orientations
R_top = eye(3);
R_down = [1, 0, 0;...
          0, -1, 0;...
          0, 0, -1];
R_left = [-1, 0, 0;...
           0, 0, 1;...
           0, 1, 0];
R_right = [1, 0, 0;...
           0, 0, 1;...
           0, -1, 0];

% 2nd body: left hip clamp rotor
model.joint{2} = revoluteJoint();
model.joint{2}.jointAxis = 'z';
model.Xtree{2} = plux(R_down, [0, 126, -26]*1E-3);
model.I{2} = rotorInertia;

% 3rd body: left hip clamp
model.joint{3} = revoluteJoint();
model.joint{3}.jointAxis = 'z';
model.Xtree{3} = plux(R_top, [0, 126, -87]*1E-3);
model.I{3} = mcI(inertial.mass{2},...
                 inertial.CoM{2},...
                 inertial.Ic{2});

% 4th body: hip differential rotor 1
model.joint{4} = revoluteJoint();
model.joint{4}.jointAxis = 'z';
model.Xtree{4} = plux(R_left, [0, 40, 0]*1E-3);
model.I{4} = rotorInertia;

% 5th body: hip differential rotor 2
model.joint{5} = revoluteJoint();
model.joint{5}.jointAxis = 'z';
model.Xtree{5} = plux(R_right, [0, -40, 0]*1E-3);
model.I{5} = rotorInertia;

% 6th body: gimbal link
model.joint{6} = revoluteJoint();
model.joint{6}.jointAxis = 'x';
model.Xtree{6} = plux(R_top, [0, 0, -142.5]*1E-3);
model.I{6} = mcI(inertial.mass{3},...
                 inertial.CoM{3},...
                 inertial.Ic{3});

% 7th body: thigh link
model.joint{7} = revoluteJoint();
model.joint{7}.jointAxis = 'y';
model.Xtree{7} = plux(R_top, [0,0,0]*1E-3);
model.I{7} = mcI(inertial.mass{4},...
                 inertial.CoM{4},...
                 inertial.Ic{4});

% 8th body: knee-ankle differential rotor 1
model.joint{8} = revoluteJoint();
model.joint{8}.jointAxis = 'z';
model.Xtree{8} = plux(R_right, [0, 26.55, 0]*1E-3);
model.I{8} = rotorInertia;

% 9th body: knee-ankle differential rotor 2
model.joint{9} = revoluteJoint();
model.joint{9}.jointAxis = 'z';
model.Xtree{9} = plux(R_left, [0, -26.55, 0]*1E-3);
model.I{9} = rotorInertia;

% 10th body: shank link
model.joint{10} = revoluteJoint();
model.joint{10}.jointAxis = 'y';
model.Xtree{10} = plux(R_top, [0, 0, -226.8]*1E-3);
model.I{10} = mcI(inertial.mass{5},...
                 inertial.CoM{5},...
                 inertial.Ic{5});

% 11th body: foot link
model.joint{11} = revoluteJoint();
model.joint{11}.jointAxis = 'y';
model.Xtree{11} = plux(R_top, [0, 0, -260]*1E-3);
model.I{11} = mcI(inertial.mass{6},...
                 inertial.CoM{6},...
                 inertial.Ic{6});

% 12th body: right hip clamp rotor
model.joint{12} = revoluteJoint();
model.joint{12}.jointAxis = 'z';
model.Xtree{12} = plux(R_down, [0, -126, -26]*1E-3);
model.I{12} = rotorInertia;

% 13th body: right hip clamp
model.joint{13} = revoluteJoint();
model.joint{13}.jointAxis = 'z';
model.Xtree{13} = plux(R_top, [0, -126, -87]*1E-3);
model.I{13} = mcI(inertial.mass{2},...
                 inertial.CoM{2},...
                 inertial.Ic{2});

% 14th body: hip differential rotor 1
model.joint{14} = revoluteJoint();
model.joint{14}.jointAxis = 'z';
model.Xtree{14} = plux(R_left, [0, 40, 0]*1E-3);
model.I{14} = rotorInertia;

% 15th body: hip differential rotor 2
model.joint{15} = revoluteJoint();
model.joint{15}.jointAxis = 'z';
model.Xtree{15} = plux(R_right, [0, -40, 0]*1E-3);
model.I{15} = rotorInertia;

% 16th body: gimbal link
model.joint{16} = revoluteJoint();
model.joint{16}.jointAxis = 'x';
model.Xtree{16} = plux(R_top, [0, 0, -142.5]*1E-3);
model.I{16} = mcI(inertial.mass{3},...
                 inertial.CoM{3},...
                 inertial.Ic{3});

% 17th body: thigh link
model.joint{17} = revoluteJoint();
model.joint{17}.jointAxis = 'y';
model.Xtree{17} = plux(R_top, [0,0,0]*1E-3);
model.I{17} = mcI(inertial.mass{4},...
                 inertial.CoM{4},...
                 inertial.Ic{4});

% 18th body: knee-ankle differential rotor 1
model.joint{18} = revoluteJoint();
model.joint{18}.jointAxis = 'z';
model.Xtree{18} = plux(R_right, [0, 26.55, 0]*1E-3);
model.I{18} = rotorInertia;

% 19th body: knee-ankle differential rotor 2
model.joint{19} = revoluteJoint();
model.joint{19}.jointAxis = 'z';
model.Xtree{19} = plux(R_left, [0, -26.55, 0]*1E-3);
model.I{19} = rotorInertia;

% 20th body: shank link
model.joint{20} = revoluteJoint();
model.joint{20}.jointAxis = 'y';
model.Xtree{20} = plux(R_top, [0, 0, -226.8]*1E-3);
model.I{20} = mcI(inertial.mass{5},...
                 inertial.CoM{5},...
                 inertial.Ic{5});

% 21th body: foot link
model.joint{21} = revoluteJoint();
model.joint{21}.jointAxis = 'y';
model.Xtree{21} = plux(R_top, [0, 0, -260]*1E-3);
model.I{21} = mcI(inertial.mass{6},...
                 inertial.CoM{6},...
                 inertial.Ic{6});

model = model.postProcessModel();

end