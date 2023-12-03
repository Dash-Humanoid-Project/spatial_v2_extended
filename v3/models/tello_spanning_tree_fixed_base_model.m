function model = tello_spanning_tree_fixed_base_model(gravity)

import casadi.*

if nargin == 1
    model.gravity = [0; 0; gravity];
end

model = RBD_model();

model.NB = 20;
model.parent = [0, 0, 2, 2, 2, 5, 6, 6, 6, 9, ...
                0, 0, 12, 12, 12, 15, 16, 16, 16, 19];

% import tello inertial parameters from solidworks-generated data files
inertial = fcn_import_params_solidworks('tello_params/inertial_*.txt');

% calculate rotor inertia
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

% 1st body: left hip clamp rotor
model.joint{1} = revoluteJoint();
model.joint{1}.jointAxis = 'z';
model.Xtree{1} = plux(R_down, [0, 126, -26]*1E-3);
model.I{1} = rotorInertia;

% 2nd body: left hip clamp
model.joint{2} = revoluteJoint();
model.joint{2}.jointAxis = 'z';
model.Xtree{2} = plux(R_top, [0, 126, -87]*1E-3);
model.I{2} = mcI(inertial.mass{2},...
                 inertial.CoM{2},...
                 inertial.Ic{2});

% 3rd body: hip differential rotor 1
model.joint{3} = revoluteJoint();
model.joint{3}.jointAxis = 'z';
model.Xtree{3} = plux(R_left, [0, 40, 0]*1E-3);
model.I{3} = rotorInertia;

% 4th body: hip differential rotor 2
model.joint{4} = revoluteJoint();
model.joint{4}.jointAxis = 'z';
model.Xtree{4} = plux(R_right, [0, -40, 0]*1E-3);
model.I{4} = rotorInertia;

% 5th body: gimbal link
model.joint{5} = revoluteJoint();
model.joint{5}.jointAxis = 'x';
model.Xtree{5} = plux(R_top, [0, 0, -142.5]*1E-3);
model.I{5} = mcI(inertial.mass{3},...
                 inertial.CoM{3},...
                 inertial.Ic{3});

% 6th body: thigh link
model.joint{6} = revoluteJoint();
model.joint{6}.jointAxis = 'y';
model.Xtree{6} = plux(R_top, [0,0,0]*1E-3);
model.I{6} = mcI(inertial.mass{4},...
                 inertial.CoM{4},...
                 inertial.Ic{4});

% 7th body: knee-ankle differential rotor 1
model.joint{7} = revoluteJoint();
model.joint{7}.jointAxis = 'z';
model.Xtree{7} = plux(R_right, [0, 26.55, 0]*1E-3);
model.I{7} = rotorInertia;

% 8th body: knee-ankle differential rotor 2
model.joint{8} = revoluteJoint();
model.joint{8}.jointAxis = 'z';
model.Xtree{8} = plux(R_left, [0, -26.55, 0]*1E-3);
model.I{8} = rotorInertia;

% 9th body: shank link
model.joint{9} = revoluteJoint();
model.joint{9}.jointAxis = 'y';
model.Xtree{9} = plux(R_top, [0, 0, -226.8]*1E-3);
model.I{9} = mcI(inertial.mass{5},...
                 inertial.CoM{5},...
                 inertial.Ic{5});

% 10th body: foot link
model.joint{10} = revoluteJoint();
model.joint{10}.jointAxis = 'y';
model.Xtree{10} = plux(R_top, [0, 0, -260]*1E-3);
model.I{10} = mcI(inertial.mass{6},...
                 inertial.CoM{6},...
                 inertial.Ic{6});

% 11th body: right hip clamp rotor
model.joint{11} = revoluteJoint();
model.joint{11}.jointAxis = 'z';
model.Xtree{11} = plux(R_down, [0, -126, -26]*1E-3);
model.I{11} = rotorInertia;

% 12th body: right hip clamp
model.joint{12} = revoluteJoint();
model.joint{12}.jointAxis = 'z';
model.Xtree{12} = plux(R_top, [0, -126, -87]*1E-3);
model.I{12} = mcI(inertial.mass{2},...
                 inertial.CoM{2},...
                 inertial.Ic{2});

% 13th body: hip differential rotor 1
model.joint{13} = revoluteJoint();
model.joint{13}.jointAxis = 'z';
model.Xtree{13} = plux(R_left, [0, 40, 0]*1E-3);
model.I{13} = rotorInertia;

% 14th body: hip differential rotor 2
model.joint{14} = revoluteJoint();
model.joint{14}.jointAxis = 'z';
model.Xtree{14} = plux(R_right, [0, -40, 0]*1E-3);
model.I{14} = rotorInertia;

% 15th body: gimbal link
model.joint{15} = revoluteJoint();
model.joint{15}.jointAxis = 'x';
model.Xtree{15} = plux(R_top, [0, 0, -142.5]*1E-3);
model.I{15} = mcI(inertial.mass{3},...
                 inertial.CoM{3},...
                 inertial.Ic{3});

% 16th body: thigh link
model.joint{16} = revoluteJoint();
model.joint{16}.jointAxis = 'y';
model.Xtree{16} = plux(R_top, [0,0,0]*1E-3);
model.I{16} = mcI(inertial.mass{4},...
                 inertial.CoM{4},...
                 inertial.Ic{4});

% 17th body: knee-ankle differential rotor 1
model.joint{17} = revoluteJoint();
model.joint{17}.jointAxis = 'z';
model.Xtree{17} = plux(R_right, [0, 26.55, 0]*1E-3);
model.I{17} = rotorInertia;

% 18th body: knee-ankle differential rotor 2
model.joint{18} = revoluteJoint();
model.joint{18}.jointAxis = 'z';
model.Xtree{18} = plux(R_left, [0, -26.55, 0]*1E-3);
model.I{18} = rotorInertia;

% 19th body: shank link
model.joint{19} = revoluteJoint();
model.joint{19}.jointAxis = 'y';
model.Xtree{19} = plux(R_top, [0, 0, -226.8]*1E-3);
model.I{19} = mcI(inertial.mass{5},...
                 inertial.CoM{5},...
                 inertial.Ic{5});

% 20th body: foot link
model.joint{20} = revoluteJoint();
model.joint{20}.jointAxis = 'y';
model.Xtree{20} = plux(R_top, [0, 0, -260]*1E-3);
model.I{20} = mcI(inertial.mass{6},...
                 inertial.CoM{6},...
                 inertial.Ic{6});

model = model.postProcessModel();

end