clear all;
close all;
clc;

import casadi.*

T = 1;   % End time
N = 100; % Number of control intervals
dt = T/N;

% System states
p = SX.sym('p',3,1); % object position
Rt  = SX.sym('Rt',3,3); % translational Frenet-Serret frame
x = [p;Rt(:)];

% System controls (invariants)
i1 = SX.sym('i1'); % object translation speed
i2 = SX.sym('i2'); % curvature speed translational Frenet-Serret
i3 = SX.sym('i3'); % torsion speed translational Frenet-Serret
u = [i1 ; i2 ; i3];

% State dynamics equations of the form: dx/dt = f(x,u,t)
dp = Rt*[i1;0;0];
dRt = Rt*skew([i3;i2;0]);

rhs = [dp; dRt(:)];

% Define ordinary differential equations
ode_simp = Function('ode_simp',{x,u},{rhs});


R0 = eye(3);
X = [0;0;0;R0(:)];
traj=[X(1:3)];
for k=1:N
    X = full(rk4(ode_simp,dt,X,[2;-10;-8]));
    traj=[traj X(1:3)];
end

figure
hold on
plot3(traj(1,:),traj(2,:),traj(3,:),'b-')
axis equal