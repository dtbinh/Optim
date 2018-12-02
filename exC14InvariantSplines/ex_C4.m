clear all;
close all;
clc;

import casadi.*

T = 1;   % End time
N = 100; % Number of control intervals
dt = T/N;
t = linspace(0,T,N+1); % time vector

fit = load('fit');
U_ref = fit.U_sol;

meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];


% System states
p = SX.sym('p',3,1); % object position
Rt  = SX.sym('Rt' ,3,3); % translational Frenet-Serret frame
x = [p;Rt(:)];

% System controls (invariants)
i1 = SX.sym('i1'); % object translation speed
i2 = SX.sym('i2'); % curvature speed translational Frenet-Serret
i3 = SX.sym('i3'); % torsion speed translational Frenet-Serret
u = [i1 ; i2 ; i3];
nu = size(u,1);

% State dynamics equations of the form: dx/dt = f(x,u,t)
dRt = Rt*skew([i3;i2;0]);
dp = Rt*[i1;0;0];

rhs = [dp;dRt(:)];

% Define ordinary differential equations
ode_simp = Function('ode_simp',{x,u},{rhs});

opti = casadi.Opti();

% Create decision variables and parameters for multipleshooting
p = cell(1,N+1);
Rt = cell(1,N+1);
X = cell(1,N+1);

for k=1:N+1
    % System states
    p{k} = opti.variable(3,1); % object position
    Rt{k}  = opti.variable(3,3); % translational Frenet-Serret frame
    
    X{k} =  [p{k};vec(Rt{k})];
end

U = opti.variable(nu,N);

opti.subject_to(U(1,:)>=0); % Can only move forward

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
opti.subject_to(Rt{1}'*Rt{1} == eye(3));

P_start = [0;1;0];
P_end = [1;1;1];

opti.subject_to(p{1}==P_start);
opti.subject_to(p{end}==P_end);

% Dynamic constraints
for k=1:N
    % Integrate current state to obtain next state
    Xk_end = rk4(ode_simp,dt,X{k},U(:,k));
    
    % Gap closing constraint
    opti.subject_to(Xk_end==X{k+1});
    
end

% Construct objective
objective = 0;
for k=1:N
    e = U(:,k) - U_ref(:,k); % position error
    objective = objective + e'*e;
end

opti.set_initial(U, U_ref);

% Initialize states
for k=1:N+1
    opti.set_initial(Rt{k}, eye(3));
end

opti.minimize(objective);
opti.solver('ipopt');

% Solve the NLP
sol = opti.solve();
%%

figure
hold on
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:),'b-')

traj = sol.value([p{:}]);
plot3(traj(1,:),traj(2,:),traj(3,:),'ro')
plot3(P_start(1),P_start(2),P_start(3),'kx')
plot3(P_end(1),P_end(2),P_end(3),'ks')
axis equal

view([-76 14])
figure
hold on
plot(sol.value(U)')
plot(U_ref')
