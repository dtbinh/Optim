clear all;
close all;
clc;
% addpath('C:\Users\Wolf\Documents\2de Master\Optimization\optispline_windows_matlabR2014a_v0.1')
import casadi.*
import splines.*

%%
T = 1;   % End time
N = 50; % Number of control intervals
dt = T/N;
t = linspace(0,T,N+1); % time vector

meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];

opti = splines.OptiSpline();

% Set up splines
d = 10; % number of knots
L = T; % spline domain [0 , L ]
n = 3; % degree of splines
Bl  = BSplineBasis([0 , L], d, n);

% System states
p  = SX.sym('p',3,1); % object position
Rt = SX.sym('Rt',3,3); % translational Frenet-Serret frame
x = [p; Rt(:)];

% System controls (invariants)
i1 = SX.sym('i1'); % object translation speed
i2 = SX.sym('i2'); % curvature speed translational Frenet-Serret
i3 = SX.sym('i3'); % torsion speed translational Frenet-Serret
u = [i1 ; i2 ; i3];
nu = size(u,1);

% State dynamics equations of the form: dx/dt = f(x,u,t)
dRt = Rt*skew([i3;i2;0]);
dp = Rt*[i1;0;0];

rhs = [dp; dRt(:)];

% Define ordinary differential equations
ode_simp = Function('ode_simp',{x,u},{rhs});


% System controls (invariants)
u = opti.Function(Bl, [3,1]);

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

opti.subject_to(u(1)>=0);

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
opti.subject_to(Rt{1}'*Rt{1} == eye(3));

% Dynamic constraints
for k=1:N
    % Integrate current state to obtain next state
    Xk_end = rk4(ode_simp,dt,X{k},u.eval(k*dt));
    
    % Gap closing constraint
    opti.subject_to(Xk_end==X{k+1});
    
end

% Construct objective
objective_fit = 0;
for k=1:N+1
    e = p{k} - meas_pos(:,k); % position error
    objective_fit = objective_fit + e'*e;
end

objective_reg = 0;
for k=1:N-1
    e = u.eval(t(k+1)) - u.eval(t(k));
    objective_reg = objective_reg + 1e-4*e'*e;
end

% Initialize states
for k=1:N+1
    opti.set_initial(Rt{k}, eye(3));
    opti.set_initial(p{k}, meas_pos(:,k));
end

opti.minimize(objective_fit + objective_reg);
opti.solver('ipopt');

% Solve the NLP
sol = opti.solve();
%%

figure
hold on
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:),'b-')

traj = sol.value([p{:}]);
plot3(traj(1,:),traj(2,:),traj(3,:),'ro')
axis equal

view([-76 14])
figure
u_sol = sol.value(u);
plot(u_sol.list_eval(t))

u_sol_C3 = u_sol.list_eval(t);


%% Save results
save('Ex_C3spline','u_sol_C3','traj');




%% ========================================================================
%   Ex_C4Splines
%  ========================================================================

% fit = load('fit');
u_ref = u_sol;

% T = 1;   % End time
% N = 100; % Number of control intervals
% dt = T/N;
% t = linspace(0,T,N+1); % time vector

% NOTE: not used to fit here, only for comparison in figures.
% Actual fit is on spline solution from ex_C3Splines (u_ref)
meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];

opti = splines.OptiSpline();

% Set up splines
% d = 15; % number of knots
% L = T; % spline domain [0 , L ]
% n = 5; % degree of splines
% Bl  = BSplineBasis([0 , L], d, n);

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


% System controls (invariants)
u = opti.Function(Bl, [3,1]);

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

opti.subject_to(u(1)>=0); % Can only move forward

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
opti.subject_to(Rt{1}'*Rt{1} == eye(3));

P_start = [0;1;0];
P_end = [1;1;1];

opti.subject_to(p{1}==P_start);
opti.subject_to(p{end}==P_end);

% Dynamic constraints
for k=1:N
    % Integrate current state to obtain next state
    Xk_end = rk4(ode_simp,dt,X{k},u.eval(k*dt));
    
    % Gap closing constraint
    opti.subject_to(Xk_end==X{k+1});
    
end

% Construct objective
objective = 0;
for k=1:N
    e = u.eval(t(k)) - u_ref.eval(t(k)); % invariant error
    objective = objective + e'*e;
end

opti.set_initial(u, u_ref);

% Initialize states
for k=1:N+1
    opti.set_initial(Rt{k}, eye(3));
end

opti.minimize(objective);
opti.solver('ipopt');

% Solve the NLP
sol = opti.solve();

u_sol = sol.value(u);
u_sol_C4 = u_sol.list_eval(t);
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
plot(u_sol.list_eval(t))
plot(u_ref.list_eval(t))

%% Save results
save('Ex_C4spline','u_sol_C4','traj');