clear all;
% close all;
clc;
import casadi.*

%%
T = 1;   % End time
N = 40; % Number of control intervals
dt = T/N;
t = linspace(0,T,N+1); % time vector

fit = load('ExC3_s');
U_ref = fit.U_sol_s;

meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];
% meas_pos = [t;0.1*t.*sin(8*pi*t);0.1*t.*cos(8*pi*t)];


% System states
p = SX.sym('p',3,1); % object position
Rt  = SX.sym('Rt' ,3,3); % translational Frenet-Serret frame
x = [p;Rt(:)];

% System controls (invariants)
i1t = SX.sym('i1t'); % object translation speed
i2s = SX.sym('i2s'); % curvature speed translational Frenet-Serret
i3s = SX.sym('i3s'); % torsion speed translational Frenet-Serret
u = [i1t ; i2s ; i3s];
nu = size(u,1);

% State dynamics equations of the form: dx/dt = f(x,u,t)
dRt = i1t*Rt*skew([i3s;i2s;0]);
dp = i1t*Rt*[1;0;0];

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
% P_end = [2;2;2];


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
    e = U(:,k) - U_ref(:,k); % invariant error
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
traj = sol.value([p{:}]);
plot3(traj(1,:),traj(2,:),traj(3,:), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:),'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
plot3(P_start(1),P_start(2),P_start(3),'kx')
plot3(P_end(1),P_end(2),P_end(3),'ks')
axis equal
legend('new trajectory', 'reference')
view([-76 14])

u_new = full(sol.value(U))';
U_ref = U_ref';

s1 = zeros(1,length(t));
s1(1) = 0;
for i = 2:length(t)-1
   s1(i) = dt*u_new(i,1) + s1(i-1);
end

s2 = zeros(1,length(t));
s2(1) = 0;
for i = 2:length(t)-1
   s2(i) = dt*U_ref(i,1) + s2(i-1);
end

figure
hold on
plot(s1(1:end-1), ones(1,length(u_new(:,1))), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.0)
plot(s2(1:end-1), ones(1,length(U_ref(:,1))), 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.0)
% plot(s1(1:end-1), u_new(:,1), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.0)
% plot(s2(1:end-1), U_ref(:,1), 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.0)
plot(s1(1:end-1), u_new(:,2:3), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.0)
plot(s2(1:end-1), U_ref(:,2:3), 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.0)
legend('new trajectory', 'reference')
xlim([s1(1) s1(end-1)])


