clear all;
close all;
clc;
% addpath('C:\Users\Wolf/casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

%% frenet serret invariants as function of s
T = 1;   % End time
N = 40; % Number of control intervals
dt = T/N;
t = linspace(0,T,N+1); % time vector

% System states
p  = SX.sym('p',3,1); % object position
Rt = SX.sym('Rt',3,3); % translational Frenet-Serret frame
x = [p; Rt(:)];

% System controls (invariants)
i1t = SX.sym('i1t'); % object translation speed
i2s = SX.sym('i2s'); % curvature speed translational Frenet-Serret
i3s = SX.sym('i3s'); % torsion speed translational Frenet-Serret
u = [i1t ; i2s ; i3s];
nu = size(u,1);

% State dynamics equations of the form: dx/dt = f(x,u,t)
dRt = i1t*Rt*skew([i3s;i2s;0]);
dp = i1t*Rt*[1;0;0];  % i1s = 1

rhs = [dp; dRt(:)];

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

opti.subject_to(U(1,:)>=0);

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
opti.subject_to(Rt{1}'*Rt{1} == eye(3));

% Dynamic constraints
for k=1:N
    % Integrate current state to obtain next state
    Xk_end = rk4(ode_simp,dt,X{k},U(:,k));
%     Xk_end = rk4(ode_simp,dt,X{k},U(:,k));

    
    % Gap closing constraint
    opti.subject_to(Xk_end==X{k+1});
    
end

% Construct objective
objective_fit = 0;
meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];
% meas_pos = [t;0.1*t.*sin(8*pi*t);0.1*t.*cos(8*pi*t)];


for k=1:N+1
    e = p{k} - meas_pos(:,k); % position error
    objective_fit = objective_fit + e'*e;
end

objective_reg = 0;
for k=1:N-1
    e = U(:,k+1) - U(:,k);
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

sol.value(objective_fit)

figure
hold on
traj = sol.value([p{:}]);
plot3(traj(1,:),traj(2,:),traj(3,:),'ko')
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
legend('invariant fit', 'reference')
axis equal
view([-76 14])

U_sol_s = sol.value(U);
save('ExC3_s','U_sol_s');


%% do same again for t

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

opti.subject_to(U(1,:)>=0);

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
opti.subject_to(Rt{1}'*Rt{1} == eye(3));

% Dynamic constraints
for k=1:N
    % Integrate current state to obtain next state
    Xk_end = rk4(ode_simp,dt,X{k},U(:,k));
%     Xk_end = rk4(ode_simp,dt,X{k},U(:,k));

    
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
    e = U(:,k+1) - U(:,k);
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
sol.value(objective_fit)
% time solution:
U_sol = sol.value(U);
% geometric solution:
u_full = full(U_sol_s);

s = zeros(1,length(t));
s(1) = 0;
for i = 2:length(t)-1
   s(i) = dt*u_full(1,i) + s(i-1);
end


% figure
% hold on
% traj = sol.value([p{:}]);
% plot3(traj(1,:),traj(2,:),traj(3,:),'ko')
% plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
% legend('spline invariants', 'measurement')
% axis equal
% view([-76 14])


figure
hold on
plot(s(1:end-1), ones(1,length(u_full(1,:))), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
plot(s(1:end-1), u_full(2,:), 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
plot(s(1:end-1), u_full(3,:), 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',2.5)
legend('i1s', 'i2s', 'i3s')
xlim([s(1) s(end-1)])

save('ExC3','U_sol');
%% plot s as function of t
% figure()
% plot(t(1:end-1),t(1:end-1).*full(U_sol_s(1,:)))

% figure()
% plot(t(1:end-1), U_sol')
% legend('i1', 'i2', 'i3')
% 
% figure()
% plot(t(1:end-1).*U_sol_s(1,:), U_sol_s')
% legend('i1_s', 'i2_s', 'i3_s')
% 
% figure()
% hold on
% plot(t(1:end-1), U_sol')
% plot(t(1:end-1), ones(1,length(U_sol_s(:,2:3)')))
% plot(t(1:end-1), (U_sol_s(:,2:3)'))
% legend('i1', 'i2', 'i3', 'i1_s', 'i2_s', 'i3_s')

