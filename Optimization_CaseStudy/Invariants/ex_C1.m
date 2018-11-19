clear all;
close all;
clc;

T = 1;   % End time
N = 100; % Number of control intervals
dt = T/N;
t = linspace(0,T,N+1); % time vector

meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];

% states: position;velocity;acceleration
% controls: jerk
ode = @(x,u) [x(4:end);u];

opti = casadi.Opti();

% Create decision variables and parameters for multipleshooting
p = cell(1,N+1);
v = cell(1,N+1);
a = cell(1,N+1);
X = cell(1,N+1);

for k=1:N+1
    % System states
    p{k} = opti.variable(3,1); % object position
    v{k} = opti.variable(3,1); % object velocity
    a{k} = opti.variable(3,1); % object acceleration
    
    X{k} =  [p{k};v{k};a{k}];
end

U = opti.variable(3,N);

P_start = [0;0;0];
P_end = [1;1;1];
%P_end = meas_pos(:,end);

opti.subject_to(p{1}==P_start);
opti.subject_to(p{end}==P_end);

% Dynamic constraints
for k=1:N
    % Integrate current state to obtain next state
    Xk_end = rk4(ode,dt,X{k},U(:,k));
    
    % Gap closing constraint
    opti.subject_to(Xk_end==X{k+1});
    
end

% Construct objective
objective = 0;
for k=1:N+1
    e = p{k} - meas_pos(:,k); % position error
    objective = objective + e'*e;
end


% Construct objective
for k=1:N
    e = U(:,k); % jerk
    objective = objective + 1e-7*e'*e;
end

% Initialize states
for k=1:N+1
    opti.set_initial(p{k}, meas_pos(:,k));
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
