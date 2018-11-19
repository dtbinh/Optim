% Car A-->B single shooting

clear all
close all
clc

N = 100;  % control discretization
m = 500;  % vehicle mass [kg]

% Create optimization environment
opti = casadi.Opti();

% declare variables
x1 = opti.variable(4,1); % initial state
u = opti.variable(2,N); % control
T = opti.variable(); % motion time

% ODE right hand side
ode = @(x,u)[x(3); x(4); u(1); u(2)];

% input constraints
F_min = -2500.;
F_max = 2500.;

% path constraints
v_min = -10.;
v_max = 10.;

% initial and terminal constraints
x_init = [0.; 0.; 0.; 0.];
x_final = [10.; 10.; 0.; 0.];

opti.subject_to(x1==x_init);
x = x1;
for k=1:N   
   % shooting constraint
   h = T/N;
   k1 = ode(x,       u(:,k));    %= u(:,k);
   k2 = ode(x+h/2*k1,u(:,k));
   k3 = ode(x+h/2*k2,u(:,k));
   k4 = ode(x+h*k3,  u(:,k));
   x = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% path constraint
opti.subject_to(v_min<= x(3,:) <= v_max);
opti.subject_to(v_min<= x(4,:) <= v_max);
opti.subject_to(F_min <= m*u(1,:) <= F_max);
opti.subject_to(F_min <= m*u(2,:) <= F_max);
opti.subject_to(T >= 0);

opti.subject_to(x==x_final);

% set initial guess
opti.set_initial(T, 5); % seconds

% Objective function
opti.minimize(T);

opti.solver('ipopt');
% solve optimization problem
sol = opti.solve();
