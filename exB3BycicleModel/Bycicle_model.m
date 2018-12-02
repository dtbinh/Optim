% Car A-->B multiple shooting

clear all
close all
clc
addpath('C:\Users\Wolf/casadi-windows-matlabR2016a-v3.4.5')
import casadi.*
%%
N = 20; % Control discretization
m = 80; % vehicle mass [kg]
L = 1.5;

% Create optimization environment
opti = casadi.Opti();

% declare variables
% states
x = opti.variable(N+1,1); % position
y = opti.variable(N+1,1);
theta = opti.variable(N+1,1); % angle
v = opti.variable(N+1,1); % velocity
% controls
a = opti.variable(N,1); % acceleration
delta = opti.variable(N,1);

T = opti.variable(1); % motion time

% ODE rhs function
ode = @(x,u)[x(4).*cos(x(3)); x(4).*sin(x(3)); x(4)/L*tan(u(2)) ; u(1)];  % =[xdot, ydot, thetadot, vdot]

% input constraints
F_min = -100.;
F_max = 100.;
delta_min = -pi/4;
delta_max = pi/4;

% Path constraints
v_min = -5.;
v_max = 5.;

% Initial and terminal constraints
x_init = [0., 0., pi, 0.];
x_final = [10., 10., 0., 0.];

% Construct all constraints

for k=1:N
   xk      = [x(k); y(k); theta(k); v(k)];
   xk_plus = [x(k+1); y(k+1); theta(k+1); v(k+1)];
   
   % shooting constraint
   xf = rk4(ode,T/N,xk,[a(k), delta(k)]);
   opti.subject_to(xk_plus==xf);
end

% path constraint
opti.subject_to(v_min <= v <= v_max);
opti.subject_to(F_min <= m*a <= F_max);
opti.subject_to(delta_min <= delta <= delta_max);
opti.subject_to(T >= 0);

opti.subject_to({x(1)==x_init(1), x(end)==x_final(1), y(1)==x_init(2), y(end)==x_final(2)});
opti.subject_to({theta(1)==x_init(3), theta(end)==x_final(3), v(1)==x_init(4), v(end)==x_final(4)}); 

% set initial guess
opti.set_initial(T, 10); % seconds

% Objective function
opti.minimize(T);

opti.solver('ipopt')
% solve optimization problem
sol = opti.solve();

% retrieve the solution
posx_opt = sol.value(x);
posy_opt = sol.value(y);
theta_opt = sol.value(theta);
vel_opt = sol.value(v);
a_opt = sol.value(a);
delta_opt = sol.value(delta);
T_opt = sol.value(T);

% time grid for printing
tgrid = linspace(0,T_opt, N+1);

%%
figure;
title('State trajectories')
subplot(3,1,1)
hold on
plot(tgrid, posx_opt, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',1.2)
plot(tgrid, posy_opt, 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',1.2)
xlabel('time [s]')
ylabel('position [m]')

subplot(3,1,2)
plot(tgrid, vel_opt, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',1.2)
xlabel('time [s]')
ylabel('velocity [m/s]')

subplot(3,1,3)
stairs(tgrid(1:end-1), m*a_opt, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',1.2)
xlabel('time [s]')
ylabel('acceleration [m/s^2]')

figure
stairs(tgrid(1:end-1), delta_opt, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.0)
xlabel('time [s]')
ylabel('delta')

figure
hold on
plot(posx_opt,posy_opt, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
xlabel('position-x [m]')
ylabel('position-y [m]')
title('Top view')
axis equal

disp(strcat('Optimal motion time: ' , num2str(sol.value(T)), ' s'));
