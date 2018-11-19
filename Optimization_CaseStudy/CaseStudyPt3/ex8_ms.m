% single pendulum on a cart

clear all;
close all;
clc;

import casadi.*

N = 200; % Control discretization

opti = casadi.Opti();
% declare variables
% states 
x = opti.variable(4,N+1);
pos    = x(1,:);
theta  = x(2,:);
dpos   = x(3,:);
dtheta = x(4,:);

% controls
u =  opti.variable(N,1); % force
% time
T =  opti.variable(1,1);  % motion time

g = 9.81;  % gravitation
l = 0.2;  % pendulum length
m = 1; % pendulum mass
mcart = 0.5;

% ODE rhs function
ode = @(x,u)[x(3); x(4) ;((u+m*l*x(4)*x(4)*sin(x(2))-m*g*sin(x(2))*cos(x(2)))/(mcart+m-m*cos(x(2))*cos(x(2))));  g/l*sin(x(2))-cos(x(2))*((u+m*l*x(4)*x(4)*sin(x(2))-m*g*sin(x(2))*cos(x(2)))/(l*(mcart+m-m*cos(x(2))*cos(x(2)))))];  % =xdot
% represents:
% ddpos = ((u+m*l*x(4)*x(4)*sin(x(2))-m*g*sin(x(2))*cos(x(2)))/(mcart+m-m*cos(x(2))*cos(x(2))));
% ode = @(x,u)[x(3); x(4) ;ddpos;  g/l*sin(x(2))-cos(x(2))*ddpos];  % =xdot
% but required to write in full, because x and u are only defined locally

% state limits
pos_lim = 3;
dpos_lim = 2;

% input limits
u_lim = 10;

% Initial and terminal constraints
x_init = [0;pi;0;0];  % stable
x_final = [0;0;0;0];  % swing up

for k=1:N
   % shooting constraint
   xf = rk4(ode,T/N,x(:,k),u(k));
   opti.subject_to(x(:,k+1)==xf);
end

% path constraint
opti.subject_to(-pos_lim <= pos <= pos_lim);
opti.subject_to( -u_lim <= u <= u_lim);
opti.subject_to(T >= 0);

opti.subject_to(x(:,1)==x_init);
opti.subject_to(x(:,end)==x_final);

% set initial guess
opti.set_initial(theta,linspace(pi,0,N+1));
opti.set_initial(T,1);

opti.minimize(T); %sum(u.^2);

% solve optimization problem
opti.solver('ipopt')

sol = opti.solve();
% retrieve the solution

pos_opt = sol.value(pos);
theta_opt = sol.value(theta);
dpos_opt = sol.value(dpos);
dtheta_opt = sol.value(dtheta);

u_opt = sol.value(u);
T_opt = sol.value(T);

% time grid for printing
tgrid = linspace(0,T_opt, N+1);

figure;
subplot(3,1,1)
hold on
plot(tgrid, theta_opt, 'b')
plot(tgrid, pos_opt, 'b')
legend('theta [rad]','pos [m]')
xlabel('Time [s]')
subplot(3,1,2)
hold on
plot(tgrid, dtheta_opt, 'b')
plot(tgrid, dpos_opt, 'b')
legend('dtheta [rad/s]','dpos [m/s]')
xlabel('Time [s]')
subplot(3,1,3)
stairs(tgrid(1:end-1), u_opt, 'b')
legend('u [m/s^2]')
xlabel('Time [s]')

% make animation
figure(10)
axes('nextplot','replacechildren')
xlim([-0.6 0.6])
ylim([-l-0.1 l+0.1])
line([x_init(1),x_init(1)+l*sin(theta_opt(k))],[0,0+l*cos(theta_opt(1))],'LineWidth',3)
hold on
plot(pos_opt(k)+l*sin(theta_opt(k)),0+l*cos(theta_opt(k)),'bo')
plot([-0.6,0.6],[0,0],'LineWidth',2, 'color','red')
xlim([-0.7,0.7])
ylim([-0.7,0.7])
axis equal
for k = 1:N
    line([pos_opt(k),pos_opt(k)+l*sin(theta_opt(k))],[0,0+l*cos(theta_opt(k))],'LineWidth',3)
    plot(pos_opt(k)+l*sin(theta_opt(k)),0+l*cos(theta_opt(k)),'bo')
    plot([-0.6,0.6],[0,0],'LineWidth',2, 'color','red')
    drawnow
    pause(0.1)
    cla
end
