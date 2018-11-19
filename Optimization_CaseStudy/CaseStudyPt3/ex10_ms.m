% double pendulum on a cart

clear all;
close all;
clc;

N = 200; % Control discretization

opti = casadi.Opti();
import casadi.*

mcart = 0.1;
m1 = 1;
m2 = 1;
l1 = 0.3;
l2 = 0.3;
g = 9.81;

% declare variables
% states
x = opti.variable(6,N+1);
pos     = x(1,:);  % position
theta1  = x(2,:);  % angle1
theta2  = x(3,:);  % angle2
dpos    = x(4,:);
dtheta1 = x(5,:);
dtheta2 = x(6,:);

% controls
u = opti.variable(N,1); % acceleration, ddpos
% time
T = opti.variable(1,1);  % motion time

% Pendulum dynamics
ode = double_pendulum_ode(mcart,m1,m2,l1,l2,g);

% state limits
pos_lim = 2;

% input limits
u_lim = 50;

% Initial and terminal constraints
x_init = [0;pi;pi;0;0;0];  % stable
x_final = [0;0;0;0;0;0];  % swing up

T_s = SX.sym('T');

for k=1:N
   % shooting constraint
   xf = rk4(ode,T/N,x(:,k),u(k));
   opti.subject_to(x(:,k+1)==xf);
end

% path constraint
opti.subject_to(-pos_lim <= pos <= pos_lim);
opti.subject_to(-u_lim <= u <= u_lim);
opti.subject_to(T >= 0);

opti.subject_to(x(:,1)==x_init);
opti.subject_to(x(:,end)==x_final);

opti.minimize(T);

% set initial guess
opti.set_initial(theta1,linspace(pi,0,N+1));
opti.set_initial(theta2,linspace(pi,0,N+1));
opti.set_initial(T,1);

opti.callback(@(i)   plot( linspace(0,opti.debug.value(T), N+1) , opti.debug.value(x(1:3,:)) )    )

% solve optimization problem
opti.solver('ipopt');

sol = opti.solve();

% retrieve the solution
pos_opt = sol.value(pos);
dpos_opt = sol.value(dpos);
theta1_opt = sol.value(theta1);
dtheta1_opt = sol.value(dtheta1);
theta2_opt = sol.value(theta2);
dtheta2_opt = sol.value(dtheta2);
u_opt = sol.value(u);
T_opt = sol.value(T);

% time grid for printing
tgrid = linspace(0,T_opt, N+1);

figure;
subplot(3,1,1)
hold on
plot(tgrid, pos_opt, 'b')
plot(tgrid, theta1_opt, 'r')
plot(tgrid, theta2_opt, 'k')
legend('pos [m]','theta1 [rad]','theta2 [rad]')
xlabel('Time [s]')
subplot(3,1,2)
hold on
plot(tgrid, dpos_opt, 'b')
plot(tgrid, theta1_opt, 'r')
plot(tgrid, theta2_opt, 'k')
legend('dpos [m/s]','dtheta1 [rad/s]','dtheta2 [rad/s]')
xlabel('Time [s]')
subplot(3,1,3)
stairs(tgrid(1:end-1), u_opt, 'b')
legend('u [m/s^2]')
xlabel('Time [s]')

% make animation

x0 = pos_opt;
x1 = x0+l1*sin(theta1_opt);
y1 = l1*cos(theta1_opt);
x2 = x1+l2*sin(theta2_opt);
y2 = y1+l2*cos(theta2_opt);

figure(10)
axes('nextplot','replacechildren')
xlim([-0.6 0.6])
ylim([-l1-l2-0.1 l1+l2+0.1])
hold on
line([x0(1),x1(1),x2(1)],[0,y1(1),y2(1)],'LineWidth',3);
plot(x1(1),y1(1),'bo')
plot(x2(1),y2(1),'bo')
plot([-0.6,0.6],[0,0],'LineWidth',2, 'color','red')
xlim([-3,3])
ylim([-0.7,0.7])
axis equal
for k = 1:N+1
    line([x0(k),x1(k),x2(k)],[0,y1(k),y2(k)],'LineWidth',3);
    plot(x1(k),y1(k),'bo')
    plot(x2(k),y2(k),'bo')
    plot([-0.6,0.6],[0,0],'LineWidth',2, 'color','red')
    drawnow
    pause(0.1)
    cla
end

% simscape validation

open_system('double_pendulum');

double_pendulum([],[],[],'compile');

s = rand(6,1);

% In the simulink model, the order is [pos;dpos;alpha1;dalpha1;alpha2;dalpha2]
% with alpha1 = theta1 and alpha2 = theta2-theta1
r = double_pendulum(0,[s(1);s(4);s(2);s(5);s(3)-s(2);s(6)-s(5);0;0],[],'derivs');
ref=[r(1);r(3);r(3)+r(5);r(2);r(4);r(4)+r(6)];
ours = full(ode(s,0));

mismatch = norm(ref-ours)
assert(mismatch<1e-12)

double_pendulum(0,[],[],'term')

sim('double_pendulum')

