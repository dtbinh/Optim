% double pendulum with a fixed base

clear all;
close all;
clc;

N = 100; % Control discretization

opti = casadi.Opti();

% declare variables

% states
x = opti.variable(4,N+1);
theta1 = x(1,:);  % angle1
theta2 = x(2,:);  % angle2
dtheta1 = x(3,:);
dtheta2 = x(4,:);

% controls
u = opti.variable(N,1); % ddtheta1
% time
T = opti.variable(1,1);  % motion time

c = 0;  % damping term
g = 9.81;  % gravitation
m1 = 1;
m2 = 1.5;
l1 = 0.2;  % pendulum length
l2 = 0.3;  % pendulum length

% ODE rhs function
ode = @(x,u)[x(3); x(4) ; u; 1/(m2*l2)*(-m2*g*sin(x(2))-m2*l1*u*cos(x(1)-x(2))+m2*l1*x(3)^2*sin(x(1)-x(2)))];  % =xdot

% state limits

% input limits
u_min = -20;
u_max = 20;

% Initial and terminal constraints
x_init = [pi;pi;0;0];  % stable
x_final = [0;0;0;0];  % swing up

for k=1:N
   % shooting constraint
   xf = rk4(ode,T/N,x(:,k),u(k));
   opti.subject_to(x(:,k+1)==xf);
end

% path constraint
% opti.subject_to(pos_min <= pos <= pos_max);
opti.subject_to(u_min <= u <= u_max);
opti.subject_to(T >= 0);

opti.subject_to(x(:,1)==x_init);
opti.subject_to(x(:,end)==x_final);

opti.minimize(T);

% set initial guess
opti.set_initial(T,10);

opti.solver('ipopt');

opti.callback(@(i)    plot( linspace(0,opti.debug.value(T), N+1) , opti.debug.value(x(1:2,:)) )   )

sol = opti.solve();


% retrieve the solution
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
plot(tgrid, theta1_opt, 'b')
plot(tgrid, theta2_opt, 'r')
legend('theta1 [rad]','theta2 [rad]')
xlabel('Time [s]')
subplot(3,1,2)
hold on
plot(tgrid, dtheta1_opt, 'b')
plot(tgrid, dtheta2_opt, 'r')
legend('dtheta1 [rad/s]','dtheta2 [rad/s]')
xlabel('Time [s]')
subplot(3,1,3)
stairs(tgrid(1:end-1), u_opt, 'b')
legend('u [rad/s^2]')
xlabel('Time [s]')

% make animation
figure(10)
axes('nextplot','replacechildren')
xlim([-0.6 0.6])
ylim([-l1-l2-0.1 l1+l2+0.1])
line([0,0+l1*sin(theta1_opt(1))],[0,0+l1*cos(theta1_opt(1))],'LineWidth',3);
hold on
line([0+l1*sin(theta1_opt(1)),l1*sin(theta1_opt(1))+l2*sin(theta2_opt(1))],[0+l1*cos(theta1_opt(1)),l1*cos(theta1_opt(1))+l2*cos(theta2_opt(1))],'LineWidth',3);
plot(0+l1*sin(theta1_opt(1)),0+l1*cos(theta1_opt(1)),'bo')
plot(0+l1*sin(theta1_opt(1))+l2*sin(theta2_opt(1)),0+l1*cos(theta1_opt(1))+l2*cos(theta2_opt(1)),'bo')
plot([-0.6,0.6],[0,0],'LineWidth',2, 'color','red')
xlim([-0.7,0.7])
ylim([-0.7,0.7])
axis equal
for k = 1:N
    line([0,0+l1*sin(theta1_opt(k))],[0,0+l1*cos(theta1_opt(k))],'LineWidth',3);
    line([0+l1*sin(theta1_opt(k)),l1*sin(theta1_opt(k))+l2*sin(theta2_opt(k))],[0+l1*cos(theta1_opt(k)),l1*cos(theta1_opt(k))+l2*cos(theta2_opt(k))],'LineWidth',3);
    plot(0+l1*sin(theta1_opt(k)),0+l1*cos(theta1_opt(k)),'bo')
    plot(0+l1*sin(theta1_opt(k))+l2*sin(theta2_opt(k)),0+l1*cos(theta1_opt(k))+l2*cos(theta2_opt(k)),'bo')
    plot([-0.6,0.6],[0,0],'LineWidth',2, 'color','red')
    drawnow
    pause(0.1)
    cla
end
