close all; clear all; clc;

import splines.*
import casadi.*
Opti;

%This example calculates the time-optimal motion trajectory for moving a
%vehicle from A to B.

%Get user input
start = [0,0];
goal = [10,10]; 

% obstacle info
obs_pos = [5,5];  % [5,5.05]
obs_rad = 1;
obs_pos1 = [3,2];
obs_rad1 = 1;
obs_pos2 = [6,7];
obs_rad2 = 1;

%Define constants
%Note: it is also possible to impose separate speed limits for x and y-axis
vmin = -10;  % speed limits [m/s]
vmax =  10;
m = 500;  % vehicle mass [kg]
Fmin = -2500;  % minimum force[N]
Fmax =  2500;  % maximum force [N] 

opti = splines.OptiSpline();

%Set up splines
d   = 3;      % spline degree
L   = 1;      % range [0 , L]
n   = 11;     % number of knots
Bl  = BSplineBasis([0 , L], d, n);  %make B-spline basis
x   = opti.Function(Bl, [1,2]);   %variable: movement spline. Shape/Size: 1 by 2 (= [x , y])
vx  = x(1).derivative(1);           %velocity spline
vy  = x(2).derivative(1);
ax  = x(1).derivative(2);           %acceleration spline
ay  = x(2).derivative(2);
T   = opti.variable(1, 1); %variable: movement time

%In this problem there are two variables present.
%On the one hand there is the motion time T. Since we want to move the tool
%as fast as possible from A to B, you try to minimize this one
%On the other hand there is the movement spline, this one represents the
%motion trajectory as a function of time. The example considers a movement
%in the 2D plane, so x is of size 1 by 2. x represents the motion trajectory
%that brings us from A to B in minimal time.

%Initial system state
x0   = 0;   %initial tool position
y0   = 0;
vx0  = 0;   %initial tool speed
vy0  = 0;

%Initial guess
t0 = linspace(0,1,Bl.dimension);
x0_guess = linspace(x0,goal(1),Bl.dimension);
x1_guess = linspace(y0,goal(2),Bl.dimension);
%try a straight line between start and end as init_guess
x0_guess = Function.linear(t0, x0_guess);
x1_guess = Function.linear(t0, x1_guess);
%or easier:
x0_guess = Function.linear([0,1], [x0,goal(1)]);
x1_guess = Function.linear([0,1], [y0,goal(2)]);
Tvar = sqrt((goal(1)-x0)^2+(goal(2)-y0)^2)/(vmax/2); %suppose you move at half speed on average
%Note: this is not a good guess in each case

%Initial plot
figure(1)               %create figure
plot(x0,y0,'ro')        %mark start
hold on
plot(goal(1),goal(2),'ro')  %mark goal

%Set up optimization problem
opti.minimize(T);  %objective function: motion time
    
%Define constraints: initial and final conditions, kinematic constraints motion time constraint

opti.subject_to({x(1).eval(0) == x0,    x(2).eval(0) == y0});  %initial positions
opti.subject_to({x(2).eval(1) == goal(1), x(1).eval(1) == goal(2)}); %goal position
opti.subject_to({vx.eval(0)   == vx0,   vy.eval(0)   == vy0}); %initial speed
opti.subject_to({vy.eval(1)   == 0,     vx.eval(1)   == 0}); %final speed
opti.subject_to(T   * vmin <= vx <= T   * vmax); %velocity constraints
opti.subject_to(T   * vmin <= vy <= T   * vmax);
opti.subject_to(T^2   * Fmin <= m*ax <= T^2   * Fmax); %acceleration constraints
opti.subject_to(T^2   * Fmin <= m*ay <= T^2   * Fmax);
opti.subject_to((x(1)-obs_pos(1))^2 + (x(2)-obs_pos(2))^2 >= obs_rad^2);
opti.subject_to(T >= 0); %motion time, is always positive
%Note: the final speed is set as zero, we want to stop at the goal position

%Note: the T and T^2 in the velocity and acceleration constraints show up
%because the problem description is nondimensional. We describe everything
%on a basis of [0,1], meaning that everything is expressed in a dimensionless
%time 'tau' which is [0,1]. This means that x is a function of tau, or: x(tau)
%So the real time 't' = tau * T. 
%When we want to find the velocity this becomes dx(tau)/dt = dx(tau)/dtau * dtau/dt 
%And since tau = t / T, dtau/dt = 1/T 
%So: v(t) = 1/T * dx/dtau. Then we set v(t) = 1/T * dx/dtau <= vmax
%After repeating this procedure the T^2 shows up in the acceleration constraint.
opti.solver('ipopt',struct('ipopt',struct('max_iter',2000)));

opti.set_initial(T, Tvar); % seconds
opti.set_initial(x, [x0_guess,x1_guess]);

%Solve optimization problem
sol = opti.solve();

%The actual variables of this problem are the motion time T and the
%coefficients of the spline x.

%Process solution
T = sol.value(T);
x = sol.value(x);
vx = x(1).derivative(1);  %obtain velocity
vy = x(2).derivative(1);
ax = x(1).derivative(2);  %obtain velocity
ay = x(2).derivative(2);

%Plot solutions
s = 0:0.01:1; %create time vector: go from 0 to 1 in steps of 0.01s
figure(1)     %select right figure
vehicle_path = plot(x(1).list_eval(s),x(2).list_eval(s), 'b');%plot calculated path
hold on
arc = 0:0.01:2*pi;
circle = plot(obs_pos(1)+obs_rad*cos(arc), obs_pos(2)+obs_rad*sin(arc),'r');
% circle = plot(obs_pos1(1)+obs_rad1*cos(arc), obs_pos1(2)+obs_rad1*sin(arc),'r');
% circle = plot(obs_pos2(1)+obs_rad2*cos(arc), obs_pos2(2)+obs_rad2*sin(arc),'r');
figure(2)
subplot(211)
plot(s * T, vx.list_eval(s) / T) %  / T to convert to real speed [m/s]
hold on
plot(s * T, vmax*ones(size(s,2)),'r-')
plot(s * T, vmin*ones(size(s,2)),'r-')
subplot(212)
plot(s * T, vy.list_eval(s) / T, 'b')
hold on
plot(s * T, vmax*ones(size(s,2)),'r-')
plot(s * T, vmin*ones(size(s,2)),'r-')
figure(3)
subplot(311)
plot(s * T, m*ax.list_eval(s) / T^2)
hold on
plot(s * T, Fmax*ones(size(s,2)),'r-')
plot(s * T, Fmin*ones(size(s,2)),'r-')
subplot(312)
plot(s * T, m*ay.list_eval(s) / T^2, 'b')
hold on
plot(s * T, Fmax*ones(size(s,2)),'r-')
plot(s * T, Fmin*ones(size(s,2)),'r-')
