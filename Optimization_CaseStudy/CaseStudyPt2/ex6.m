% Car A-->B splines

close all; clear all; clc;

import splines.*

%This example calculates the time-optimal motion trajectory for moving a
%tool from A to B. Only the kinematic constraints of the machine are
%considered, i.e. limits on speed and acceleration.

%Get user input
xgoal = 10; %target position of the tool
ygoal = 10;

%Define constants
%Note: it is also possible to impose separate speed limits for x and y-axis
vmin = -10;%speed limits [m/s]
vmax =  10;
m = 500; % vehicle mass [kg]
Fmin = -2500; %minimum force[N]
Fmax =  2500; %maximum force [N] 

opti = OptiSpline();

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
x0_guess = linspace(x0,xgoal,Bl.dimension);
x1_guess = linspace(y0,ygoal,Bl.dimension);
%try a straight line between start and end as init_guess
x0_guess = Function.linear(t0, x0_guess);
x1_guess = Function.linear(t0, x1_guess);
%or easier:
x0_guess = Function.linear([0,1], [x0,xgoal]);
x1_guess = Function.linear([0,1], [y0,ygoal]);
Tvar = sqrt((xgoal-x0)^2+(ygoal-y0)^2)/(vmax/2); %suppose you move at half speed on average
%Note: this is not a good guess in each case

%Initial plot
figure(1)               %create figure
hold on
plot(x0,y0,'ro')        %mark start
plot(xgoal,ygoal,'ro')  %mark goal

%Set up optimization problem
opti.minimize(T);  %objective function: motion time
    
%Define constraints: initial and final conditions, kinematic constraints motion time constraint
opti.subject_to({x(1).eval(0) == x0,    x(2).eval(0) == y0});  %initial positions
opti.subject_to({x(2).eval(1) == xgoal, x(1).eval(1) == ygoal}); %goal position
opti.subject_to({vx.eval(0)   == vx0,   vy.eval(0)   == vy0}); %initial speed
opti.subject_to({vy.eval(1)   == 0,     vx.eval(1)   == 0}); %final speed
opti.subject_to(T   * vmin <= vx <= T   * vmax); %velocity constraints
opti.subject_to(T   * vmin <= vy <= T   * vmax);
opti.subject_to(T^2 * Fmin <= m*ax <= T^2 * Fmax); %acceleration constraints
opti.subject_to(T^2 * Fmin <= m*ay <= T^2 * Fmax);
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

opti.solver('ipopt',struct('ipopt',struct('tol',1e-12)));

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
axis equal
figure(2)
subplot(211)
plot(s * T, vx.list_eval(s) / T) %  / T to convert to real speed [m/s]
hold on
plot(s * T, vy.list_eval(s) / T, 'b')
plot(s * T, vmax*ones(size(s,2)),'r-')
plot(s * T, vmin*ones(size(s,2)),'r-')
subplot(212)
plot(s * T, m*ax.list_eval(s) / T^2)
hold on
plot(s * T, m*ay.list_eval(s) / T^2, 'b')
plot(s * T, Fmax*ones(size(s,2)),'r-')
plot(s * T, Fmin*ones(size(s,2)),'r-')
