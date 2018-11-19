function ode = double_pendulum_ode( mcart, m1, m2, l1, l2, g )
% 
% Returns an ode-right hand side evaluating function:
%   ode(x,u)
%
% x -- states: pos theta1 theta2 dpod dtheta1 dtheta2
% u -- control: 
%
% Example:
% ode = double_pendulum_ode(mcart,m1,m2,l1,l2,g);
% ode([0;pi;pi;0;0;0],0)

import casadi.*

% symbolic description of Lagrangian
pos_s = SX.sym('pos');
theta1_s = SX.sym('theta1');
theta2_s = SX.sym('theta2');
dpos_s = SX.sym('dpos');
dtheta1_s = SX.sym('dtheta1');
dtheta2_s = SX.sym('dtheta2');

q = [pos_s;theta1_s;theta2_s];
dq = [dpos_s;dtheta1_s;dtheta2_s];
ddq = SX.sym('ddq',size(q));

x0 = pos_s;
y0 = 0;
x1 = x0+l1*sin(theta1_s);
y1 = y0+l1*cos(theta1_s);
x2 = x1+l2*sin(theta2_s);
y2 = y1+l2*cos(theta2_s);

E_pot = m1*g*y1+m2*g*y2;

v0 = jtimes([x0;y0],q,dq);
v1 = jtimes([x1;y1],q,dq);
v2 = jtimes([x2;y2],q,dq);

E_kin = 1/2*mcart*(v0'*v0)+1/2*m1*(v1'*v1) + 1/2*m2*(v2'*v2);

% Lagrange function for translation part
Lag = E_kin - E_pot;

u_s = SX.sym('u');

% Lagrange equations: eq(q,dq,ddq,u) == 0
eq = jtimes(gradient(Lag,dq),[q;dq],[dq;ddq]) - gradient(Lag,q) - [u_s;0;0];
% jtimes(f,x,y) = jacobian(f,x)*y

% Write ddq as function of q,dq,u
ddqsol = -jacobian(eq,ddq)\substitute(eq,ddq,0);

% ODE rhs function
ode = casadi.Function('ode',{[q;dq],u_s},{[dq;ddqsol]},{'x','u'},{'xdot'});

end

