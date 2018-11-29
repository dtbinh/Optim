clear all;
close all;
clc;

import splines.*
import casadi.*


T = 1;   % End time
N = 100; % Number of control intervals
dt = T/N;
t = linspace(0,T,N+1); % time vector

meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];


opti = splines.OptiSpline();

% Set up splines
d = 3;
L = 1;
n = 10;
Bl  = BSplineBasis([0 , L], d, n);
p   = opti.Function(Bl, [3,1]);
dp  = p.derivative(1);
Rt  = opti.Function(Bl, [3,3]);
dRt = Rt.derivative(1);

% u(1) = i1, u(2) = i2, u(3) = i3
u = opti.Function(Bl, [3,1]);

% State dynamics equations of the form: dx/dt = f(x,u,t)
opti.subject_to(dRt == 1/T*Rt*[0    0    u(2);
                               0    0   -u(3);
                              -u(2) u(3)  0  ]);
opti.subject_to(dp == 1/T*Rt*[u(1);0;0]);


opti.subject_to(u(1)>=0);
opti.subject_to(Rt.eval(0)'*Rt.eval(0) == eye(3)); % FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)


% Construct objective
% objective_fit = 0;
% e = opti.variable(3, N+1);
% for k=1:N+1
%     opti.subject_to(e(:,k) == p.eval(t(k)) - meas_pos(:,k)); % position error
%     objective_fit = objective_fit + e(:,k)'*e(:,k);
% end
% 
% objective_reg = 0;
% reg = opti.variable(3,N+1);
% for k=1:N-1
%     opti.subject_to(reg(:,k) == u.eval(t(k+1)) - u.eval(t(k)));
%     objective_reg = objective_reg + 1e-4*reg(:,k)'*reg(:,k);
% end

objective_fit = 0;
for k=1:N+1
    e= p.eval(t(k)) - meas_pos(:,k); % position error
    objective_fit = objective_fit + e'*e;
end

objective_reg = 0;
% reg = opti.variable(3,N+1);
for k=1:N-1
    reg= u.eval(t(k+1)) - u.eval(t(k));
    objective_reg = objective_reg + 1e-4*reg'*reg;
end

opti.minimize(objective_fit + objective_reg);
opti.solver('ipopt');

% Solve the NLP
sol = opti.solve();
%%

sol.value(objective_fit)

figure
hold on
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:),'b-')

traj = sol.value([p{:}]);
plot3(traj(1,:),traj(2,:),traj(3,:),'ro')
axis equal

view([-76 14])
figure('Name','Invariants solution')
plot(sol.value(U)')


U_sol = sol.value(U);
save('fit','U_sol');
