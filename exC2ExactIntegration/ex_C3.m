clear all;
close all;
clc;

import casadi.*

T = 1;   % End time
N = 100; % Number of control intervals
dt = T/N;
t = linspace(0,T,N+1); % time vector

meas_pos = [t;0.1*t.*sin(4*pi*t);0.1*t.*cos(4*pi*t)];

% System states
p  = SX.sym('p',3,1); % object position
Rt = SX.sym('Rt',3,3); % translational Frenet-Serret frame
x = [p; Rt(:)];

% System controls (invariants)
i1 = SX.sym('i1'); % object translation speed
i2 = SX.sym('i2'); % curvature speed translational Frenet-Serret
i3 = SX.sym('i3'); % torsion speed translational Frenet-Serret
u = [i1 ; i2 ; i3];
nu = size(u,1);

% State dynamics equations of the form: dx/dt = f(x,u,t)
dRt = Rt*skew([i3;i2;0]);
dp = Rt*[i1;0;0];

rhs = [dp; dRt(:)];

% Define ordinary differential equations
ode_simp = Function('ode_simp',{x,u},{rhs});

opti = casadi.Opti();

% Create decision variables and parameters for multipleshooting
p = cell(1,N+1);
Rt = cell(1,N+1);
X = cell(1,N+1);

for k=1:N+1
    % System states
    p{k} = opti.variable(3,1); % object position
    Rt{k}  = opti.variable(3,3); % translational Frenet-Serret frame
    
    X{k} =  [p{k};vec(Rt{k})];
end

U = opti.variable(nu,N);

opti.subject_to(U(1,:)>=0);

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
opti.subject_to(Rt{1}'*Rt{1} == eye(3));

% Dynamic constraints
for k=1:N
    % Integrate current state to obtain next state
    Xk_end = rk4(ode_simp,dt,X{k},U(:,k));
    
    % Gap closing constraint
    opti.subject_to(Xk_end==X{k+1});
    
end

% Construct objective
objective_fit = 0;
for k=1:N+1
    e = p{k} - meas_pos(:,k); % position error
    objective_fit = objective_fit + e'*e;
end

objective_reg = 0;
for k=1:N-1
    e = U(:,k+1) - U(:,k);
    objective_reg = objective_reg + 1e-4*e'*e;
end

% Initialize states
for k=1:N+1
    opti.set_initial(Rt{k}, eye(3));
    opti.set_initial(p{k}, meas_pos(:,k));
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
figure
plot(sol.value(U)')

U_sol = sol.value(U);
save('fit','U_sol');

%% ========================================================================
%  Do the same with exact integration instead of RK4
%  ========================================================================

% 'x' for 'exact integration'
import casadi.*

Tx = 1;   % End time
Nx = 100; % Number of control intervals
dtx = Tx/Nx;
tx = linspace(0,Tx,Nx+1); % time vector

meas_posx = [tx;0.1*tx.*sin(4*pi*tx);0.1*tx.*cos(4*pi*tx)];

% System states
px  = SX.sym('px',3,1); % object position
Rtx = SX.sym('Rtx',3,3); % translational Frenet-Serret frame
xx = [px; Rtx(:)];

% System controls (invariants)
i1x = SX.sym('i1x'); % object translation speed
i2x = SX.sym('i2x'); % curvature speed translational Frenet-Serret
i3x = SX.sym('i3x'); % torsion speed translational Frenet-Serret
ux = [i1x ; i2x ; i3x];
nux = size(ux,1);

% State dynamics equations of the form: dx/dt = f(x,u,t)




dRtx = Rtx*skew([i3x;i2x;0]);
dpx = Rtx*[i1x;0;0];
%variable coefficients in the diff eq:
H = [0  0  0  i1x 0  0  0   0  0  0   0  0  ;
     0  0  0  0  i1x 0  0   0  0  0   0  0  ;
     0  0  0  0  0  i1x 0   0  0  0   0  0  ;
     0  0  0  0  0   0  0   0  0 -i2x 0  0  ;
     0  0  0  0  0   0  0   0  0  0 -i2x 0  ;
     0  0  0  0  0   0  0   0  0  0   0 -i2x;
     0  0  0  0  0   0  0   0  0  i3x 0  0  ;
     0  0  0  0  0   0  0   0  0  0  i3x 0  ;
     0  0  0  0  0   0  0   0  0  0   0  i3x;
     0  0  0  i2x 0  0 -i3x 0  0  0   0  0  ;
     0  0  0  0  i2x 0  0 -i3x 0  0   0  0  ;
     0  0  0  0  0  i2x 0  0 -i3x 0   0  0  ];


rhsx = [dpx; dRtx(:)]; % = H*[px; Rtx(:)]

% Define ordinary differential equations
% ode_simpx = Function('ode_simpx',{xx,ux},{rhsx});
% removed this because ode : Xdot = H*X. Need H for analytical integration.
Hsimp = Function('Hsimp',{x,u},{rhs});

optix = casadi.Opti();

% Create decision variables and parameters for multipleshooting
px = cell(1,Nx+1);
Rtx = cell(1,Nx+1);
Xx = cell(1,Nx+1);

for k=1:Nx+1
    % System states
    px{k} = optix.variable(3,1); % object position
    Rtx{k}  = optix.variable(3,3); % translational Frenet-Serret frame
    
    Xx{k} =  [px{k};vec(Rtx{k})];
end

Ux = optix.variable(nux,Nx);

optix.subject_to(Ux(1,:)>=0);

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
optix.subject_to(Rtx{1}'*Rtx{1} == eye(3));

% Dynamic constraints
for k=1:Nx
    % Integrate current state to obtain next state
    Xk_endx = integrate(dtx,Xx{k},Ux(:,k));
    
    % Gap closing constraint
    optix.subject_to(Xk_endx==Xx{k+1});
    
end

% Construct objective
objective_fitx = 0;
for k=1:Nx+1
    ex = px{k} - meas_posx(:,k); % position error
    objective_fitx = objective_fitx + ex'*ex;
end

objective_regx = 0;
for k=1:Nx-1
    ex = Ux(:,k+1) - Ux(:,k);
    objective_regx = objective_regx + 1e-4*ex'*ex;
end

% Initialize states
for k=1:Nx+1
    optix.set_initial(Rtx{k}, eye(3));
    optix.set_initial(px{k}, meas_posx(:,k));
end

optix.minimize(objective_fitx + objective_regx);
optix.solver('ipopt');

% Solve the NLP
solx = optix.solve();
%%

solx.value(objective_fitx)

figure()
hold on
% plot3(meas_posx(1,:),meas_posx(2,:),meas_posx(3,:),'b-') identical
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:),'g-')

trajx = solx.value([px{:}]);
plot3(trajx(1,:),trajx(2,:),trajx(3,:),'ro')
plot3(traj(1,:),traj(2,:),traj(3,:),'bo')
axis equal
legend('pos_{meas}', 'RK4 integration', 'exact integration')

view([-76 14])
figure
plot(solx.value(Ux)')

U_solx = solx.value(Ux);
save('fitx','U_solx');