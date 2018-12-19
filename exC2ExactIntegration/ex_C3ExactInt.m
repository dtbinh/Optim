%% ========================================================================
%  Do the same with exact integration instead of RK4
%  ========================================================================
clear all
close all
%% 'x' for 'exact integration'

import casadi.*

Tx = 1;  % End time
Nx = 50; % Number of control intervals
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

%variable coefficients in the diff eq X'=h*X:
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
Hsimp = Function('Hsimp',{ux},{H});

opti = casadi.Opti();

% Create decision variables and parameters for multipleshooting
px = cell(1,Nx+1);
Rtx = cell(1,Nx+1);
Xx = cell(1,Nx+1);

for k=1:Nx+1
    % System states
    px{k} = opti.variable(3,1); % object position
    Rtx{k}  = opti.variable(3,3); % translational Frenet-Serret frame
    
    Xx{k} =  [px{k};vec(Rtx{k})];
end

Ux = opti.variable(nux,Nx);

opti.subject_to(Ux(1,:)>=0);

% FS_frame - constrain to be orthogonal (only needed for one timestep, property is propagated by integrator)
opti.subject_to(Rtx{1}'*Rtx{1} == eye(3));

% Dynamic constraints
for k=1:Nx
    % Integrate current state to obtain next state
    Xk_endx = integrate(Hsimp,dtx,Xx{k},Ux(:,k));
    
    % Gap closing constraint
    opti.subject_to(Xk_endx==Xx{k+1});
    
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
    objective_regx = objective_regx + 1e-5*(ex'*ex);
end

% Initialize states
for k=1:Nx+1
    opti.set_initial(Rtx{k}, eye(3));
    opti.set_initial(px{k}, meas_posx(:,k));
end

opti.minimize(objective_fitx + objective_regx);
opti.solver('ipopt');

% Solve the NLP
solx = opti.solve();
%%

solx.value(objective_fitx)

figure()
hold on
% plot3(meas_posx(1,:),meas_posx(2,:),meas_posx(3,:),'b-') identical
plot3(meas_posx(1,:),meas_posx(2,:),meas_posx(3,:),'g-')

trajx = solx.value([px{:}]);
plot3(trajx(1,:),trajx(2,:),trajx(3,:),'ro-')
axis equal
legend('pos_{meas}', 'exact integration')

view([-76 14])
figure
plot(solx.value(Ux)')
U_solx = solx.value(Ux);


%% Save the result
save('ExC3Exact','U_solx','trajx');

