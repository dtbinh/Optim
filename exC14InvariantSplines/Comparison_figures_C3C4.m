clear
close all

%% Load data

C3ms = load('Ex_C3');
C4ms = load('Ex_C4');
C3sp = load('Ex_C3Spline');
C4sp = load('Ex_C4Spline');

meas_pos = C3ms.meas_pos;


%% Ex C3
% compare reference VS ms VS spline
figure('Name','Ex_C3: Comparison Position - reference VS ms VS spline')
hold on
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:),'b-')
plot3(C3ms.traj(1,:),C3ms.traj(2,:),C3ms.traj(3,:),'ro')
plot3(C3sp.traj(1,:),C3sp.traj(2,:),C3sp.traj(3,:),'go')
legend('Reference','Multiple shooting','Splines')

figure('Name','Ex_C3: Comparison Invariants - ms VS spline')
hold on
plot(C3ms.U_sol')
plot(C3sp.u_sol_C3)
legend('i1 - Multiple shooting',...
       'i2 - Multiple shooting',...
       'i3 - Multiple shooting',...
       'i1 - Splines',...
       'i2 - Splines',...
       'i3 - Splines')

% differences 
figure('Name','Ex_C3: Difference Invariants - ms VS spline')
plot(C3ms.U_sol' - C3sp.u_sol_C3(1:end-1,:))


%% Ex C4
% compare reference VS ms VS spline
figure('Name','Ex_C4: Comparison Position - reference VS ms VS spline')
hold on
plot3(meas_pos(1,:),meas_pos(2,:),meas_pos(3,:),'b-')
plot3(C4ms.traj(1,:),C4sp.traj(2,:),C4ms.traj(3,:),'ro')
plot3(C4sp.traj(1,:),C4sp.traj(2,:),C4sp.traj(3,:),'go')
legend('Reference','Multiple shooting','Splines')

figure('Name','Ex_C4: Comparison Invariants - ms VS spline')
hold on
plot(C4ms.U_sol')
plot(C4sp.u_sol_C4)
legend('i1 - Multiple shooting',...
       'i2 - Multiple shooting',...
       'i3 - Multiple shooting',...
       'i1 - Splines',...
       'i2 - Splines',...
       'i3 - Splines')

% differences 
figure('Name','Ex_C4: Difference Invariants - ms VS spline')
plot(C4ms.U_sol' - C4sp.u_sol_C4(1:end-1,:))