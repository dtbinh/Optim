close all
clear all

%% load results

exact = load('ExC3Exact');
RK4 = load('ExC3');


%% Plot both results and the difference
diff = exact.U_solx - RK4.U_sol;

figure('Name','Simulation results Exact - RK4')
plot3(exact.trajx(1,:),exact.trajx(2,:),exact.trajx(3,:),'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
hold on
plot3(RK4.traj(1,:),RK4.traj(2,:),RK4.traj(3,:),'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
legend('Exact','Runge Kutta 4')

figure('Name','Difference Exact - RK4')
plot3(diff(1,:),diff(2,:),diff(3,:),'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)

figure('Name','Exact - RK4 - componentwise')
hold on
plot(exact.trajx(1,:),'LineWidth',2.)
plot(exact.trajx(2,:),'LineWidth',2.)
plot(exact.trajx(3,:),'LineWidth',2.)
plot(RK4.traj(1,:),'LineWidth',2.)
plot(RK4.traj(2,:),'LineWidth',2.)
plot(RK4.traj(3,:),'LineWidth',2.)
legend('Exact x', 'Exact y', 'Exact z', 'RK4 x', 'RK4 y', 'RK4 z')


figure('Name','Difference Exact - RK4 - Componentwise')
hold on
plot(diff(1,:)./exact.U_solx(1,:),'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',2.5)
plot(diff(2,:)./exact.U_solx(2,:),'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
plot(diff(3,:)./exact.U_solx(3,:),'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
legend({'Difference i_1', 'Difference i_2', 'Difference i_3'},'FontSize',14)
xlabel('Time [s]')
