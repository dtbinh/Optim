clear all
close all

%% load files
exact = load('ExC2Exact_int');
RK4 = load('ExC2RK4_int');


%% Plot both results and the difference

diff = exact.traj - RK4.traj;

figure('Name','Simulation results')
plot3(exact.traj(1,:),exact.traj(2,:),exact.traj(3,:),'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
hold on
plot3(RK4.traj(1,:),RK4.traj(2,:),RK4.traj(3,:),'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
legend('Exact','Runge Kutta 4')

figure('Name','Difference Exact - RK4')
plot3(diff(1,:),diff(2,:),diff(3,:), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)

figure('Name','Exact - RK4 - componentwise')
hold on
% title('Exact - RK4 - componentwise')
plot(exact.traj(1,:),'LineWidth',2.)
plot(exact.traj(2,:),'LineWidth',2.)
plot(exact.traj(3,:),'LineWidth',2.)
plot(RK4.traj(1,:),'LineWidth',2.)
plot(RK4.traj(2,:),'LineWidth',2.)
plot(RK4.traj(3,:),'LineWidth',2.)
legend('Exact x', 'Exact y', 'Exact z', 'RK4 x', 'RK4 y', 'RK4 z')


figure('Name','Difference Exact - RK4 - componentwise')
hold on
% title('Difference Exact - RK4 - componentwise')
plot(diff(1,:),'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',2.5)
plot(diff(2,:),'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
plot(diff(3,:),'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
legend('Difference x', 'Difference y', 'Difference z')


