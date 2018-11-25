clear all
close all

%% load files
exact = load('exact_int');
RK4 = load('RK4_int');


%% Plot both results and the difference

diff = exact.traj - RK4.traj;

figure('Name','Simulation results')
plot3(exact.traj(1,:),exact.traj(2,:),exact.traj(3,:),'b-')
hold on
plot3(RK4.traj(1,:),RK4.traj(2,:),RK4.traj(3,:),'r-')
legend('Exact','Runge Kutta 4')

figure('Name','Difference Exact - RK4')
plot3(diff(1,:),diff(2,:),diff(3,:),'b-')

figure('Name','Exact - RK4 - componentwise')
hold on
plot(exact.traj(1,:))
plot(exact.traj(2,:))
plot(exact.traj(3,:))
plot(RK4.traj(1,:))
plot(RK4.traj(2,:))
plot(RK4.traj(3,:))
legend('Exact x', 'Exact y', 'Exact z', 'RK4 x', 'RK4 y', 'RK4 z')


figure('Name','Difference Exact - RK4 - componentwise')
hold on
plot(diff(1,:))
plot(diff(2,:))
plot(diff(3,:))
legend('Exact x', 'Exact y', 'Exact z', 'RK4 x', 'RK4 y', 'RK4 z')


