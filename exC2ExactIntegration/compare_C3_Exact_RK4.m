close all
clear all

%% load results

exact = load('fitx');
RK4 = load('fit');


%% Plot both results and the difference
diff = exact.trajx - RK4.traj;

figure('Name','Simulation results Exact - RK4')
plot3(exact.trajx(1,:),exact.trajx(2,:),exact.trajx(3,:),'b-')
hold on
plot3(RK4.traj(1,:),RK4.traj(2,:),RK4.traj(3,:),'r-')
legend('Exact','Runge Kutta 4')

figure('Name','Difference Exact - RK4')
plot3(diff(1,:),diff(2,:),diff(3,:),'b-')

figure('Name','Exact - RK4 - componentwise')
hold on
plot(exact.trajx(1,:))
plot(exact.trajx(2,:))
plot(exact.trajx(3,:))
plot(RK4.traj(1,:))
plot(RK4.traj(2,:))
plot(RK4.traj(3,:))
legend('Exact x', 'Exact y', 'Exact z', 'RK4 x', 'RK4 y', 'RK4 z')


figure('Name','Difference Exact - RK4 - Componentwise')
hold on
plot(diff(1,:))
plot(diff(2,:))
plot(diff(3,:))
legend('Difference x', 'Difference y', 'Difference z')

