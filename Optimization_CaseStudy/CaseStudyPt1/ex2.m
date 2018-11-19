close all
clc
clear all

% Model definition
% statespace: [position;velocity]
ode = @(x) [x(2);-x(1)];
refsol = @(t) cos(t);

% Final time of integration
Tf  = 10;


% Intial condition
x0 = [1;0];

N = 50;
ts = linspace(0,Tf,N+1);
h = Tf/N;

% Euler integration
x_euler = zeros(2,N+1);
x_euler(:,1) = x0;
for k=1:N
    x_euler(:,k+1) = x_euler(:,k)+ode(x_euler(:,k))*h;
end

% Runge-kutta integration
x_rk4 = zeros(2,N+1);
x_rk4(:,1) = x0;
for k=1:N
  k1 = ode(x_rk4(:,k));
  k2 = ode(x_rk4(:,k)+h/2*k1);
  k3 = ode(x_rk4(:,k)+h/2*k2);
  k4 = ode(x_rk4(:,k)+h*k3  );
  x_rk4(:,k+1) = x_rk4(:,k) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

figure
hold on
plot(ts,refsol(ts),'k-')
plot(ts,x_euler(1,:),'rx')
plot(ts,x_rk4(1,:),'bo')
legend('ref','euler','rk4')
title('Forward simulation')
xlabel('Time [s]')

% Sensitivity to amount of control intervals

Ns = ceil(logspace(1,5));

e_euler = [];
e_rk4 = [];

for N=Ns
    h = Tf/N;
    % Euler integration
    x = x0;
    for k=1:N
        x = x+ode(x)*h;
    end
    
    e_euler = [e_euler x(1)-refsol(Tf)];

    % Runge-kutta integration
    x = x0;
    for k=1:N
      k1 = ode(x);
      k2 = ode(x+h/2*k1);
      k3 = ode(x+h/2*k2);
      k4 = ode(x+h*k3  );
      x = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end

    e_rk4 = [e_rk4 x(1)-refsol(Tf)];
        
end

figure()
loglog(Ns, abs(e_euler),'linewidth',3)
hold on
loglog(Ns, abs(e_rk4),'linewidth',3)
title(['Simulation error @ Tf = ' num2str(Tf) ' s'])
xlabel('Amount of steps')
ylabel('|| e(Tf) ||')
legend('euler','rk4')
grid on

set(gcf,'paperunits','centimeters');
set(gcf,'paperorientation','landscape');
set(gcf,'papersize',[29.7 21]/2);
set(gcf,'paperposition',[1 1 28.7 20]/2);
print('../../seminar-text/figures/integrators','-dpdf')
