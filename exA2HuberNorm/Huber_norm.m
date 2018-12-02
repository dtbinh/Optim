clear all
close all

% addpath('C:\Users\Wolf/casadi-windows-matlabR2016a-v3.4.5')
import casadi.*

%% Data
data = [1 2-sqrt(2)/2 2 3 2;
        0 sqrt(2)/2   1 0 -1];
data_perturbed = [1 2-sqrt(2)/2 3.5 3 2;
                  0 sqrt(2)/2   1   0 -1];
% 
% x_data = data(1,:);
% y_data = data(2,:);

x_data = data_perturbed(1,:);
y_data = data_perturbed(2,:);

radius = 1;
%% Exercise A.1
% cost: e1² + e2² + e3² + e4² + e5² = ||e||² where e = [e1 e2 e3 e4 e5]^T
opti = casadi.Opti();
%% Huber norm
si = opti.variable(1,5);
p = opti.variable(1,2);
opti.minimize(sum(si));
opti.set_initial(p, [1.5;0.5])

% x_c = opti.variable(1,1);
% y_c = opti.variable(1,1);
cost = 0;
delta = 0.2; %good value for fitting
% delta = 0.5; %good value for plotting
x = SX.sym('error');
% func = piecewise(x<-delta, delta*(-x-0.5*delta), -delta<=x<=delta, 0.5*x^2, x>delta, delta*(x-0.5*delta));
func = (x<-delta)*(delta*(-x-0.5*delta)) + (-delta<=x)*(x<=delta)*(0.5*x^2) + (x>delta)*(delta*(x-0.5*delta));
huber_cost = Function('huber_cost',{x},{func});

% Plot huber function
y_value = -1:0.01:1;
x_value = -1:0.01:1;
for i = 1:length(x_value)
    y_value(i) = full(huber_cost(x_value(i)));
end

figure('Name', 'Huber function')
hold on
plot(x_value, 0.5*x_value.^2, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',2.5)
plot(x_value, delta*(abs(x_value) - 0.5*delta), 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
plot(x_value, y_value, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)


% Set up cost function
for i = 1:length(x_data)
error = sqrt((x_data(i) - p(1))^2 + (y_data(i) - p(2))^2) - 1;
cost = cost + huber_cost(error);
opti.subject_to(-si(i)<=huber_cost(error))
opti.subject_to(huber_cost(error)<=si(i))
end

opti.solver('ipopt');
sol = opti.solve();
sol.value(p)

alpha = 0:pi/100:2*pi;
p_hub = sol.value(p);
x_c_hub = p_hub(1);
y_c_hub = p_hub(2);
x_circle_hub = x_c_hub + radius * cos(alpha);
y_circle_hub = y_c_hub + radius * sin(alpha);

%% L2 norm
% coordinate of center is optimization variable.
x_c = opti.variable(1,1);
y_c = opti.variable(1,1);
% radius is 1. Could also be taken as variable to optimize.
cost = sum((sqrt((x_data - x_c).^2 + (y_data - y_c).^2) - radius).^2);

opti.minimize(cost);
opti.solver('ipopt');
sol = opti.solve();

% Plot
alpha = 0:pi/100:2*pi;
x_c_L2 = sol.value(x_c);
y_c_L2 = sol.value(y_c);
x_circle_L2 = x_c_L2 + radius * cos(alpha);
y_circle_L2 = y_c_L2 + radius * sin(alpha);

%% L1 norm

si = opti.variable(1,5);
p = opti.variable(1,2);
opti.minimize(sum(si));
opti.subject_to(-si<=sqrt((x_data - p(1)).^2 + (y_data - p(2)).^2) - 1)
opti.subject_to(sqrt((x_data - p(1)).^2 + (y_data - p(2)).^2) -1<=si)
opti.set_initial(p, [1.5;0.5])
opti.solver('ipopt');
sol = opti.solve();
sol.value(p);
sol.value(si);

alpha = 0:pi/100:2*pi;
p_L1 = sol.value(p);
x_c_L1 = p_L1(1);
y_c_L1 = p_L1(2);
x_circle_L1 = x_c_L1 + radius * cos(alpha);
y_circle_L1 = y_c_L1 + radius * sin(alpha);

%% Plot circles to compare methods
figure('Name', 'CSomparison between L1 norm, L2 norm and Huber norm')
hold on
plot(x_circle_L1, y_circle_L1, 'Color', [0.3010, 0.7450, 0.9330], 'LineWidth',2.5)
plot(x_circle_L2, y_circle_L2, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',2.5)
plot(x_circle_hub, y_circle_hub, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2.5)
scatter(x_data, y_data, 40.0, [0.8500, 0.3250, 0.0980], 'LineWidth',2.5)	
legend('L1 norm', 'L2 norm', 'Huber norm', 'data points')
xlim([0.0 4])
ylim([-1.5 1.5])

p_L1 
P_L2 = [x_c_L2, y_c_L2]
p_hub

