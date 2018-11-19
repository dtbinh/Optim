% Example splines

clear all
close all
clc

import splines.*

% Make a random B-spline with 3 coefficients and plot it.
% Also plot the basis functions

% Note: feel free to play with the degree, knots, coefficients,... and see
% what happens

degree = 3
n_knots = 5
L = 1  % domain

B = BSplineBasis([0,L], n_knots, degree); %knots 0 , 0.5 , 1 ; degree 1
coeffs = randn(B.dimension(),1); %generate random coefficient values
s = splines.Function(B, coeffs); %define spline

t = [0:0.01:1]; %time vector
y = s.list_eval(t); %evaluate spline in every t
figure(1)
plot(t,y) %plot spline

b = B.list_eval(t); %evaluate basis functions in every t
b = full(b); %sparse to full matrix

figure(2)
hold all
for i = 1:size(b,2)
    plot(t,b(:,i))
end

