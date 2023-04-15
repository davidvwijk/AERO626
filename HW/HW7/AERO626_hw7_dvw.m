%% AERO 626 Homework 7
%
%   Texas A&M University
%   Aerospace Engineering
%   van Wijk, David

close all; clear; clc;
plot_flag = true;

% Seed to reproduce results
% rng(10)

xaxis_sz = 20; yaxis_sz = 20; legend_sz = 18;

%% Part A - Initialization and Data Processing

% Load in the provided data
data = load('data_HW07.mat'); 
x0 = data.x0;
xk = data.xk(:,1:20); 
zk = data.zk(:,1:20);

% Define the dynamics and measurement Jacobians (constant)
Fx = [1 1; 0 1];
Hx = [1 0];

% Process noise and measurement noise covariances
Pww = diag([.01,.01]);
Pvv = 1;

% Define weights, means and covariances for the components of the GM
w1 = .2; m1x0 = [-2.5; -.5]; P1xx0 = diag([.2,.1]);
w2 = .3; m2x0 = [-1;    .2]; P2xx0 = diag([.25,.05]);
w3 = .1; m3x0 = [-.5;  -.3]; P3xx0 = diag([.2,.1]);
w4 = .4; m4x0 = [1.2;    1]; P4xx0 = diag([1,.3]);

% Sample an initial guess first by selecting the component of the GM and
% then sampling that Gaussian
sampleGaussian = rand;
if sampleGaussian < w1
    mx0 = m1x0 + diag(sqrt(P1xx0)).*randn(2,1);
elseif sampleGaussian < w1 + w2
    mx0 = m2x0 + diag(sqrt(P2xx0)).*randn(2,1);
elseif sampleGaussian < w1 + w2 + w3
    mx0 = m3x0 + diag(sqrt(P3xx0)).*randn(2,1);
else
    mx0 = m4x0 + diag(sqrt(P4xx0)).*randn(2,1);
end


for i = 1:length(zk)






end

