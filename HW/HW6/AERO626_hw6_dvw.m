%% AERO 626 Homework 6 - Problem 2
%
%   Texas A&M University
%   Aerospace Engineering
%   van Wijk, David

close all; clear; clc;
plot_flag = true;

% Seed to reproduce results
rng(10)

%% Initialization and Data Generation

Pww  = .01^2;
Pvv  = .02;
mx0  = 1.5;
Pxx0 = .15^2;
numPts = 500;

% Generate random truth
x0 = mx0 + rand*sqrt(Pxx0);

x_full = recursiveProp(x0,Pww,numPts);
z_full = measurementFun(x_full,Pvv,numPts);

%% Part A
xaxis_sz = 20; yaxis_sz = 20; legend_sz = 18;
plotPartA(x_full,z_full,xaxis_sz,yaxis_sz,legend_sz)

%% Functions

function plotPartA(x,z,xaxis_sz,yaxis_sz,legend_sz)
    measx = 1:length(x);
    figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
    title('\textbf{True Trajectory vs. Measurements: Part (a)}','Fontsize',25,'interpreter','latex')
    b1 = plot(measx,x,"Color",'b','LineWidth',2);
    a1 = scatter(measx(2:end),z(2:end),'filled','MarkerFaceColor','r');
    ylabel('State','Fontsize',yaxis_sz,'interpreter','latex')
    xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
    legendtxt = {'True Trajectory', 'Measurements'};
    legend([b1 a1],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','northwest')
end

function [z] = measurementFun(x,Pvv,numMeasurements)
% Generate measurments for k > 0 

z = nan(numMeasurements,1);
for i = 2:numMeasurements
    xk = x(i);
    vk = rand*sqrt(Pvv);
    zk = .5*sin(2*xk) + vk;
    z(i,1) = zk;
end
end

function [x] = recursiveProp(x0,Pww,numPts)
% Recursively propagate the state 

x = zeros(numPts,1);
x(1) = x0;
for i = 2:numPts
    xkminus1 = x(i-1);
    wkminus1 = rand*sqrt(Pww);
    xk = xkminus1 - .01*sin(xkminus1) + wkminus1;
    x(i) = xk;
end
end