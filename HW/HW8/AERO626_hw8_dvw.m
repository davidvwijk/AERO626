%% AERO 626 Homework 8
%
%   Texas A&M University
%   Aerospace Engineering
%   van Wijk, David

close all; clear; clc;

plot_flag = false;

% Seed to reproduce results
% rng(100)

xaxis_sz = 22; yaxis_sz = 22; legend_sz = 20;
data = load('data_HW08.mat');

%% Part 1 - EKF Implementation

tic

Pww  = .01^2;
Pvv  = .02;
mx0  = 1.5;
Pxx0 = .15^2;
numPts = length(data.xk);

% Propagate the truth
x0 = data.x0;
x_truth = data.xk';

opts = odeset('AbsTol',1e-9,'RelTol',1e-9);

% NOTATION:
%   tkm1   = t_{k-1}        time at the (k-1)th time
%   mxkm1  = m_{x,k-1}^{+}  posterior mean at the (k-1)th time
%   Pxxkm1 = P_{xx,k-1}^{+} posterior covariance at the (k-1)th time
%   tk     = t_{k}          time at the kth time
%   mxkm   = m_{x,k}^{-}    prior mean at the kth time
%   Pxxkm  = P_{xx,k}^{-}   prior covariance at the kth time
%   mxkp   = m_{x,k}^{+}    posterior mean at the kth time
%   Pxxkp  = P_{xx,k}^{+}   posterior covariance at the kth time
%   Pzzkm  = P_{zz,k}^{-}   innovation covariance at the kth time
%   Pvvk   = P_{vv,k}       measurement noise covariance at the kth time

% declare storage space for saving state estimation error information
xcount  = 0;
xstore  = nan(1,2*numPts-1);
kxstore = nan(1,2*numPts-1);
mxstore = nan(1,2*numPts-1);
exstore = nan(1,2*numPts-1);
sxstore = nan(1,2*numPts-1);

% declare storage space
zcount  = 0;
zstore  = nan(1,numPts);
kzstore = nan(1,numPts);
mzstore = nan(1,numPts);
ezstore = nan(1,numPts);
dzstore = nan(1,numPts);
dxstore = nan(1,numPts);
szstore = nan(1,numPts);
svstore = nan(1,numPts);

% store initial data
xcount            = xcount + 1;
kxstore(:,xcount) = 0;
xstore(:,xcount)  = x0;
mxstore(:,xcount) = mx0;
sxstore(:,xcount) = sqrt(diag(Pxx0));
exstore(:,xcount) = data.x0 - mx0;

% measurement noise
Hv = 1;

% process noise
Fw = 1;

% initialize time, mean, and covariance for the EKF
mxkm1  = mx0;
Pxxkm1 = Pxx0;

% loop over the number of data points
for k = 1:numPts
    zk   = data.zk(k);     % current measurement to process
    Pvvk = Pvv;            % measurement noise covariance

    % unpack the truth -- this cannot be used in the filter, only for analysis
    xk = x_truth(k,:);       % true state

    % propagate the mean and covariance
    mxkm = recursivePropSingle(mxkm1,0);
    Fx = stateJacobianMean(mxkm1);
    Pxxkm = Fx*Pxxkm1*Fx' + Fw*Pww*Fw';

    % store a priori state information for analysis
    xcount            = xcount + 1;
    xstore(:,xcount)  = xk;
    kxstore(:,xcount) = k;
    mxstore(:,xcount) = mxkm;
    exstore(:,xcount) = xk - mxkm;
    sxstore(:,xcount) = sqrt(diag(Pxxkm));

    % compute the estimated measurement
    mzkm = measurementFunSingle(mxkm,0);

    % compute the measurement Jacobian
    Hxk = measurementJacobianMean(mxkm);

    % update the mean and covariance
    Pxzkm = Pxxkm*Hxk';
    Pzzkm = Hxk*Pxxkm*Hxk' + Hv*Pvvk*Hv';
    Kk = Pxzkm/Pzzkm;
    mxkp = mxkm + Kk*(zk - mzkm);
    Pxxkp = Pxxkm - Pxzkm*Kk' - Kk*(Pxzkm)' + Kk*(Pzzkm)*Kk';

    % store a posteriori state information for analysis
    xcount            = xcount + 1;
    xstore(:,xcount)  = xk;
    kxstore(:,xcount) = k;
    mxstore(:,xcount) = mxkp;
    exstore(:,xcount) = xk - mxkp;
    sxstore(:,xcount) = sqrt(diag(Pxxkp));

    % store measurement information for analysis
    zcount            = zcount + 1;
    zstore(:,zcount)  = zk;
    kzstore(:,zcount) = k;
    mzstore(:,zcount) = mzkm;
    ezstore(:,zcount) = zk - mzkm;
    dzstore(:,zcount) = (zk - mzkm)'*(Pzzkm\(zk - mzkm));
    dxstore(:,zcount) = (xk - mxkp)'*(Pxxkp\(xk - mxkp));
    szstore(:,zcount) = sqrt(diag(Pzzkm));
    svstore(:,zcount) = sqrt(diag(Pvvk));

    % cycle the time, mean, and covariance for the next step of the EKF
    mxkm1  = mxkp;
    Pxxkm1 = Pxxkp;
end

toc

% dxstoreEKF = dxstore;
% rmsEKF     = sqrt((sum(ezstore(1:end-1).^2))/length(ezstore(1:end-1)));
MAE_EKF  = mean(abs(exstore));
RMSE_EKF = rms(exstore);

if plot_flag
    % Plot Innovations
    plotPart1_Innovations(ezstore,szstore,svstore,xaxis_sz,yaxis_sz,legend_sz) %#ok<*UNRCH>

    % Plot Estimation Error
    plotPart1_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
end

%% Part 2 - Bootstrap Particle Filter (no resampling)

tic

% declare storage space for saving state estimation error information
xcount    = 0;
xstore    = nan(1,numPts+1);
kxstore   = nan(1,numPts+1);
mxstore   = nan(1,numPts+1);
exstore   = nan(1,numPts+1);
sxstore   = nan(1,numPts+1);
Neffstore = nan(1,numPts+1);

N = 10000;

% sample N particles from Gaussian (mx0, Pxx0)
[particles_m, weights_m]   = sampleInitialParticles(mx0,Pxx0,N);
[mxk0] = calcWMeanParticles(particles_m,weights_m,N);

% calc Neff
[Neff]  = calcNeff(weights_m);

% store initial data
xcount              = xcount + 1;
kxstore(:,xcount)   = 0;
xstore(:,xcount)    = x0;
mxstore(:,xcount)   = mxk0;
sxstore(:,xcount)   = sqrt(diag(Pxx0));
exstore(:,xcount)   = data.x0 - mxk0;
Neffstore(:,xcount) = Neff;

% process all measurements
for k = 1:numPts
    zk   = data.zk(k);     % current measurement to process

    % unpack the truth -- this cannot be used in the filter, only for analysis
    xk = x_truth(k,:);       % true state

    % propagate particles
    [particles_p] = propParticles(particles_m,Pww);

    % calculate the weights based on true measurements vs. expected
    [weights_p] = calcWeights(particles_p,weights_m,zk,Pvv);
    
    [Neff] = calcNeff(weights_p);

    [mxkp] = calcWMeanParticles(particles_p,weights_p,N);

    % posterior storing
    xcount            = xcount + 1;
    xstore(:,xcount)  = xk;
    kxstore(:,xcount) = k;
    mxstore(:,xcount) = mxkp;
    exstore(:,xcount) = xk - mxkp;
    sxstore(:,xcount) = sqrt(sum(weights_p.*(particles_p-mxkp).^2));
    Neffstore(:,xcount) = Neff;

    % initialize for next step
    particles_m = particles_p;
    weights_m   = weights_p;
    
end

toc

MAE_PF1  = mean(abs(exstore));
RMSE_PF1 = rms(exstore);

if plot_flag
    plotPart2_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
    plotPart2_Neff(Neffstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
end

%% Part 3 - Bootstrap Particle Filter (with resampling)

tic

% declare storage space for saving state estimation error information
xcount    = 0;
xstore    = nan(1,numPts+1);
kxstore   = nan(1,numPts+1);
mxstore   = nan(1,numPts+1);
exstore   = nan(1,numPts+1);
sxstore   = nan(1,numPts+1);
Neffstore = nan(1,numPts+1);

N = 10000;

% sample N particles from Gaussian (mx0, Pxx0)
[particles_m, weights_m]   = sampleInitialParticles(mx0,Pxx0,N);
[mxk0] = calcWMeanParticles(particles_m,weights_m,N);

% calc Neff
[Neff]  = calcNeff(weights_m);
Nthresh = .1*N;

% store initial data
xcount              = xcount + 1;
kxstore(:,xcount)   = 0;
xstore(:,xcount)    = x0;
mxstore(:,xcount)   = mxk0;
sxstore(:,xcount)   = sqrt(diag(Pxx0));
exstore(:,xcount)   = data.x0 - mxk0;
Neffstore(:,xcount) = Neff;
resampling_counter  = 0;

% process all measurements
for k = 1:numPts
    zk   = data.zk(k);     % current measurement to process

    % unpack the truth -- this cannot be used in the filter, only for analysis
    xk = x_truth(k,:);       % true state

    % propagate particles
    [particles_p] = propParticles(particles_m,Pww);

    % calculate the weights based on true measurements vs. expected
    [weights_p] = calcWeights(particles_p,weights_m,zk,Pvv);
    
    [Neff] = calcNeff(weights_p);
    if Neff < Nthresh
        [particles_p,weights_p] = basicResampling(particles_p,weights_p,N);
        resampling_counter = resampling_counter + 1;
    end

    [mxkp] = calcWMeanParticles(particles_p,weights_p,N);

    % posterior storing
    xcount            = xcount + 1;
    xstore(:,xcount)  = xk;
    kxstore(:,xcount) = k;
    mxstore(:,xcount) = mxkp;
    exstore(:,xcount) = xk - mxkp;
    sxstore(:,xcount) = sqrt(sum(weights_p.*(particles_p-mxkp).^2));
    Neffstore(:,xcount) = Neff;

    % initialize for next step
    particles_m = particles_p;
    weights_m   = weights_p;
    
end

toc

MAE_PF2  = mean(abs(exstore));
RMSE_PF2 = rms(exstore);

disp(' ')
disp(['Resampling was performed ', num2str(resampling_counter), ' times'])

if plot_flag
    plotPart3_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
    plotPart3_Neff(Neffstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
end

%% Part 4 - Filter Comparison

disp('EKF: RMSE  |  MAE')
disp([RMSE_EKF MAE_EKF])

disp('PF1: RMSE  |  MAE')
disp([RMSE_PF1 MAE_PF1])

disp('PF2: RMSE  |  MAE')
disp([RMSE_PF2 MAE_PF2])

%% Dynamics Functions

function [Fx] = stateJacobianMean(xkminus1)
Fx = 1 - .01*cos(xkminus1);
end

function [Hx] = measurementJacobianMean(xk)
Hx = cos(2*xk);
end

function [zk] = measurementFunSingle(xk,Pvv)
% Generate measurment given xk

vk = randn*sqrt(Pvv);
zk = .5*sin(2*xk) + vk;

end

function [xk] = recursivePropSingle(x0,Pww)
% Recursively propagate the state

xkminus1 = x0;
wkminus1 = randn*sqrt(Pww);
xk = xkminus1 - .01*sin(xkminus1) + wkminus1;

end

%% Particle Filter Functions

function [particles,weights] = sampleInitialParticles(mx0,Pxx0,N)

particles = zeros(1,N);
for i = 1:N
    x_i = mx0 + randn*sqrt(Pxx0);
    particles(:,i) = x_i;
end
weights = ones(1,N)/N;

end

function [particles_p] = propParticles(particles_m,Pww)

N = length(particles_m);
particles_p = zeros(size(particles_m));
for i = 1:N
    particles_p(:,i) = recursivePropSingle(particles_m(:,i),Pww);
end

end

function [weights_p] = calcWeights(particles_p,weights_m,zk,Pvv)
    
N = length(weights_m);
weights_p = zeros(size(weights_m));
for i = 1:N
    xk       = particles_p(:,i);
    zk_tilde = measurementFunSingle(xk,0);
    wkm1     = weights_m(:,i);
    pw       = abs(2*pi*Pvv)^(-.5)*exp(-.5*(zk - zk_tilde)'*Pvv^-1*(zk - zk_tilde));
    % pw       = normpdf(zk - zk_tilde,0,sqrt(Pvv));
    weights_p(:,i) = wkm1*pw;
end
weights_p = weights_p/sum(weights_p);

end

function [mxk] = calcWMeanParticles(particles,weights,N)

mxk = 0;
for i = 1:N
    xk = particles(:,i);
    wk = weights(:,i);
    mxk = mxk + xk*wk;
end

end

function [Neff] = calcNeff(weights)

Neff = 1/sum(weights.^2);

end

function [particles,weights] = basicResampling(particles_k,weights_k,N)

c = cumsum(weights_k);
particles = zeros(size(particles_k));
for i = 1:N
    u = rand;
    j = find(u <= c,1,'first');
    particles(:,i) = particles_k(:,j);
end
weights = ones(1,N)/N;

end

%% Plotting Functions

function plotPart1_Innovations(ezstore,szstore,svstore,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:length(ezstore);
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];
ez_scatter_opts = {100,'x','r'};

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{Innovations vs. Measurement Number using EKF}','Fontsize',25,'interpreter','latex')
a1 = scatter(measx,ezstore(1,:),ez_scatter_opts{:});
b1 = plot(measx,std_plot*szstore(1,:),'-','Color','k','LineWidth',2,'MarkerSize',20);
b2 = plot(measx,std_plot*svstore(1,:),'-','Color',[.7 .7 .7],'LineWidth',2,'MarkerSize',20);
plot(measx,-std_plot*szstore(1,:),'-','Color','k','LineWidth',2,'MarkerSize',20);
plot(measx,-std_plot*svstore(1,:),'-','Color',[.7 .7 .7],'LineWidth',2,'MarkerSize',20);
ylabel('State Innovation','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'State Innov.', [txt  ' Innov. cov.'],[txt ' Meas. noise cov.']};
legend([a1 b1 b2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')
end

function plotPart1_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:numPts;
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];

err_line_opts = {'-','LineWidth',2};
std_line_opts = {'-','LineWidth',2,'Color','k'};

measx1 = sort([measx measx]); measx1 = [0 measx1];

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{State Estimation Error vs. Measurement Number using EKF}','Fontsize',25,'interpreter','latex')
a1 = plot(measx1,exstore(1,:),err_line_opts{:});
a2 = plot(measx1,std_plot*sxstore(1,:),std_line_opts{:});
plot(measx1,-std_plot*sxstore(1,:),std_line_opts{:})
ylabel('State','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'Error',txt};
legend([a1 a2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')

end

function plotPart2_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:numPts;
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];

err_line_opts = {'-','LineWidth',2};
std_line_opts = {'-','LineWidth',2,'Color','k'};

measx1 = [0 measx];

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{State Estimation Error vs. Measurement Number using Particle Filter with No Resampling}','Fontsize',25,'interpreter','latex')
a1 = plot(measx1,exstore(1,:),err_line_opts{:});
a2 = plot(measx1,std_plot*sxstore(1,:),std_line_opts{:});
plot(measx1,-std_plot*sxstore(1,:),std_line_opts{:})
ylabel('State','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'Error',txt};
legend([a1 a2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')

end

function plotPart2_Neff(Neffstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:numPts;

err_line_opts = {'-','LineWidth',2};

measx1 = [0 measx];

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{$\hat{N}_{eff}$ vs. Measurement Number using Particle Filter with No Resampling}','Fontsize',25,'interpreter','latex')
a1 = plot(measx1,Neffstore(1,:),err_line_opts{:});
ylabel('Particles','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legend(a1,'$\hat{N}_{eff}$ ','FontSize',legend_sz,'interpreter','latex','location','southeast')

end

function plotPart3_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:numPts;
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];

err_line_opts = {'-','LineWidth',2};
std_line_opts = {'-','LineWidth',2,'Color','k'};

measx1 = [0 measx];

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{State Estimation Error vs. Measurement Number using Particle Filter with Basic Resampling}','Fontsize',25,'interpreter','latex')
a1 = plot(measx1,exstore(1,:),err_line_opts{:});
a2 = plot(measx1,std_plot*sxstore(1,:),std_line_opts{:});
plot(measx1,-std_plot*sxstore(1,:),std_line_opts{:})
ylabel('State','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'Error',txt};
legend([a1 a2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')

end

function plotPart3_Neff(Neffstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:numPts;

err_line_opts = {'-','LineWidth',2};

measx1 = [0 measx];

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{$\hat{N}_{eff}$ vs. Measurement Number using Particle Filter with Basic Resampling}','Fontsize',25,'interpreter','latex')
a1 = plot(measx1,Neffstore(1,:),err_line_opts{:});
ylabel('Particles','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legend(a1,'$\hat{N}_{eff}$ ','FontSize',legend_sz,'interpreter','latex','location','southeast')

end
