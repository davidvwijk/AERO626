%% AERO 626 Homework 6 - Problem 2
%
%   Texas A&M University
%   Aerospace Engineering
%   van Wijk, David

close all; clear; clc;
plot_flag = true;

% Seed to reproduce results
rng(10)

xaxis_sz = 20; yaxis_sz = 20; legend_sz = 18;

%% Part A - Initialization and Data Generation

Pww  = .01^2;
Pvv  = .02;
mx0  = 1.5;
Pxx0 = .15^2;
numPts = 500;

% Generate random truth
x0 = mx0 + rand*sqrt(Pxx0);

x_truth = recursivePropFull(x0,Pww,numPts);
z_full = measurementFunFull(x_truth,Pvv,numPts);

if plot_flag
    plotPartA(x_truth,z_full,xaxis_sz,yaxis_sz,legend_sz) %#ok<*UNRCH>
end

%% Part B - EKF Implementation

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
szstore = nan(1,numPts);
svstore = nan(1,numPts);

% store initial data
xcount            = xcount + 1;
kxstore(:,xcount) = 0;
mxstore(:,xcount) = mx0;
sxstore(:,xcount) = sqrt(diag(Pxx0));

% measurement noise
Hv = 1;

% process noise
Fw = 1;

% initialize time, mean, and covariance for the EKF
mxkm1  = mx0;
Pxxkm1 = Pxx0;

% loop over the number of data points
for k = 2:numPts
    zk   = z_full(k,:);     % current measurement to process
    Pvvk = Pvv;             % measurement noise covariance

    % unpack the truth -- this cannot be used in the filter, only for analysis
    xk = x_truth(k,:);       % true state

    % propagate the mean and covariance
    mxkm = recursivePropSingle(mxkm1,Pww);
    Fx = stateJacobianMean(mxkm1);
    Pxxkm = Fx*Pxxkm1*Fx' + Fw*Pww*Fw;

    % store a priori state information for analysis
    xcount            = xcount + 1;
    xstore(:,xcount)  = xk;
    kxstore(:,xcount) = k;
    mxstore(:,xcount) = mxkm;
    exstore(:,xcount) = xk - mxkm;
    sxstore(:,xcount) = sqrt(diag(Pxxkm));

    % compute the estimated measurement
    mzkm = measurementFunSingle(mxkm,Pvv);

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
    szstore(:,zcount) = sqrt(diag(Pzzkm));
    svstore(:,zcount) = sqrt(diag(Pvvk));

    % cycle the time, mean, and covariance for the next step of the EKF
    mxkm1  = mxkp;
    Pxxkm1 = Pxxkp;
end

if plot_flag
    % Plot Innovations
    plotPartB_Innovations(ezstore,szstore,svstore,xaxis_sz,yaxis_sz,legend_sz)

    % Plot Estimation Error
    plotPartB_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
end

%% Part C - UKF Implementation

% UKF parameters
params = {[1, 2, 0],[0.5, 2, 2]};

for i=1:length(params)
    params_curr = params(1,i);
    alpha = params_curr{1}(1);
    kappa = params_curr{1}(2);
    beta = params_curr{1}(3);

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
    szstore = nan(1,numPts);
    svstore = nan(1,numPts);

    % store initial data
    xcount            = xcount + 1;
    kxstore(:,xcount) = 0;
    mxstore(:,xcount) = mx0;
    sxstore(:,xcount) = sqrt(diag(Pxx0));

    % measurement noise
    Hv = 1;

    % process noise
    Fw = 1;

    % initialize time, mean, and covariance for the EKF
    mxkm1  = mx0;
    Pxxkm1 = Pxx0;

    % loop over the number of data points
    for k = 2:numPts
        zk   = z_full(k,:);     % current measurement to process
        Pvvk = Pvv;             % measurement noise covariance

        % unpack the truth -- this cannot be used in the filter, only for analysis
        xk = x_truth(k,:);       % true state

        % propagate the mean and covariance
        [mxkm, Pxxkm, ~] = UT(mxkm1,Pxxkm1,@recursivePropSingle,alpha,kappa,beta,Pww);
        Pxxkm = Pxxkm + Pww;

        % store a priori state information for analysis
        xcount            = xcount + 1;
        xstore(:,xcount)  = xk;
        kxstore(:,xcount) = k;
        mxstore(:,xcount) = mxkm;
        exstore(:,xcount) = xk - mxkm;
        sxstore(:,xcount) = sqrt(diag(Pxxkm));

        % update the mean and covariance
        [mzkm, Pzzkm, Pxzkm] = UT(mxkm,Pxxkm,@measurementFunSingle,alpha,kappa,beta,Pvv);
        Pzzkm = Pzzkm + Pvvk;

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
        szstore(:,zcount) = sqrt(diag(Pzzkm));
        svstore(:,zcount) = sqrt(diag(Pvvk));

        % cycle the time, mean, and covariance for the next step of the EKF
        mxkm1  = mxkp;
        Pxxkm1 = Pxxkp;
    end

    if plot_flag
        params_plot = [alpha, kappa, beta];
        % Plot innovations
        plotPartC_Innovations(ezstore,szstore,svstore,xaxis_sz,yaxis_sz,legend_sz,params_plot)

        % Plot state estimation error
        plotPartC_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz,params_plot)
    end
end
%% UKF Functions

function [m_y,Pyy,Pxy] = UT(m_x,Pxx,gfun,alpha,kappa,beta,P_noise)
n = length(m_x);
lambda = alpha^2*(n+kappa) - n;
Sxx = chol(Pxx)';
x_bar = repmat(m_x,1,2*n+1) + sqrt(n+lambda)*[zeros(n),Sxx,-Sxx];
omega_m = [lambda/(n+lambda); (1/(2*(n+lambda)))*ones(2*n,1)];
omega_c = [lambda/(n+lambda) + (1-alpha^2+beta); (1/(2*(n+lambda)))*ones(2*n,1)];

m_y = zeros(1,n);
for i = 1:(2*n+1)
    y_i = gfun(x_bar(:,i),P_noise);
    w_i_m = omega_m(i,:);
    m_y = m_y + w_i_m * y_i;
end
Pyy = zeros(n,n);
Pxy = zeros(n,n);
for i = 1:(2*n+1)
    y_i = gfun(x_bar(:,i),P_noise);
    x_i = x_bar(:,i);
    w_i_c  = omega_c(i,:);
    Pyy = Pyy + w_i_c*(y_i - m_y)*(y_i - m_y)';
    Pxy = Pxy + w_i_c*(x_i - m_x)*(y_i - m_y)';
end
end

%% Dynamics Functions

function [Fx] = stateJacobianMean(xkminus1)
Fx = 1 - .01*cos(xkminus1);
end

function [Hx] = measurementJacobianMean(xk)
Hx = cos(2*xk);
end

function [z] = measurementFunFull(x,Pvv,numMeasurements)
% Generate measurments for k > 0, for all measurements

z = nan(numMeasurements,1);
for i = 2:numMeasurements
    xk = x(i);
    vk = rand*sqrt(Pvv);
    zk = .5*sin(2*xk) + vk;
    z(i,1) = zk;
end
end

function [zk] = measurementFunSingle(xk,Pvv)
% Generate measurments for k > 0, for all measurements

vk = rand*sqrt(Pvv);
zk = .5*sin(2*xk) + vk;

end

function [x] = recursivePropFull(x0,Pww,numPts)
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

function [xk] = recursivePropSingle(x0,Pww)
% Recursively propagate the state for all measurements

xkminus1 = x0;
wkminus1 = rand*sqrt(Pww);
xk = xkminus1 - .01*sin(xkminus1) + wkminus1;

end

%% Plotting Functions

function plotPartA(x,z,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:length(x);
figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{True Trajectory vs. Measurements}','Fontsize',25,'interpreter','latex')
b1 = plot(measx,x,"Color",'b','LineWidth',2);
a1 = scatter(measx(2:end),z(2:end),'filled','MarkerFaceColor','r');
ylabel('State','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'True Trajectory', 'Measurements'};
legend([b1 a1],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','northwest')
end

function plotPartB_Innovations(ezstore,szstore,svstore,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:length(ezstore);
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{Innovations versus Measurement Number using EKF}','Fontsize',25,'interpreter','latex')
a1 = plot(measx,ezstore(1,:),'-','Color','b','LineWidth',2,'MarkerSize',20);
b1 = plot(measx,std_plot*szstore(1,:),'-','Color','k','LineWidth',2,'MarkerSize',20);
b2 = plot(measx,std_plot*svstore(1,:),'-','Color',[.7 .7 .7],'LineWidth',2,'MarkerSize',20);
plot(measx,-std_plot*szstore(1,:),'-','Color','k','LineWidth',2,'MarkerSize',20);
plot(measx,-std_plot*svstore(1,:),'-','Color',[.7 .7 .7],'LineWidth',2,'MarkerSize',20);
ylabel('State Innovation','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'State Innov.', [txt  ' Innov. cov.'],[txt ' Meas. noise cov.']};
legend([a1 b1 b2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')
end

function plotPartB_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz)
measx = 1:numPts;
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];

err_line_opts = {'-','LineWidth',1.3};
std_line_opts = {'-','LineWidth',1.3,'Color','k'};

measx1 = sort([measx measx]);

figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
title('\textbf{State Estimation Error versus Measurement Number using EKF}','Fontsize',25,'interpreter','latex')
a1 = plot(measx1(3:end),exstore(1,2:end),err_line_opts{:});
a2 = plot(measx1(3:end),std_plot*sxstore(1,2:end),std_line_opts{:});
plot(measx1(3:end),-std_plot*sxstore(1,2:end),std_line_opts{:})
ylabel('State','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'Error',txt};
legend([a1 a2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')

end

function plotPartC_Innovations(ezstore,szstore,svstore,xaxis_sz,yaxis_sz,legend_sz,params)
measx = 1:length(ezstore);
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];

a = params(1); k = params(2); b = params(3);
figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
titletxt = ['\textbf{Innovations versus Measurement Number using UFK with parameters:} '...
    '$\alpha$=',num2str(a),' $\kappa$=',num2str(k),' $\beta$=',num2str(b)];
title(titletxt,'Fontsize',25,'interpreter','latex')
a1 = plot(measx,ezstore(1,:),'-','Color','b','LineWidth',2,'MarkerSize',20);
b1 = plot(measx,std_plot*szstore(1,:),'-','Color','k','LineWidth',2,'MarkerSize',20);
b2 = plot(measx,std_plot*svstore(1,:),'-','Color',[.7 .7 .7],'LineWidth',2,'MarkerSize',20);
plot(measx,-std_plot*szstore(1,:),'-','Color','k','LineWidth',2,'MarkerSize',20);
plot(measx,-std_plot*svstore(1,:),'-','Color',[.7 .7 .7],'LineWidth',2,'MarkerSize',20);
ylabel('State Innovation','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'State Innov.', [txt  ' Innov. cov.'],[txt ' Meas. noise cov.']};
legend([a1 b1 b2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')
end

function plotPartC_EstimationError(exstore,sxstore,numPts,xaxis_sz,yaxis_sz,legend_sz,params)
measx = 1:numPts;
std_plot = 3; txt = [num2str(std_plot) '$\sigma$'];

err_line_opts = {'-','LineWidth',1.3};
std_line_opts = {'-','LineWidth',1.3,'Color','k'};

measx1 = sort([measx measx]);
a = params(1); k = params(2); b = params(3);
figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
titletxt = ['\textbf{State Estimation Error versus Measurement Number using UFK with parameters:} '...
    ' $\alpha$=',num2str(a),' $\kappa$=',num2str(k),' $\beta$=',num2str(b)];
title(titletxt,'Fontsize',25,'interpreter','latex')
a1 = plot(measx1(3:end),exstore(1,2:end),err_line_opts{:});
a2 = plot(measx1(3:end),std_plot*sxstore(1,2:end),std_line_opts{:});
plot(measx1(3:end),-std_plot*sxstore(1,2:end),std_line_opts{:})
ylabel('State','Fontsize',yaxis_sz,'interpreter','latex')
xlabel('Measurement number','Fontsize',xaxis_sz,'interpreter','latex')
legendtxt = {'Error',txt};
legend([a1 a2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')

end
