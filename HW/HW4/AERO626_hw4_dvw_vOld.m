%% AERO 626 Homework 4
%
%   Texas A&M University
%   Aerospace Engineering
%   van Wijk, David

%% Problem 2

plot_flag = true;

rng(100) % Seed to reproduce results
opts = odeset('AbsTol',1e-6,'RelTol',1e-6);

mx0 = [1; 0];        %[m], [m/s]
Pxx0 = diag([2,1]);  %[m]^2, [m/s]^2
Pvv = .1^2;          %[m]^2 constant covariance

N = 21;
tspan = linspace(0,20,N);
Fx = [0 1; -1 0];
Hx = [1 0]; 
Hv = 1;
x0 = mx0 + chol(Pxx0)'*randn(2,1);

% Store data
std_x_full = zeros(N*2-1,2);
std_innov_full = zeros(N-1,2);

% Initialize for first step
xkm1 = x0;
mxkm1 = mx0;
Pxxkm1 = Pxx0;
std_x_full(1,:) = sqrt(diag(Pxx0))';
count = 1;
for i = 2:N
    % Propagate the truth state
    [~,X] = ode45(@(t,x) prop(t,x,Fx),[tspan(i-1),tspan(i)],xkm1,opts);
    xk = X(end,1:2)';
    
    % Generate measurement
    zk = Hx*xk + chol(Pvv)'*randn;

    % Propagate step of Kalman filter
    [~,X] = ode45(@(t,x) prop(t,x,Fx),[tspan(i-1),tspan(i)],[mxkm1;Pxxkm1(:)],opts);
    mxkm = X(end,1:2)';
    Pxxkm = reshape(X(end,3:end)',2,2);
    
    % Store stuff
    count = count + 1;
    std_x_full(count,:) = sqrt(diag(Pxxkm))';

    mzkm = Hx*mxkm;
    Pxzkm = Pxxkm*Hx';
    Pzzkm = Hx*Pxxkm*Hx' + Hv*Pvv*Hv';

    % Kalman Gain
    Kk = Pxzkm/Pzzkm;
    
    % Update mean and covariance
    mxkp = mxkm + Kk*(zk - mzkm);
    Pxxkp = Pxxkm - Pxzkm*Kk' - Kk*(Pxzkm)' + Kk*(Pzzkm)*Kk';
    
    % Re-initialize for next loop
    xkm1 = xk;
    mxkm1 = mxkp;
    Pxxkm1 = Pxxkp;

    % Store stuff
    count = count + 1;
    std_x_full(count,:) = sqrt(diag(Pxxkp))';
    std_innov_full(i,:) = sqrt(diag(Pzzkm))';
end

if plot_flag 
    xaxis_sz = 20; yaxis_sz = 20; legend_sz = 18;
    
    txt = '1$\sigma$';
    figure; grid on; set(gcf, 'WindowState', 'maximized');
    subplot(2,1,1); hold on;
    title('\textbf{State Covariance vs. Time}','Fontsize',25,'interpreter','latex')
    t_mod = sort([tspan tspan]);
    a1 = plot(t_mod(2:end),std_x_full(:,1),'-','Color','k','LineWidth',2,'MarkerSize',20);
    plot(t_mod(2:end),-std_x_full(:,1),'-','Color','k','LineWidth',2,'MarkerSize',20)
    ylabel('x [$m$]','Fontsize',yaxis_sz,'interpreter','latex')
    legend(a1,txt,'FontSize',legend_sz,'interpreter','latex','location','southeast')
    subplot(2,1,2); hold on; 
    a2 = plot(t_mod(2:end),std_x_full(:,2),'-','Color','k','LineWidth',2,'MarkerSize',20);
    plot(t_mod(2:end),-std_x_full(:,2),'-','Color','k','LineWidth',2,'MarkerSize',20)
    xlabel('Time [sec]','Fontsize',xaxis_sz,'interpreter','latex')
    ylabel('v [$\frac{m}{s}$]','Fontsize',yaxis_sz,'interpreter','latex')
    legend(a2,txt,'FontSize',legend_sz,'interpreter','latex','location','southeast')

%     figure; grid on; set(gcf, 'WindowState', 'maximized');
%     subplot(2,1,1); hold on;
%     title('\textbf{Covariances vs. Time}','Fontsize',25,'interpreter','latex')
%     a1 = plot(tspan(2:end),std_innov_full(2:end,1),'-','Color','k','LineWidth',2,'MarkerSize',20);
%     plot(tspan(2:end),-std_innov_full(2:end,1),'-','Color','k','LineWidth',2,'MarkerSize',20)
%     a2 = plot(tspan(2:end),-ones(N-1,1)*sqrt(Pvv),'-','Color','b','LineWidth',2,'MarkerSize',20);
%     plot(tspan(2:end),ones(N-1,1)*sqrt(Pvv),'-','Color','b','LineWidth',2,'MarkerSize',20)
%     ylabel('x [$m$]','Fontsize',yaxis_sz,'interpreter','latex')
%     legend([a1, a2],{'1$\sigma$ innovation cov.','1$\sigma$ measurement noise'},'FontSize',legend_sz,'interpreter','latex','location','southeast')
%     subplot(2,1,2); hold on; 
%     a3 = plot(tspan(2:end),std_innov_full(2:end,2),'-','Color','k','LineWidth',2,'MarkerSize',20);
%     plot(tspan(2:end),-std_innov_full(2:end,2),'-','Color','k','LineWidth',2,'MarkerSize',20)
%     xlabel('Time [sec]','Fontsize',xaxis_sz,'interpreter','latex')
%     ylabel('v [$\frac{m}{s}$]','Fontsize',yaxis_sz,'interpreter','latex')
%     legend(a3,'1$\sigma$ innovation cov.','FontSize',legend_sz,'interpreter','latex','location','southeast')

    figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
    title('\textbf{Covariances vs. Time}','Fontsize',25,'interpreter','latex')
    a1 = plot(tspan(2:end),std_innov_full(2:end,1),'-','Color','k','LineWidth',2,'MarkerSize',20);
    plot(tspan(2:end),-std_innov_full(2:end,1),'-','Color','k','LineWidth',2,'MarkerSize',20)
    a2 = plot(tspan(2:end),-ones(N-1,1)*sqrt(Pvv),'-','Color','b','LineWidth',2,'MarkerSize',20);
    plot(tspan(2:end),ones(N-1,1)*sqrt(Pvv),'-','Color','b','LineWidth',2,'MarkerSize',20)
    ylabel('x [$m$]','Fontsize',yaxis_sz,'interpreter','latex')
    legend([a1, a2],{'1$\sigma$ innovation cov.','1$\sigma$ measurement noise'},'FontSize',legend_sz,'interpreter','latex','location','southeast')

end

%% Functions

function dx = prop(~,x,Fx)
% Propagate state or mean
dx = Fx*x(1:2);

if length(x) > 2
% Reshape cov. matrix
P = reshape(x(3:end)',2,2);
% Prop cov.
dP = Fx*P + P*Fx';
dx = [dx; dP(:)];
end
end
