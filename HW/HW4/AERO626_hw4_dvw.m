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
std_x_full = zeros(N*2,2);
std_innov_full = zeros(N,2);

count = 0;
for i = 1:N
    % Propagate the truth state
    if i == 1
        xk = x0;
    else
        [~,X] = ode45(@(t,x) prop(t,x,Fx),[tspan(i-1),tspan(i)],xkm1,opts);
        xk = X(end,1:2)';
    end
    
    % Generate measurement
    zk = Hx*xk + chol(Pvv)'*randn;

    % Propagate step of Kalman filter
    if i == 1
        mxkm = mx0;
        Pxxkm = Pxx0;
    else
        [~,X] = ode45(@(t,x) prop(t,x,Fx),[tspan(i-1),tspan(i)],[mxkm1;Pxxkm1(:)],opts);
        mxkm = X(end,1:2)';
        Pxxkm = reshape(X(end,3:end)',2,2);
    end

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

disp(['The final standard deviation of the position is: ', num2str(std_x_full(end,1)),' [m]'])
disp(['The final standard deviation of the velocity is: ', num2str(std_x_full(end,2)),' [m/s]'])
disp(['The final standard deviation of the innovation covariance is: ', num2str(std_innov_full(end,1)),' [m]'])

if plot_flag 
    xaxis_sz = 20; yaxis_sz = 20; legend_sz = 18;
    
    % Figure 1
    txt = '1$\sigma$';
    opts_pts1 = {80,'b',"o"};
    opts_pts2 = {80,'r',"o"};

    figure; grid on; set(gcf, 'WindowState', 'maximized');
    subplot(2,1,1); hold on; grid on;
    title('\textbf{State Covariance vs. Time}','Fontsize',25,'interpreter','latex')
    t_mod = sort([tspan tspan]);
    a1 = plot(t_mod,std_x_full(:,1),'-','Color','k','LineWidth',2,'MarkerSize',20);
    b1 = scatter(tspan,std_x_full(1:2:end,1),opts_pts1{:});
    b2 = scatter(tspan,std_x_full(2:2:end,1),opts_pts2{:});
    plot(t_mod,-std_x_full(:,1),'-','Color','k','LineWidth',2,'MarkerSize',20)
    scatter(tspan,-std_x_full(1:2:end,1),opts_pts1{:})
    scatter(tspan,-std_x_full(2:2:end,1),opts_pts2{:})
    ylabel('x [$m$]','Fontsize',yaxis_sz,'interpreter','latex')
    legendtxt = {txt, 'before update', 'after update'};
    legend([a1 b1 b2],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')

    subplot(2,1,2); hold on; grid on;
    a2 = plot(t_mod,std_x_full(:,2),'-','Color','k','LineWidth',2,'MarkerSize',20);
    b3 = scatter(tspan,std_x_full(1:2:end,2),opts_pts1{:});
    b4 = scatter(tspan,std_x_full(2:2:end,2),opts_pts2{:});
    scatter(tspan,-std_x_full(1:2:end,2),opts_pts1{:})
    scatter(tspan,-std_x_full(2:2:end,2),opts_pts2{:})
    plot(t_mod,-std_x_full(:,2),'-','Color','k','LineWidth',2,'MarkerSize',20)
    xlabel('Time [sec]','Fontsize',xaxis_sz,'interpreter','latex')
    ylabel('v [$\frac{m}{s}$]','Fontsize',yaxis_sz,'interpreter','latex')
    legend([a2 b3 b4],legendtxt,'FontSize',legend_sz,'interpreter','latex','location','southeast')
    
    % Figure 2
    figure; grid on; set(gcf, 'WindowState', 'maximized'); hold on;
    title('\textbf{Covariances vs. Time}','Fontsize',25,'interpreter','latex')
    a1 = plot(tspan,std_innov_full(:,1),'-','Color','k','LineWidth',2,'MarkerSize',20);
    plot(tspan,-std_innov_full(:,1),'-','Color','k','LineWidth',2,'MarkerSize',20)
    a2 = plot(tspan,-ones(N,1)*sqrt(Pvv),'-','Color','b','LineWidth',2,'MarkerSize',20);
    plot(tspan,ones(N,1)*sqrt(Pvv),'-','Color','b','LineWidth',2,'MarkerSize',20)
    ylabel('x [$m$]','Fontsize',yaxis_sz,'interpreter','latex')
    xlabel('Time [sec]','Fontsize',xaxis_sz,'interpreter','latex')
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
