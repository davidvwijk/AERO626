% Plotting
xaxis_sz = 20; yaxis_sz = 20; legend_sz = 18;

% Parameters
m = 1.5;  %[kg]
k1 = 2.5; k2 = 3.7; %[N/m]
h = 5.4;  x0 = 3.0; %[m]
v0 = 0.0; %[m/s]

x_star = [4; .2];
del_x0 = [0; 0];
Pxx0_bar = [1000 0; 0 100];
data = load('data_HW03_p1.mat');

omega = sqrt((k1 + k2)/m);
F = [0 1; -omega^2 0];
Pvv = eye(2);

n = 4;
for i = 1:n
    Delta = (Pxx0_bar^-1);
    Lambda = Delta*del_x0;
    t0 = 0;
    h_array = zeros(length(data.hw3_p1_data),2);
    for l = 1:length(data.hw3_p1_data)
        t_l = data.hw3_p1_data(l,1);
        z_l = data.hw3_p1_data(l,2:3)';
        Pvv_l = Pvv;

        Phi = [cos(omega*(t_l-t0)) (1/omega)*(sin(omega*(t_l-t0)));
            -omega*(sin(omega*(t_l-t0))) cos(omega*(t_l-t0))];

        x_l_star = Phi*x_star;

        x = x_l_star(1);
        v = x_l_star(2);
        denom = sqrt(x^2 + h^2);
        H = [x/denom, 0; v/denom - (x^2*v)/(denom^3), x/denom];
        A = H*Phi;
        B = z_l - [denom ; (x*v)/denom];
        Lambda = Lambda + (A)'*(Pvv_l^-1)*(B);
        Delta = Delta + (A)'*(Pvv_l^-1)*(A);
        h_array(l,:) = [denom (x*v)/denom];
    end
    Pxx0 = Delta^-1;
    del_x0_hat = Delta\Lambda;
    x_star = x_star + del_x0_hat;
    del_x0 = del_x0 - del_x0_hat;
    disp(['The reference state after ', num2str(i), ' iterations is:'])
    x_star
    % Plot residuals

    del_z = data.hw3_p1_data(:,2:3) - h_array;
    legend_txt = ['Iteration ',num2str(i)];

    figure(1); hold on; grid on; set(gcf, 'WindowState', 'maximized');
    xlabel('Time [sec]','Fontsize',xaxis_sz,'interpreter','latex')
    ylabel('Range Residuals [m]','Fontsize',yaxis_sz,'interpreter','latex')
    scatter(data.hw3_p1_data(:,1),del_z(:,1),120,'*','DisplayName',legend_txt)
    legend('FontSize',legend_sz,'interpreter','latex','location','southeast')

    figure(2); hold on; grid on; set(gcf, 'WindowState', 'maximized');
    xlabel('Time [sec]','Fontsize',xaxis_sz,'interpreter','latex')
    ylabel('Range-Rate Residuals [m/s]','Fontsize',yaxis_sz,'interpreter','latex')
    scatter(data.hw3_p1_data(:,1),del_z(:,2),120,'*','DisplayName',legend_txt)
    legend('FontSize',legend_sz,'interpreter','latex','location','northeast')
end