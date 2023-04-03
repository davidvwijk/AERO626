%% AERO 626 Homework #2
% Spring 2023
% David van Wijk

data = load('data_HW02.mat');
format long

% DATA PROVIDED ARE:
% T = (m x 1) array of measurement times [s]
% Z = (m x 1) array of position measurements [m]
% W = (m x 1) array of measurement weights [nd]
% R = (m x 1) array of measurement noise covariances [m^2]

%% Part A: Least-Squares Estimate of State

F = [0 1; -1 0];
H_tilde = [1 0];
H = [];
for i = 1:length(data.T)
    Phi_i = expm(F*(data.T(i) - data.T(1)));
    H_i = H_tilde*Phi_i;
    H = [H; H_i];
end
disp('Least-Squares Estimate of Initial State using matrix exponential for Phi:')
x_hat_0 = (H'*H)\(H'*data.Z)

H = [];
for i = 1:length(data.T)
    t_i = data.T(i);
    Phi_i = [cos(t_i) sin(t_i); -sin(t_i) cos(t_i)];
    H_i = H_tilde*Phi_i;
    H = [H; H_i];
end
disp('Least-Squares Estimate of Initial State using analytical solution for Phi:')
x_hat_0 = (H'*H)\(H'*data.Z)

%% Part B: Weighted Least-Squares Estimate of State

W = diag(data.W);
disp('Weighted Least-Squares Estimate of Initial State:')
x_hat_0 = (H'*W*H)\(H'*W*data.Z)

%% Part C: Weighted Least-Squares Estimate with prior information
x_bar = [1; 0];
W_bar = [3 0; 0 3];
disp('Weighted Least-Squares Estimate of Initial State using prior info:')
x_hat_0 = (H'*W*H + W_bar)\(H'*W*data.Z + W_bar*x_bar)

%% Part D: LUMVE of State

disp('LUMVE of Initial State:')
P_vv = diag(data.R);
x_hat_0 = (H'*P_vv^-1*H)\(H'*P_vv^-1*data.Z)
disp('The uncertainty in our measurement can be evaluated using the covariance matrix of the estimate:')
cov_matrix = (H'*P_vv^-1*H)^-1

