clear all
close all
clc

% THIS IS A TEMPLATE FOR THE EKF IN HW #)5
% YOU WILL NEED TO FILL IN QUITE A FEW MISSING ITEMS
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

% constants
GM = 3.986004415e5;         % [km^3/s^2]

% conversion factors
asc2deg = 1.0/3600.0;       % [arcsec] -> [deg]
deg2rad = pi/180.0;         % [deg] -> [rad]
asc2rad = asc2deg*deg2rad;  % [arcsec] -> [rad]
rad2asc = 3600.0*180.0/pi;  % [rad] -> [arcsec]
rad2deg = 180.0/pi;         % [rad] -> [deg]

% load data file
load('data_HW05','T','Z','Pvv','Robsv','Xtrue')
% DATA PROVIDED ARE:
%   T     = (m x 1) array of observation times [sec]
%   Z     = (2 x m) array of right-ascension and declination observations [arcsec]
%   Pvv   = (2 x 2 x m) array of measurement noise covariances [arcsec^2]
%   Robsv = (3 x m) array of observer position in inertial frame [km]
%   Xtrue = (6 x m) array of true object position and velocity [km] and [km/s]
% m = 1199 for this dataset
% there are 123, 118, 101, 128, 106, 137, 108, 110, 134, 134 measurements on
% each of the ten consecutive nights, respectively

% given initial conditions
t0   = 0.0;
mx0  = 
Pxx0 = 

% declare storage space for saving state estimation error information
xcount  = 0;
xstore  = nan(6,2*length(T)-1);
txstore = nan(1,2*length(T)-1);
kxstore = nan(1,2*length(T)-1);
mxstore = nan(6,2*length(T)-1);
exstore = nan(6,2*length(T)-1);
sxstore = nan(6,2*length(T)-1);

% declare storage space
zcount  = 0;
zstore  = nan(2,length(T));
tzstore = nan(1,length(T));
kzstore = nan(1,length(T));
mzstore = nan(2,length(T));
ezstore = nan(2,length(T));
dzstore = nan(1,length(T));
szstore = nan(2,length(T));
svstore = nan(2,length(T));

% store initial data
xcount            = xcount + 1;
txstore(:,xcount) = T(1);
kxstore(:,xcount) = 0;
mxstore(:,xcount) = mx0;
sxstore(:,xcount) = sqrt(diag(Pxx0));

% define the process noise mapping matrix and the process noise PSD
Fw  = 
Qww = 

% initialize time, mean, and covariance for the EKF
tkm1   = 
mxkm1  = 
Pxxkm1 = 

% loop over the number of data points
for k = 1:length(T)
    tk   = T(k);       % time of the current measurement to process [sec]
    zk   = Z(:,k);     % current measurement to process [arcsec]
    Pvvk = Pvv(:,:,k); % measurement noise covariance [arcsec^2]
    rk   = Robsv(:,k); % inertial position of the observer [km]

    % unpack the truth -- this cannot be used in the filter, only for analysis
    xk = Xtrue(:,k);   % true state of the object [km and km/s]
    
    % propagate the mean and covariance
    
    
    
    
    

    % store a priori state information for analysis
    xcount            = xcount + 1;
    xstore(:,xcount)  = xk;
    txstore(:,xcount) = tk;
    kxstore(:,xcount) = k;
    mxstore(:,xcount) = mxkm;
    exstore(:,xcount) = xk - mxkm;
    sxstore(:,xcount) = sqrt(diag(Pxxkm));
    
    % compute the estimated measurement
    %   form the relative position, then compute the expected measurement
    %   change the units of the reference measurement to [arcsec]
    
    
    
    
    
    
    
    
    
    
    % compute the measurement Jacobian
    %   change the units of the Jacobian to [arcsec]
    
          
    
    
    % update the mean and covariance
    
    
    
    
    
    
    
    % store a posteriori state information for analysis
    xcount            = xcount + 1;
    xstore(:,xcount)  = Xtrue(:,k);
    txstore(:,xcount) = tk;
    kxstore(:,xcount) = k;
    mxstore(:,xcount) = mxkp;
    exstore(:,xcount) = Xtrue(:,k) - mxkp;
    sxstore(:,xcount) = sqrt(diag(Pxxkp));
    
    % store measurement information for analysis
    zcount            = zcount + 1;
    zstore(:,zcount)  = zk;
    tzstore(:,zcount) = tk;
    kzstore(:,zcount) = k;
    mzstore(:,zcount) = mzkm;
    ezstore(:,zcount) = zk - mzkm;
    dzstore(:,zcount) = (zk - mzkm)'*(Pzzkm\(zk - mzkm));
    szstore(:,zcount) = sqrt(diag(Pzzkm));
    svstore(:,zcount) = sqrt(diag(Pvvk));
    
    % cycle the time, mean, and covariance for the next step of the EKF
    tkm1   = 
    mxkm1  = 
    Pxxkm1 = 
end