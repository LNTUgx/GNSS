clear
clc
% ------------------------------------------------------%
    % input parameters:
    % 
    % lat :  latitude in degree [-90:90] 
    % lon :  longitude in degree [-180:180] 
    % h   :  ellipsoidal height in meter 
    % doy :  day of year 
    % hod :  hour of day 
    % 
    % output parameters:
    % 
    % pres:  atmospheric pressure in hPa 
    % ZHD :  zenith hydrostatic delay in mm 
    % ZWD :  zenith wet delay in mm 
    % ZTD :  zenith total delay in mm 
% ------------------------------------------------------%

% Model height: ellipsoidal height
load("surf_height.mat"); %
sf=re(:,1);

% coefficient matrix
p_a0=load("P_surface");
w_a0=load("ZWD_surface");
p_coef1=load("htx1");
p_coef2=load("htx2");
p_coef3=load("htx3");
w_coef1=load("wtx1");
w_coef2=load("wtx2");
w_coef3=load("wtx3");
rlevel_coef=[p_a0;w_a0];
clevel_coef=[p_coef1;p_coef2;p_coef3;w_coef1;w_coef2;w_coef3];

% example
% lat  = 16.26;  % unit: degree
% lon  = -61.51; % unit: degree
% h    = 141.0;  % unit: meter  
% doy  = 20;       
% hour = 5;    

% prediction              
[pres,ZHD,ZWD,ZTD] = IGPZWD(lat,lon,h,doy,hour,sf,rlevel_coef,clevel_coef);
 


