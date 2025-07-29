clear
clc
% ------------------------------------------------------%
% EGZTD model
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
    % ZTD :  zenith total delay in mm 
% ------------------------------------------------------%

% Model height: ellipsoidal height
load("surface_height.mat"); %
sf=re(:,1);

% Coefficient matrix of the EGZTD model for six schemes
ZTD_a0 = load("ZTD_7pre");
coef1  = load("ctx1");
coef2  = load("ctx2");
coef3  = load("ptx1");
coef4  = load("ptx2");

rlevel_coef = ZTD_a0;
clevel_coef = [coef1;coef2;coef3;coef4];

% Example
lat  = 25.3;  % unit: degree
lon  = -60.7; % unit: degree
eh   = 251.0;  % unit: meter  
doy  = 120;       
hod  = 3;    

% Prediction              
ZTD6 = EGZTD(lat,lon,eh,doy,hod,sf,rlevel_coef,clevel_coef); 


