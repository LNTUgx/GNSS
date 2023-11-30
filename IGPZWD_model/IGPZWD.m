function [pres,ZHD,ZWD,ZTD] = IGPZWD(lat,lon,h,doy,hour,sf,rlevel_coef,clevel_coef)

    % coefficient matrix
    
    % surface
    a0p = rlevel_coef(1:64800,:);   %pressure
    a0w = rlevel_coef(64801:end,:); %ZWD
    
    % vertical
    %pressure
    pp1 = clevel_coef(1:64800,:);    
    pp2 = clevel_coef(64800+1:64800*2,:);
    pp3 = clevel_coef(64800*2+1:64800*3,:);
    %zwd
    ww1 = clevel_coef(64800*3+1:64800*4,:);
    ww2 = clevel_coef(64800*4+1:64800*5,:);
    ww3 = clevel_coef(64800*5+1:64800*6,:);
    
    %interpolation coefficients
    aaa    = abs(lon-(round(lon)-0.5));
    bbb    = abs(lat-(round(lat)-0.5));
    weight = [(1-aaa)*(1-bbb) bbb*(1-aaa) aaa*(1-bbb) aaa*bbb];
    num    = (179.5+(round(lon)-0.5))*180+89.5+(round(lat)-0.5)+1;
    num2   = 90+round(lat);
    
    %periodic terms
    m0 = 1;
    m1 = cos(2*acos(-1.0)*doy/365.25);
    m2 = sin(2*acos(-1.0)*doy/365.25);
    m3 = cos(4*acos(-1.0)*doy/365.25);
    m4 = sin(4*acos(-1.0)*doy/365.25);
    m5 = cos(2*acos(-1.0)*hour/24);
    m6 = sin(2*acos(-1.0)*hour/24);
    m7 = cos(4*acos(-1.0)*hour/24);
    m8 = sin(4*acos(-1.0)*hour/24);

    if (lat>=89.5)
        num=num-1;
        num2=num2-1;
    elseif(lat<=-89.5)
        num=num+1;
        num2=num2+1;
    end

    if (lon>-179.5)&&(lon<179.5)
        nn = [num,num+1,num+180,num+181];
    else
        nn = [num2+180*359,num2+180*359+1,num2,num2+1];
    end

    %v1-v4: values at four grid points
    
    % Pressure for reference layer
    L0v1 = a0p(nn(1),1) + a0p(nn(1),2)*m1 + a0p(nn(1),3)*m2 + a0p(nn(1),4)*m3 + a0p(nn(1),5)*m4+...
                          a0p(nn(1),6)*m5 + a0p(nn(1),7)*m6 + a0p(nn(1),8)*m7 + a0p(nn(1),9)*m8;

    L0v2 = a0p(nn(2),1) + a0p(nn(2),2)*m1 + a0p(nn(2),3)*m2 + a0p(nn(2),4)*m3 + a0p(nn(2),5)*m4+...
                          a0p(nn(2),6)*m5 + a0p(nn(2),7)*m6 + a0p(nn(2),8)*m7 + a0p(nn(2),9)*m8;

    L0v3 = a0p(nn(3),1) + a0p(nn(3),2)*m1 + a0p(nn(3),3)*m2 + a0p(nn(3),4)*m3 + a0p(nn(3),5)*m4+...
                          a0p(nn(3),6)*m5 + a0p(nn(3),7)*m6 + a0p(nn(3),8)*m7 + a0p(nn(3),9)*m8;

    L0v4 = a0p(nn(4),1) + a0p(nn(4),2)*m1 + a0p(nn(4),3)*m2 + a0p(nn(4),4)*m3 + a0p(nn(4),5)*m4+...
                          a0p(nn(4),6)*m5 + a0p(nn(4),7)*m6 + a0p(nn(4),8)*m7 + a0p(nn(4),9)*m8;
    
    % ZWD for reference layer
    wL0v1 = a0w(nn(1),1) + a0w(nn(1),2)*m1 + a0w(nn(1),3)*m2 + a0w(nn(1),4)*m3 + a0w(nn(1),5)*m4 + ...
                           a0w(nn(1),6)*m5 + a0w(nn(1),7)*m6 + a0w(nn(1),8)*m7 + a0w(nn(1),9)*m8;

    wL0v2 = a0w(nn(2),1) + a0w(nn(2),2)*m1 + a0w(nn(2),3)*m2 + a0w(nn(2),4)*m3 + a0w(nn(2),5)*m4 + ...
                           a0w(nn(2),6)*m5 + a0w(nn(2),7)*m6 + a0w(nn(2),8)*m7 + a0w(nn(2),9)*m8;

    wL0v3 = a0w(nn(3),1) + a0w(nn(3),2)*m1 + a0w(nn(3),3)*m2 + a0w(nn(3),4)*m3 + a0w(nn(3),5)*m4 + ...
                           a0w(nn(3),6)*m5 + a0w(nn(3),7)*m6 + a0w(nn(3),8)*m7 + a0w(nn(3),9)*m8;

    wL0v4 = a0w(nn(4),1) + a0w(nn(4),2)*m1 + a0w(nn(4),3)*m2 + a0w(nn(4),4)*m3 + a0w(nn(4),5)*m4 + ...
                           a0w(nn(4),6)*m5 + a0w(nn(4),7)*m6 + a0w(nn(4),8)*m7 + a0w(nn(4),9)*m8;

    %Hight scale factors of pressure
    %v1
    g1r1 = pp1(nn(1),1) + pp1(nn(1),2)*m1 + pp1(nn(1),3)*m2 + pp1(nn(1),4)*m3 + pp1(nn(1),5)*m4; %first order
    g1r2 = pp2(nn(1),1) + pp2(nn(1),2)*m1 + pp2(nn(1),3)*m2 + pp2(nn(1),4)*m3 + pp2(nn(1),5)*m4; %second order
    g1r3 = pp3(nn(1),1) + pp3(nn(1),2)*m1 + pp3(nn(1),3)*m2 + pp3(nn(1),4)*m3 + pp3(nn(1),5)*m4; %third order
    %v2
    g2r1 = pp1(nn(2),1) + pp1(nn(2),2)*m1 + pp1(nn(2),3)*m2 + pp1(nn(2),4)*m3 + pp1(nn(2),5)*m4;
    g2r2 = pp2(nn(2),1) + pp2(nn(2),2)*m1 + pp2(nn(2),3)*m2 + pp2(nn(2),4)*m3 + pp2(nn(2),5)*m4;
    g2r3 = pp3(nn(2),1) + pp3(nn(2),2)*m1 + pp3(nn(2),3)*m2 + pp3(nn(2),4)*m3 + pp3(nn(2),5)*m4;
    %v3
    g3r1 = pp1(nn(3),1) + pp1(nn(3),2)*m1 + pp1(nn(3),3)*m2 + pp1(nn(3),4)*m3 + pp1(nn(3),5)*m4;
    g3r2 = pp2(nn(3),1) + pp2(nn(3),2)*m1 + pp2(nn(3),3)*m2 + pp2(nn(3),4)*m3 + pp2(nn(3),5)*m4;
    g3r3 = pp3(nn(3),1) + pp3(nn(3),2)*m1 + pp3(nn(3),3)*m2 + pp3(nn(3),4)*m3 + pp3(nn(3),5)*m4; 
    %v4
    g4r1 = pp1(nn(4),1) + pp1(nn(4),2)*m1 + pp1(nn(4),3)*m2 + pp1(nn(4),4)*m3 + pp1(nn(4),5)*m4;
    g4r2 = pp2(nn(4),1) + pp2(nn(4),2)*m1 + pp2(nn(4),3)*m2 + pp2(nn(4),4)*m3 + pp2(nn(4),5)*m4;
    g4r3 = pp3(nn(4),1) + pp3(nn(4),2)*m1 + pp3(nn(4),3)*m2 + pp3(nn(4),4)*m3 + pp3(nn(4),5)*m4;  

    %Hight scale factors of ZWD
    %v1
    wg1r1 = ww1(nn(1),1) + ww1(nn(1),2)*m1 + ww1(nn(1),3)*m2 + ww1(nn(1),4)*m3 + ww1(nn(1),5)*m4;%first order
    wg1r2 = ww2(nn(1),1) + ww2(nn(1),2)*m1 + ww2(nn(1),3)*m2 + ww2(nn(1),4)*m3 + ww2(nn(1),5)*m4;%second order
    wg1r3 = ww3(nn(1),1) + ww3(nn(1),2)*m1 + ww3(nn(1),3)*m2 + ww3(nn(1),4)*m3 + ww3(nn(1),5)*m4;%third order
    %v2
    wg2r1 = ww1(nn(2),1) + ww1(nn(2),2)*m1 + ww1(nn(2),3)*m2 + ww1(nn(2),4)*m3 + ww1(nn(2),5)*m4;
    wg2r2 = ww2(nn(2),1) + ww2(nn(2),2)*m1 + ww2(nn(2),3)*m2 + ww2(nn(2),4)*m3 + ww2(nn(2),5)*m4;
    wg2r3 = ww3(nn(2),1) + ww3(nn(2),2)*m1 + ww3(nn(2),3)*m2 + ww3(nn(2),4)*m3 + ww3(nn(2),5)*m4;
    %v3
    wg3r1 = ww1(nn(3),1) + ww1(nn(3),2)*m1 + ww1(nn(3),3)*m2 + ww1(nn(3),4)*m3 + ww1(nn(3),5)*m4;
    wg3r2 = ww2(nn(3),1) + ww2(nn(3),2)*m1 + ww2(nn(3),3)*m2 + ww2(nn(3),4)*m3 + ww2(nn(3),5)*m4;
    wg3r3 = ww3(nn(3),1) + ww3(nn(3),2)*m1 + ww3(nn(3),3)*m2 + ww3(nn(3),4)*m3 + ww3(nn(3),5)*m4; 
    %v4
    wg4r1 = ww1(nn(4),1) + ww1(nn(4),2)*m1 + ww1(nn(4),3)*m2 + ww1(nn(4),4)*m3 + ww1(nn(4),5)*m4;
    wg4r2 = ww2(nn(4),1) + ww2(nn(4),2)*m1 + ww2(nn(4),3)*m2 + ww2(nn(4),4)*m3 + ww2(nn(4),5)*m4;
    wg4r3 = ww3(nn(4),1) + ww3(nn(4),2)*m1 + ww3(nn(4),3)*m2 + ww3(nn(4),4)*m3 + ww3(nn(4),5)*m4;

    %Vertical correction
    v1 = L0v1*exp(g1r1*(h-sf(nn(1))) + g1r2*(h^2-sf(nn(1))^2) + g1r3*(h^3-sf(nn(1))^3)); %
    v2 = L0v2*exp(g2r1*(h-sf(nn(2))) + g2r2*(h^2-sf(nn(2))^2) + g2r3*(h^3-sf(nn(2))^3));
    v3 = L0v3*exp(g3r1*(h-sf(nn(3))) + g3r2*(h^2-sf(nn(3))^2) + g3r3*(h^3-sf(nn(3))^3));
    v4 = L0v4*exp(g4r1*(h-sf(nn(4))) + g4r2*(h^2-sf(nn(4))^2) + g4r3*(h^3-sf(nn(4))^3));

    wv1 = wL0v1*exp(wg1r1*(h-sf(nn(1))) + wg1r2*(h^2-sf(nn(1))^2) + wg1r3*(h^3-sf(nn(1))^3)); %三阶指数外推ZTD
    wv2 = wL0v2*exp(wg2r1*(h-sf(nn(2))) + wg2r2*(h^2-sf(nn(2))^2) + wg2r3*(h^3-sf(nn(2))^3));
    wv3 = wL0v3*exp(wg3r1*(h-sf(nn(3))) + wg3r2*(h^2-sf(nn(3))^2) + wg3r3*(h^3-sf(nn(3))^3));
    wv4 = wL0v4*exp(wg4r1*(h-sf(nn(4))) + wg4r2*(h^2-sf(nn(4))^2) + wg4r3*(h^3-sf(nn(4))^3));

    %Horizontal interpolation
    if (lat>=89.5)
        pres1=(v2-v1)*(lat-89.5)+v2;
        pres2=(v4-v3)*(lat-89.5)+v4;
        pres=(pres2-pres1)*aaa+pres1;
        zwd1=(wv2-wv1)*(lat-89.5)+wv2;
        zwd2=(wv4-wv3)*(lat-89.5)+wv4;
        ZWD=(zwd2-zwd1)*aaa+zwd1;
    elseif(lat<=-89.5)
        pres1=(v1-v2)*abs(lat+89.5)+v1;
        pres2=(v3-v4)*abs(lat+89.5)+v3;
        pres=(pres2-pres1)*aaa+pres1;
        zwd1=(wv1-wv2)*abs(lat+89.5)+wv1;
        zwd2=(wv3-wv4)*abs(lat+89.5)+wv3;
        ZWD=(zwd2-zwd1)*aaa+zwd1;
    else
        pres=v1*weight(1)+v2*weight(2)+v3*weight(3)+v4*weight(4);    %pressure unit: hPa
        ZWD=wv1*weight(1)+wv2*weight(2)+wv3*weight(3)+wv4*weight(4); %zenith wet delay unit: mm
    end

    %Saastamoinen model
    ZHD=pres*2.2768/(1-0.00266*cos(2*lat*pi/180)-0.00000028*h);      %zenith hydrostatic delay unit: mm
    ZTD=ZHD++ZWD;                                                    %zenith total delay unit: mm
end
