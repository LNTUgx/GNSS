function [ZTD] = EGZTD(lat,lon,h,doy,hour,sf,rlevel_coef,clevel_coef)

    % coefficient matrix
    
    % surface ZTD
    a0ZTD = rlevel_coef;
    
    % Height scale-factor
    % Exponential function
    ef1 = clevel_coef(1:64800,:); 
    ef2 = clevel_coef(64800+1:64800*2,:);
    
    % Quadratic polynomial   
    qp1 = clevel_coef(64800*2+1:64800*3,:);
    qp2 = clevel_coef(64800*3+1:64800*4,:);
    
    % interpolation coefficients
    aaa    = abs(lon-(round(lon)-0.5));
    bbb    = abs(lat-(round(lat)-0.5));
    weight = [(1-aaa)*(1-bbb) bbb*(1-aaa) aaa*(1-bbb) aaa*bbb];
    num    = (179.5+(round(lon)-0.5))*180+89.5+(round(lat)-0.5)+1;
    num2   = 90+round(lat);
    
    % periodic terms
    % m0 = 1;
    m1 = cos(2*acos(-1.0)*doy/365.25);
    m2 = sin(2*acos(-1.0)*doy/365.25);
    m3 = cos(4*acos(-1.0)*doy/365.25);
    m4 = sin(4*acos(-1.0)*doy/365.25);
    m5 = cos(2*acos(-1.0)*hour/24);
    m6 = sin(2*acos(-1.0)*hour/24);

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
    
    % ZTD prediction for reference layer
    L0v1 = a0ZTD(nn(1),1) + a0ZTD(nn(1),2)*m1 + a0ZTD(nn(1),3)*m2 + a0ZTD(nn(1),4)*m3 + a0ZTD(nn(1),5)*m4+...
                           a0ZTD(nn(1),6)*m5 + a0ZTD(nn(1),7)*m6;

    L0v2 = a0ZTD(nn(2),1) + a0ZTD(nn(2),2)*m1 + a0ZTD(nn(2),3)*m2 + a0ZTD(nn(2),4)*m3 + a0ZTD(nn(2),5)*m4+...
                           a0ZTD(nn(2),6)*m5 + a0ZTD(nn(2),7)*m6;

    L0v3 = a0ZTD(nn(3),1) + a0ZTD(nn(3),2)*m1 + a0ZTD(nn(3),3)*m2 + a0ZTD(nn(3),4)*m3 + a0ZTD(nn(3),5)*m4+...
                           a0ZTD(nn(3),6)*m5 + a0ZTD(nn(3),7)*m6;

    L0v4 = a0ZTD(nn(4),1) + a0ZTD(nn(4),2)*m1 + a0ZTD(nn(4),3)*m2 + a0ZTD(nn(4),4)*m3 + a0ZTD(nn(4),5)*m4+...
                           a0ZTD(nn(4),6)*m5 + a0ZTD(nn(4),7)*m6;

    % HSF of ZTD
    % g1
    g1_e1 = ef1(nn(1),1) + ef1(nn(1),2)*m1 + ef1(nn(1),3)*m2 + ef1(nn(1),4)*m3 + ef1(nn(1),5)*m4;  % ¦Â1 of EFO1
    g1_e2 = ef2(nn(1),1) + ef2(nn(1),2)*m1 + ef2(nn(1),3)*m2 + ef2(nn(1),4)*m3 + ef2(nn(1),5)*m4;  % ¦Â2 of EFO1
    g1_p1 = qp1(nn(1),1) + qp1(nn(1),2)*m1 + qp1(nn(1),3)*m2 + qp1(nn(1),4)*m3 + qp1(nn(1),5)*m4;  % ¦Â1 of QP
    g1_p2 = qp2(nn(1),1) + qp2(nn(1),2)*m1 + qp2(nn(1),3)*m2 + qp2(nn(1),4)*m3 + qp2(nn(1),5)*m4;  % ¦Â2 of QP

        if (lat<60 && lat>-30) % The semi-annual term is considered at these latitude zones
            g1_e1_2 = ef1(nn(1),1) + ef1(nn(1),2)*m1 + ef1(nn(1),3)*m2 + ef1(nn(1),4)*m3 + ef1(nn(1),5)*m4;
            g1_p1_2 = qp1(nn(1),1) + qp1(nn(1),2)*m1 + qp1(nn(1),3)*m2 + qp1(nn(1),4)*m3 + qp1(nn(1),5)*m4;
            g1_e2_2 = ef2(nn(1),1) + ef2(nn(1),2)*m1 + ef2(nn(1),3)*m2 + ef2(nn(1),4)*m3 + ef2(nn(1),5)*m4;
            g1_p2_2 = qp2(nn(1),1) + qp2(nn(1),2)*m1 + qp2(nn(1),3)*m2 + qp2(nn(1),4)*m3 + qp2(nn(1),5)*m4;
        else
            g1_e1_2 = ef1(nn(1),1) + ef1(nn(1),2)*m1 + ef1(nn(1),3)*m2;   
            g1_p1_2 = qp1(nn(1),1) + qp1(nn(1),2)*m1 + qp1(nn(1),3)*m2; 
            g1_e2_2 = ef2(nn(1),1) + ef2(nn(1),2)*m1 + ef2(nn(1),3)*m2;   
            g1_p2_2 = qp2(nn(1),1) + qp2(nn(1),2)*m1 + qp2(nn(1),3)*m2; 
        end
        % The semi-annual term is not considered
            g1_e1_3 = ef1(nn(1),1) + ef1(nn(1),2)*m1 + ef1(nn(1),3)*m2;
            g1_p1_3 = qp1(nn(1),1) + qp1(nn(1),2)*m1 + qp1(nn(1),3)*m2;
            g1_e2_3 = ef2(nn(1),1) + ef2(nn(1),2)*m1 + ef2(nn(1),3)*m2;
            g1_p2_3 = qp2(nn(1),1) + qp2(nn(1),2)*m1 + qp2(nn(1),3)*m2;
   
    % g2
    g2_e1 = ef1(nn(2),1) + ef1(nn(2),2)*m1 + ef1(nn(2),3)*m2 + ef1(nn(2),4)*m3 + ef1(nn(2),5)*m4;  % ¦Â1 of EFO1
    g2_e2 = ef2(nn(2),1) + ef2(nn(2),2)*m1 + ef2(nn(2),3)*m2 + ef2(nn(2),4)*m3 + ef2(nn(2),5)*m4;  % ¦Â2 of EFO1
    g2_p1 = qp1(nn(2),1) + qp1(nn(2),2)*m1 + qp1(nn(2),3)*m2 + qp1(nn(2),4)*m3 + qp1(nn(2),5)*m4;  % ¦Â1 of QP
    g2_p2 = qp2(nn(2),1) + qp2(nn(2),2)*m1 + qp2(nn(2),3)*m2 + qp2(nn(2),4)*m3 + qp2(nn(2),5)*m4;  % ¦Â2 of QP

        if (lat<60 && lat>-30) % The semi-annual term is considered at these latitude zones
            g2_e1_2 = ef1(nn(2),1) + ef1(nn(2),2)*m1 + ef1(nn(2),3)*m2 + ef1(nn(2),4)*m3 + ef1(nn(2),5)*m4;
            g2_p1_2 = qp1(nn(2),1) + qp1(nn(2),2)*m1 + qp1(nn(2),3)*m2 + qp1(nn(2),4)*m3 + qp1(nn(2),5)*m4;
            g2_e2_2 = ef2(nn(2),1) + ef2(nn(2),2)*m1 + ef2(nn(2),3)*m2 + ef2(nn(2),4)*m3 + ef2(nn(2),5)*m4;
            g2_p2_2 = qp2(nn(2),1) + qp2(nn(2),2)*m1 + qp2(nn(2),3)*m2 + qp2(nn(2),4)*m3 + qp2(nn(2),5)*m4;
        else
            g2_e1_2 = ef1(nn(2),1) + ef1(nn(2),2)*m1 + ef1(nn(2),3)*m2;   
            g2_p1_2 = qp1(nn(2),1) + qp1(nn(2),2)*m1 + qp1(nn(2),3)*m2; 
            g2_e2_2 = ef2(nn(2),1) + ef2(nn(2),2)*m1 + ef2(nn(2),3)*m2;   
            g2_p2_2 = qp2(nn(2),1) + qp2(nn(2),2)*m1 + qp2(nn(2),3)*m2; 
        end
        % The semi-annual term is not considered
            g2_e1_3 = ef1(nn(2),1) + ef1(nn(2),2)*m1 + ef1(nn(2),3)*m2;
            g2_p1_3 = qp1(nn(2),1) + qp1(nn(2),2)*m1 + qp1(nn(2),3)*m2;
            g2_e2_3 = ef2(nn(2),1) + ef2(nn(2),2)*m1 + ef2(nn(2),3)*m2;
            g2_p2_3 = qp2(nn(2),1) + qp2(nn(2),2)*m1 + qp2(nn(2),3)*m2;
    
    % g3
    g3_e1 = ef1(nn(3),1) + ef1(nn(3),2)*m1 + ef1(nn(3),3)*m2 + ef1(nn(3),4)*m3 + ef1(nn(3),5)*m4;  % ¦Â1 of EFO1
    g3_e2 = ef2(nn(3),1) + ef2(nn(3),2)*m1 + ef2(nn(3),3)*m2 + ef2(nn(3),4)*m3 + ef2(nn(3),5)*m4;  % ¦Â2 of EFO1
    g3_p1 = qp1(nn(3),1) + qp1(nn(3),2)*m1 + qp1(nn(3),3)*m2 + qp1(nn(3),4)*m3 + qp1(nn(3),5)*m4;  % ¦Â1 of QP
    g3_p2 = qp2(nn(3),1) + qp2(nn(3),2)*m1 + qp2(nn(3),3)*m2 + qp2(nn(3),4)*m3 + qp2(nn(3),5)*m4;  % ¦Â2 of QP

        if (lat<60 && lat>-30) % The semi-annual term is considered at these latitude zones
            g3_e1_2 = ef1(nn(3),1) + ef1(nn(3),2)*m1 + ef1(nn(3),3)*m2 + ef1(nn(3),4)*m3 + ef1(nn(3),5)*m4;
            g3_p1_2 = qp1(nn(3),1) + qp1(nn(3),2)*m1 + qp1(nn(3),3)*m2 + qp1(nn(3),4)*m3 + qp1(nn(3),5)*m4;
            g3_e2_2 = ef2(nn(3),1) + ef2(nn(3),2)*m1 + ef2(nn(3),3)*m2 + ef2(nn(3),4)*m3 + ef2(nn(3),5)*m4;
            g3_p2_2 = qp2(nn(3),1) + qp2(nn(3),2)*m1 + qp2(nn(3),3)*m2 + qp2(nn(3),4)*m3 + qp2(nn(3),5)*m4;
        else
            g3_e1_2 = ef1(nn(3),1) + ef1(nn(3),2)*m1 + ef1(nn(3),3)*m2;   
            g3_p1_2 = qp1(nn(3),1) + qp1(nn(3),2)*m1 + qp1(nn(3),3)*m2; 
            g3_e2_2 = ef2(nn(3),1) + ef2(nn(3),2)*m1 + ef2(nn(3),3)*m2;   
            g3_p2_2 = qp2(nn(3),1) + qp2(nn(3),2)*m1 + qp2(nn(3),3)*m2; 
        end
        % The semi-annual term is not considered
            g3_e1_3 = ef1(nn(3),1) + ef1(nn(3),2)*m1 + ef1(nn(3),3)*m2;
            g3_p1_3 = qp1(nn(3),1) + qp1(nn(3),2)*m1 + qp1(nn(3),3)*m2;
            g3_e2_3 = ef2(nn(3),1) + ef2(nn(3),2)*m1 + ef2(nn(3),3)*m2;
            g3_p2_3 = qp2(nn(3),1) + qp2(nn(3),2)*m1 + qp2(nn(3),3)*m2;
    
    % g4
    g4_e1 = ef1(nn(4),1) + ef1(nn(4),2)*m1 + ef1(nn(4),3)*m2 + ef1(nn(4),4)*m3 + ef1(nn(4),5)*m4;  % ¦Â1 of EFO1
    g4_e2 = ef2(nn(4),1) + ef2(nn(4),2)*m1 + ef2(nn(4),3)*m2 + ef2(nn(4),4)*m3 + ef2(nn(4),5)*m4;  % ¦Â2 of EFO1
    g4_p1 = qp1(nn(4),1) + qp1(nn(4),2)*m1 + qp1(nn(4),3)*m2 + qp1(nn(4),4)*m3 + qp1(nn(4),5)*m4;  % ¦Â1 of QP
    g4_p2 = qp2(nn(4),1) + qp2(nn(4),2)*m1 + qp2(nn(4),3)*m2 + qp2(nn(4),4)*m3 + qp2(nn(4),5)*m4;  % ¦Â2 of QP

        if (lat<60 && lat>-30) % The semi-annual term is considered at these latitude zones
            g4_e1_2 = ef1(nn(4),1) + ef1(nn(4),2)*m1 + ef1(nn(4),3)*m2 + ef1(nn(4),4)*m3 + ef1(nn(4),5)*m4;
            g4_p1_2 = qp1(nn(4),1) + qp1(nn(4),2)*m1 + qp1(nn(4),3)*m2 + qp1(nn(4),4)*m3 + qp1(nn(4),5)*m4;
            g4_e2_2 = ef2(nn(4),1) + ef2(nn(4),2)*m1 + ef2(nn(4),3)*m2 + ef2(nn(4),4)*m3 + ef2(nn(4),5)*m4;
            g4_p2_2 = qp2(nn(4),1) + qp2(nn(4),2)*m1 + qp2(nn(4),3)*m2 + qp2(nn(4),4)*m3 + qp2(nn(4),5)*m4;
        else
            g4_e1_2 = ef1(nn(4),1) + ef1(nn(4),2)*m1 + ef1(nn(4),3)*m2;   
            g4_p1_2 = qp1(nn(4),1) + qp1(nn(4),2)*m1 + qp1(nn(4),3)*m2; 
            g4_e2_2 = ef2(nn(4),1) + ef2(nn(4),2)*m1 + ef2(nn(4),3)*m2;   
            g4_p2_2 = qp2(nn(4),1) + qp2(nn(4),2)*m1 + qp2(nn(4),3)*m2; 
        end
        % The semi-annual term is not considered
            g4_e1_3 = ef1(nn(4),1) + ef1(nn(4),2)*m1 + ef1(nn(4),3)*m2;
            g4_p1_3 = qp1(nn(4),1) + qp1(nn(4),2)*m1 + qp1(nn(4),3)*m2;
            g4_e2_3 = ef2(nn(4),1) + ef2(nn(4),2)*m1 + ef2(nn(4),3)*m2;
            g4_p2_3 = qp2(nn(4),1) + qp2(nn(4),2)*m1 + qp2(nn(4),3)*m2;
    
    % Vertical correction
    g1(1,:) = L0v1*exp(g1_e1*(h-sf(nn(1))) + g1_e2*(h^2-sf(nn(1))^2) );                  % version#1
    g1(2,:) = L0v1*exp(g1_e1_2*(h-sf(nn(1))) + g1_e2_2*(h^2-sf(nn(1))^2) );                % version#2
    g1(3,:) = L0v1*exp(g1_e1_3*(h-sf(nn(1))) + g1_e2_3*(h^2-sf(nn(1))^2) );                % version#3
    g1(4,:) = g1_p1*(h-sf(nn(1))) + g1_p2*(h^2-sf(nn(1))^2) + L0v1;                      % version#4
    g1(5,:) = g1_p1_2*(h-sf(nn(1))) + g1_p2_2*(h^2-sf(nn(1))^2) + L0v1;                    % version#5
    g1(6,:) = g1_p1_3*(h-sf(nn(1))) + g1_p2_3*(h^2-sf(nn(1))^2) + L0v1;                    % version#6

    g2(1,:) = L0v2*exp(g2_e1*(h-sf(nn(2))) + g2_e2*(h^2-sf(nn(2))^2) );                  % version#1
    g2(2,:) = L0v2*exp(g2_e1_2*(h-sf(nn(2))) + g2_e2_2*(h^2-sf(nn(2))^2) );                % version#2
    g2(3,:) = L0v2*exp(g2_e1_3*(h-sf(nn(2))) + g2_e2_3*(h^2-sf(nn(2))^2) );                % version#3
    g2(4,:) = g2_p1*(h-sf(nn(2))) + g2_p2*(h^2-sf(nn(2))^2) + L0v2;                      % version#4
    g2(5,:) = g2_p1_2*(h-sf(nn(2))) + g2_p2_2*(h^2-sf(nn(2))^2) + L0v2;                    % version#5
    g2(6,:) = g2_p1_3*(h-sf(nn(2))) + g2_p2_3*(h^2-sf(nn(2))^2) + L0v2;                    % version#6
    
    g3(1,:) = L0v3*exp(g3_e1*(h-sf(nn(3))) + g3_e2*(h^2-sf(nn(3))^2) );                  % version#1
    g3(2,:) = L0v3*exp(g3_e1_2*(h-sf(nn(3))) + g3_e2_2*(h^2-sf(nn(3))^2) );                % version#2
    g3(3,:) = L0v3*exp(g3_e1_3*(h-sf(nn(3))) + g3_e2_3*(h^2-sf(nn(3))^2) );                % version#3
    g3(4,:) = g3_p1*(h-sf(nn(3))) + g3_p2*(h^2-sf(nn(3))^2) + L0v3;                      % version#4
    g3(5,:) = g3_p1_2*(h-sf(nn(3))) + g3_p2_2*(h^2-sf(nn(3))^2) + L0v3;                    % version#5
    g3(6,:) = g3_p1_3*(h-sf(nn(3))) + g3_p2_3*(h^2-sf(nn(3))^2) + L0v3;                    % version#6
    
    g4(1,:) = L0v4*exp(g4_e1*(h-sf(nn(4))) + g4_e2*(h^2-sf(nn(4))^2) );                  % version#1
    g4(2,:) = L0v4*exp(g4_e1_2*(h-sf(nn(4))) + g4_e2_2*(h^2-sf(nn(4))^2) );                % version#2
    g4(3,:) = L0v4*exp(g4_e1_3*(h-sf(nn(4))) + g4_e2_3*(h^2-sf(nn(4))^2) );                % version#3
    g4(4,:) = g4_p1*(h-sf(nn(4))) + g4_p2*(h^2-sf(nn(4))^2) + L0v4;                      % version#4
    g4(5,:) = g4_p1_2*(h-sf(nn(4))) + g4_p2_2*(h^2-sf(nn(4))^2) + L0v4;                    % version#5
    g4(6,:) = g4_p1_3*(h-sf(nn(4))) + g4_p2_3*(h^2-sf(nn(4))^2) + L0v4;                    % version#6
    
    % Horizontal interpolation
    if (lat>=89.5)
        ztd1=(g2-g1)*(lat-89.5)+g2;
        ztd2=(g4-g3)*(lat-89.5)+g4;
        ZTD=(ztd2-ztd1)*aaa+ztd1;
    elseif(lat<=-89.5)
        ztd1=(g1-g2)*abs(lat+89.5)+g1;
        ztd2=(g3-g4)*abs(lat+89.5)+g3;
        ZTD=(ztd2-ztd1)*aaa+ztd1;
    else
        ZTD=g1*weight(1)+g2*weight(2)+g3*weight(3)+g4*weight(4);    % ZTD unit: mm
    end
end
