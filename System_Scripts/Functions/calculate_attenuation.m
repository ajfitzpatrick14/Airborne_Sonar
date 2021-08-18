function [water_atten_dBm, air_atten_dBm] = calculate_attenuation(freq,water_temp,salinity,air_temp,humidity)

    %%%%% Calculation of Water Attenuation %%%%%
    
    % Parameters
    D = 50;
    T = water_temp;
    Tk = 273+T;
    S = salinity;
    pH = 8;

    % Calculation of equation parameters
    c = 1412+3.21*T+1.19*S+0.0167*D;

    A1 = 8.86/c*10^(0.78*pH-5);
    P1 = 1;
    F1 = 2.8*(S/35)^0.5*10^(4-1245/Tk);

    A2 = 21.44*S/c*(1+0.025*T);
    P2 = 1-1.37*1e-4*D+6.2*1e-9*D^2;
    F2 = (8.17*10^(8-1990/Tk))/(1+0.0018*(S-35));

    if T<=20
        A3 = 4.937*1e-4-2.59*1e-5*T+9.11*1e-7*T^2-1.50*1e-8*T^3;
    else
        A3 = 3.964*1e-4-1.146*1e-5*T+1.45*1e-7*T^2-6.5*1e-10*T^3;
    end
    P3 = 1- 3.38*1e-5*D+4.9*1e-10*D^2;

    % Calculation of total absorption (Francios-Garrison)
    F = freq/1e3;
    alpha_dBkm = A1*P1*F1*F^2/(F^2+F1^2) + A2*P2*F2*F^2/(F^2+F2^2) + A3*P3*F^2;

    water_atten_dBm = alpha_dBkm/1e3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% Calculation of Air Attenuation %%%%%
    
    % Parameters
    T_air = air_temp;
    hrar=humidity; 

    % Calculations
    T_0 = 293.15;
    T_01 = 273.16;
    T = T_air + 273.15; 
    p_s0 = 1;

    psat = p_s0*10^(-6.8346*(T_01/T)^1.261 + 4.6151);
    h = p_s0*(hrar)*(psat/p_s0);
    F_rO = 1/p_s0*(24 + 4.04*10^4*h*(0.02+h)/(0.391+h));
    F_rN = 1/p_s0*(T_0/T)^(1/2)*( 9 + 280*h*exp(-4.17*((T_0/T)^(1/3)-1)) );
    alpha_ps= 100*freq^2/p_s0.*( 1.84*10^(-11)*(T/T_0)^(1/2)...
        + (T/T_0)^(-5/2)*(0.01275*exp(-2239.1/T)/(F_rO + freq^2/F_rO)...
        + 0.1068*exp(-3352/T)/(F_rN + freq^2/F_rN) ) );
    a_ps_ar = alpha_ps*20/log(10);
    
    air_atten_dBm = a_ps_ar/100; % dB per meter
    
end


