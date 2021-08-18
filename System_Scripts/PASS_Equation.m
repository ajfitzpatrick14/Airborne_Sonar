%% Photoacoustic Airborne SONAR System (PASS) Range Equation
% Author: Aidan Fitzpatrick

% Note: This script can be used to analyze the maximum imaging depth of PASS.

% User-Defined Parameters:
    % P0: Laser peak power in Watts
    % laser_wavelength_nm: Laser wavelength in nm (can take 'Optimum')
    % freq: Acoustic frequency in Hz (i.e. laser intensity modulation frequency)
    % NEP: Noise Equivalent Pressure of the CMUT receiver in Pascal (noisefloor*sqrt(bandwidth))
    % TS: Target strength in dB
    % PCG: Pulse-Compression Gain (currently not exploited, i.e. = 0 dB)
    % N: Total number of CMUTs in the array
    % CMUT_height: Height of the CMUT receivers above water in meters
    
% References:
% Media Properties
    % Feistel, Rainer. "APPENDIX TO EOLSS ARTICLE 02-03-07 THERMODYNAMIC PROPERTIES OF SEAWATER."
% Optical Absorption Coefficient
    % Hale, George M., and Marvin R. Querry. "Optical constants of water in the 200-nm to 200-?m wavelength region." Applied optics 12.3 (1973): 555-563.
% Acoustic Attenuation
    % Francois, R. E., and G. R. Garrison. "Sound absorption based on ocean measurements. Part II: Boric acid contribution and equation for total absorption." The Journal of the Acoustical Society of America 72.6 (1982): 1879-1890.
    % Bass, Henry E., et al. "Atmospheric absorption of sound: Further developments." The Journal of the Acoustical Society of America 97.1 (1995): 680-683.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% User-Defined Parameters
P0 = 1e3:100:100e3;
laser_wavelength_nm = 'Optimum';
freq = 50e3;
NEP = 100e-6;
TS = -40:0.1:20;
PCG = 0;
N = 1024; 
CMUT_height = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Functions')

%%%% Media Parameters %%%%
% Thermal Expansion Coefficient of Water
beta = 2.5e-4;
% Specific Heat Capacity of Water
Cp = 4000;
% Optical Transmission Coefficient
T = 0.98;
% Speed-of-sound, mass density, and acoustic impedance
cw = 1500;
ca = 343;
rhow = 1000;
rhoa = 1.225;
Zair = rhoa*ca;
Zwater = rhow*cw;
% Acoustic Wavenumber in Water
k = 2*pi*freq/cw;
% Optical Absorption Coefficient
u = calculate_absorption(laser_wavelength_nm,k); % assumes pure water
% Acoustic Transmission Loss from Water to Air
tau = 1-((Zair-Zwater)/(Zair+Zwater))^2;
% Acoustic Attenuation in dB/m in Water and Air
% Function Arguments: frequency, water temp. (celcius), water salinity, air temp. (celcius), air humidity percentage
[alpha_w, alpha_a] = calculate_attenuation(freq,10,35,20,60); 

maxDepth = zeros(length(TS),length(P0));
for t = 1:length(TS)
    for p = 1:length(P0)
        for z = 1:1:400
            % Source Level (SL)
            P = T*2*pi*freq*beta*P0(p)*u*k/(2*Cp*(u^2+k^2));
            SL = 10*log10(P^2/Zwater);
            % Transmission Loss - Incident (TL_i)
            TL_i = alpha_w*z+10*log10(z^2);
            % Transmission Loss - Reflected (TL_r)  
            TL_r = alpha_w*z+alpha_a*CMUT_height+10*log10(z^2)-10*log10(tau);
            % Noise Level (NL)
            NL = 10*log10(NEP^2/Zair);
            % Array Gain (AG)
            AG = 10*log10(N);
            % SNR Calculation
            SNR = SL - TL_i + TS(t) - TL_r - NL + AG + PCG;
            if SNR>10
                continue;
            else
                maxDepth(t,p) = z;
                break;
            end
        end   
    end
end
 
figure;
contourf((maxDepth),25)
set(gca,'fontsize',18)
h = colorbar;
set(get(h,'label'),'string','Maximum Imaging Depth: D (m)','fontsize',20);
ylabel('Target Strength: TS (dB)','fontsize',20)
xlabel('Laser Power (kW)','fontsize',20)
colormap(jet(25)) 
xticks([91 191 291 391 491 591 691 791 891 991])
xticklabels([10 20 30 40 50 60 70 80 90 100])
yticks([1 21 41 61 81 101 121 141]*4.95)
yticklabels([-40 -30 -20 -10 0 10 20])
caxis([0 250])

%%

clear;
clc;

% User-Defined Parameters
P0 = 50e3;
laser_wavelength_nm = 'Optimum';
freq = 50e3;
NEP = 100e-6;
TS = [-20 -10 0 10 20];
maxDepth = 100;
PCG = 0;
N = 1024; 
CMUT_height = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Functions')

%%%% Media Parameters %%%%
% Thermal Expansion Coefficient of Water
beta = 2.5e-4;
% Specific Heat Capacity of Water
Cp = 4000;
% Optical Transmission Coefficient
T = 0.98;
% Speed-of-sound, mass density, and acoustic impedance
cw = 1500;
ca = 343;
rhow = 1000;
rhoa = 1.225;
Zair = rhoa*ca;
Zwater = rhow*cw;
% Acoustic Wavenumber in Water
k = 2*pi*freq/cw;
% Optical Absorption Coefficient
u = calculate_absorption(laser_wavelength_nm,k); % assumes pure water
% Acoustic Transmission Loss from Water to Air
tau = 1-((Zair-Zwater)/(Zair+Zwater))^2;
% Acoustic Attenuation in dB/m in Water and Air
% Function Arguments: frequency, water temp. (celcius), water salinity, air temp. (celcius), air humidity percentage
[alpha_w, alpha_a] = calculate_attenuation(freq,10,35,20,60); 

depth = 1:2.5:maxDepth;
SNR = zeros(length(TS),length(depth));
figure;
markers = ['o','+','*','d','s'];
hold on;
for t = 1:length(TS)
    for z = 1:length(depth)  
        % Source Level (SL)
        P = T*2*pi*freq*beta*P0*u*k/(2*Cp*(u^2+k^2));
        SL = 10*log10(P^2/Zwater);
        % Transmission Loss - Incident (TL_i)
        TL_i = alpha_w*depth(z)+10*log10(depth(z)^2);
        % Transmission Loss - Reflected (TL_r)  
        TL_r = alpha_w*depth(z)+alpha_a*CMUT_height+10*log10(depth(z)^2)-10*log10(tau);
        % Noise Level (NL)
        NL = 10*log10(NEP^2/Zair);
        % Array Gain (AG)
        AG = 10*log10(N);
        % SNR Calculation
        SNR(t,z) = SL - TL_i + TS(t) - TL_r - NL + AG + PCG;
    end
    legend_entry = sprintf('TS = %d dB',TS(t));
    plot(depth,SNR(t,:),'linewidth',2, 'DisplayName',legend_entry,'Marker',markers(t),'MarkerSize',8)
end

critical_point = 10*ones(1,length(depth));
plot(depth,critical_point, 'k','linewidth',1.5, 'DisplayName', 'Detection Threshold', 'linestyle','--')
legend('show','fontsize',18)
ylabel('SNR (dB)', 'fontsize',20)
xlabel('Target Depth (m)', 'fontsize',20)
set(gca,'fontsize',18)
%title('SNR vs. Target Depth','fontsize',20)
grid on;
box on;
xlim([1 maxDepth])
ylim([-20 100])
