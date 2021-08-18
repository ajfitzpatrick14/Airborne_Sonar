%% k-Wave Simulation for Surface Mapping Accuracy
clear
clc

% Using a GPU?
GPU = false;

% Create the k-Wave Grid
points_per_wavelength = 4;
f_max = 71e3;
c_min = 340;
dx = c_min/(points_per_wavelength*f_max);
Nx = 608;
Nz = 608;
x_size = Nx*dx;
z_size = Nz*dx;
kgrid = kWaveGrid(Nz,dx,Nx,dx);

% Define Medium Parameters
height = 0.2;
Ca = 340;
pa = 1.225;
Cw = 1500;
pw = 1000;
medium.sound_speed = Cw*ones(Nz,Nx);
medium.density = pw*ones(Nz,Nx);

% Define Targets
disc1_z_pos = round((0.40+height)/dx);    
disc1_x_pos = round(Nx/2+0.1/dx);        
disc2_z_pos = round((0.48+height)/dx);    
disc2_x_pos = round(Nx/2);        
disc3_z_pos = round((0.42+height)/dx);    
disc3_x_pos = round(Nx/2-0.1/dx);         
disc_radius = round(0.013/dx);    
disc_1 = makeDisc(Nz, Nx, disc1_z_pos, disc1_x_pos, disc_radius);
disc_2 = makeDisc(Nz, Nx, disc2_z_pos, disc2_x_pos, disc_radius);
disc_3 = makeDisc(Nz, Nx, disc3_z_pos, disc3_x_pos, disc_radius);     

% Surface Generation - REPLACE WITH SURFACE_GENERATOR.m FUNCTION
height = 0.2;
lambda = 0.05:0.02:0.2;
k = 2*pi./lambda;
theta = 2*pi*rand(1,length(k));
phases = 2*pi*rand(1,length(k));
% 2-D wave
[Y,X] = meshgrid(-x_size/2:x_size/Nx:x_size/2,-x_size/2:x_size/Nx:x_size/2);
H2 = zeros(size(X));
for m = 1:length(k)
        H2 = H2 + lambda(m)*cos(k(m)*X*cos(theta(m))+k(m)*Y*sin(theta(m))+phases(m));
end    
% Create surface map of desired amplitude
amplitude = 0.025;
surface_map = rescale(H2,height-amplitude,height+amplitude);
surface_map = surface_map(1:Nx,Nx/2); 

% incude ocean waves
for j = 1:Nx
   i_dx = round(surface_map(j)/dx);
   for i = 1:Nz
       if i<i_dx
           medium.sound_speed(i,j) = Ca;
           medium.density(i,j) = pa;
       end
   end
end
SOSmap = medium.sound_speed;
medium.sound_speed(disc_1==1) = 2000;                                                         
medium.density(disc_1==1) = 5000;                  
medium.sound_speed(disc_2==1) = 2000;                                                            
medium.density(disc_2==1) = 5000;                 
medium.sound_speed(disc_3==1) = 2000;                                                            
medium.density(disc_3==1) = 5000;  

figure;
imagesc(medium.sound_speed)

% Define Source
j = Nx/2;
source_ind = round(surface_map(j)/dx);
source.p_mask = zeros(Nz, Nx);
source.p_mask(source_ind+2:source_ind+7,Nx/2) = 1;
Fs = 25e6;
freq = 71e3;
te = 0:1/Fs:3/freq;
ti = 0:1/Fs:2/freq;
t = 0:1/Fs:2e-3-1/Fs;
kgrid.t_array = t;
signal = 100*[sin(2*pi*freq*te) zeros(1,round(0.5*Fs/freq)) 1.059*sin(2*pi*freq*ti)];
signal = [signal zeros(1,length(t)-length(signal))];
filtered_signal = filterTimeSeries(kgrid, medium, signal);
source.p = filtered_signal;
figure; plot(t,filtered_signal)

% Define Sensor
sensor.mask = zeros(Nz,Nx);
rec_mask = 2:4:Nx;
sensor.mask(1,rec_mask) = 1;

% Run Simulation
input_args = {'PMLInside', false, 'DataCast', 'single', 'DataRecast', false,...
'PlotSim', true, 'PlotLayout', true, 'PlotPML', false};
if GPU
    %run the simulation using GPU
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
else
    %run the simulation using CPU
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
end

% CMUT Impulse Response
Q = 28; 
wc= 2*pi*71e3;
delw= wc/Q;
H = tf([delw 0], [1 delw wc^2]);

% Apply CMUT Response
signals = cast(sensor_data,'double');
CMUTdata = zeros(length(t),size(signals,1));
for ri = 1:size(signals,1)
    CMUTdata(:,ri) = lsim(H,signals(ri,:),t);
end

% Matched Filter Data
Fs_down = 2e6;
matched_filter = lsim(H,filtered_signal,t);
CMUTdata_in = preprocess(CMUTdata,Fs,Fs_down,matched_filter);


%%% Image Reconstruction %%%
addpath(genpath('Image_Reconstruction'))

% Define Computational Domain
array_size = 0.5;
sizeZ = z_size;
sizeX = x_size;
bistatic = [1 0 0 source_ind*dx];

% Run the Reconstruction Algorithm
array_size = min(array_size,sizeX);
mask = round((Nx/2-array_size/2/dx)/4):round((Nx/2+array_size/2/dx)/4);
rec_mask2 = rec_mask(mask);
CMUTdata_in2 = zeros(size(CMUTdata_in,1),Nx);
CMUTdata_in2(:,rec_mask2) = CMUTdata_in(:,mask);
Reconstructed_Image_2D = GPWSAR_2D(CMUTdata_in2, Fs_down, bistatic, sizeX, sizeZ, SOSmap,[0 sizeZ]);
final_image = Reconstructed_Image_2D;
   
% Display the Simulation Environment
cmap = [1 1 1; 0.8706 0.9216 0.9647; 0 0 0; 1 0.1765 0];
SOSmap(SOSmap==340) = 0;
SOSmap(SOSmap==1500) = 1;
SOSmap(medium.sound_speed==2000) = 2;
SOSmap(1:source_ind,Nx/2-2:Nx/2+2) = 3;
figure; 
imagesc(SOSmap)
colormap(cmap)
set(gca,'fontsize',18)
xlabel('X (cm)','fontsize',20)
ylabel('Z (cm)','fontsize',20)
xticks([57 138 219 300 381 462 543])
xticklabels([-30 -20 -10 0 10 20 30])
yticks([1 81 162 243 324 405 486 567])
yticklabels([0 10 20 30 40 50 60 70])
axis image;

% Display the Reconstructed Image
final_image(1:source_ind+round(0.03/dx),:) = 0;
im2Display = abs(final_image).^1.4;
im2Display = im2Display./(max(im2Display(source_ind+round(0.1/dx):end,:),[],'all'));
figure;
imagesc(im2Display)
set(gca,'fontsize',18)
axis image;
xlabel('X (cm)','fontsize',20)
ylabel('Z (cm)','fontsize',20)
xticks([57 138 219 300 381 462 543])
xticklabels([-30 -20 -10 0 10 20 30])
yticks([1 81 162 243 324 405 486 567])
yticklabels([0 10 20 30 40 50 60 70])
caxis([0 1]) 
cbh = colorbar;
cbh.Ticks = linspace(0, 1,6); 

