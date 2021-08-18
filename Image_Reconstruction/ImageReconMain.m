%% Image Reconstruction of Experimental & Simulated Data

clear;
clc;

% Add Data Files to Path
addpath(genpath('Data'))
addpath(genpath('Data'))
addpath(genpath('Image_Reconstruction'))

% User-Defined Inputs
in = input('Enter Case Number: ');
% 1 = 'S' target in Hydrostatic Conditions
% 2 = 'U' target in Hydrostatic Conditions
% 3 = 'U' target in Hydrodynamic Conditions
% 4 = simulated with small waves
% 5 = simulated with large waves

switch in
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 1
        dataFile = 'S_data_static.mat';
        load(dataFile);
        
        %%% File Contents %%%
        % CMUTdata: matched filtered acoustic measurement data
        % Fs: sampling frequency of measurement data
        % bistatic: [1, laser offset in x, laser offset in y, cmut height]
        % sizeX: size of reconstruction grid in X
        % sizeY: size of reconstruction grid in Y
        % sizeZ: size of reconstruction grid in Z
        % SOSmap: spatial distribution of the speed of sound
        % depth_range: depth range to reconstruct
        
        %%% Image Reconstruction Algorithm %%%
        Reconstructed_Image = GPWSAR_3D(CMUTdata,Fs,bistatic,sizeX,sizeY,sizeZ,SOSmap,depth_range);
        
        %%% Figure Display %%%
        depth2view = 0.17;
        dr = sizeY/size(CMUTdata,3);
        ind = round((bistatic(4)+depth2view)/dr);
        im2Display = squeeze(sum(abs(Reconstructed_Image(ind-2:ind+2,:,:)),1)).^1.4;
        im2Display = rescale(imresize(im2Display,4*size(im2Display)));
        figure; 
        imagesc(im2Display(10:180,15:170))
        colormap(copper)
        axis image
        cc = colorbar;
        cc.YTick = [0.2 0.4 0.6 0.8 1];
        xticks([15 36 57 78 99 120 141])
        xticklabels([-75 -50 -25 0 25 50 75])
        yticks([1 22 43 64 85 106 127 148 169])
        yticklabels([-100 -75 -50 -25 0 25 50 75 100])
        set(gca,'fontsize',18)
        ylabel('Y (mm)','fontsize',20)
        xlabel('X (mm)','fontsize',20)
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2
        dataFile = 'U_data_static.mat';
        load(dataFile);
        
        %%% File Contents %%%
        % CMUTdata: matched filtered acoustic measurement data
        % Fs: sampling frequency of measurement data
        % bistatic: [1, hyd. offset in x, hyd. offset in y, cmut height]
        % sizeX: size of reconstruction grid in X
        % sizeY: size of reconstruction grid in Y
        % sizeZ: size of reconstruction grid in Z
        % SOSmap: spatial distribution of the speed of sound
        % depth_range: depth range to reconstruct   
        
        %%% Image Reconstruction Algorithm %%%
        Reconstructed_Image = GPWSAR_3D(CMUTdata,Fs,bistatic,sizeX,sizeY,sizeZ,SOSmap,depth_range);
        
        %%% Figure Display %%%
        depth2View = 0.17;
        dr = sizeY/size(CMUTdata,3);
        ind = round((bistatic(4)-0.04+depth2View)/dr); % hyd. 4 cm underwater
        im2Display = (squeeze(sum(abs(Reconstructed_Image(ind-2:ind+2,:,:)),1)).^1.4)';
        im2Display = rescale(imresize(im2Display,2*size(im2Display)));
        figure; 
        imagesc(im2Display(1:120,16:end))
        colormap(copper)
        cc = colorbar;
        cc.YTick = [0.2 0.4 0.6 0.8 1];
        axis image
        yticks([7 33 60 87 113])
        yticklabels([-8 -4 0 4 8])
        xticks([8 41 75 108 141 175 208])
        xticklabels([-12 -8 -4 0 4 8 12])
        set(gca,'fontsize',18)
        ylabel('Y (cm)','fontsize',20)
        xlabel('X (cm)','fontsize',20) 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        dataFile = 'U_data_dynamic2.mat';
        load(dataFile);
        
        %%% File Contents %%%
        % CMUTdata: matched filtered acoustic measurement data
        % Fs: sampling frequency of measurement data
        % bistatic: [1, hyd. offset in x, hyd. offset in y, cmut height]
        % sizeX: size of reconstruction grid in X
        % sizeY: size of reconstruction grid in Y
        % sizeZ: size of reconstruction grid in Z
        % depth_range: depth range to reconstruct
        % ind: index of the depth slice that corresponds to 18cm depth
        
        %%% Loop over measurements for image reconstruction %%%
        num = 1;
        Image2Save = zeros(217,200,size(CMUTdata,2)*size(CMUTdata,3));
        for xi = 1:size(CMUTscan,2)
            
            if xi == 1
                ti = tic();
            end
            
            % Load Surface Maps
            file2load = sprintf('surface_maps_rowX%d.mat',xi);
            load(file2load);
            surface_maps_smooth = surface_maps_smooth+0.005;
        
            for yi = 1:size(CMUTscan,3)
                
                % Convert Surface Map to Spatial Distribution of SoS
                dr = 340/71e3/4;
                SOS = [340 1500];
                SOSmap = SOS(1)*ones(round(sizeZ/dr),size(surface_maps_smooth,2),size(surface_maps_smooth,3));
                for k = floor(min(surface_maps_smooth(:))/dr)-5:ceil(max(surface_maps_smooth(:))/dr)+5
                    for i = 1:size(SOSmap,2)
                        for j = 1:size(SOSmap,3)
                            if k>surface_maps_smooth(yi,i,j)/dr
                                SOSmap(k,i,j) = SOS(2);
                            end
                        end
                    end
                end
                SOSmap(k:end,:,:) = SOS(2);
                
                % Which CMUT measurement to use
                CMUTdata_in = zeros(size(CMUTdata,1),size(SOSmap,2),size(SOSmap,3));
                CMUTlocX = 2+4*(xi-1);
                CMUTlocY = 2+4*(yi-1);
                CMUTdata_in(:,CMUTlocX,CMUTlocY) = CMUTdata(:,xi,yi);
                
                % Image Reconstruction Algorithm
                Reconstructed_Image = GPWSAR_3D(CMUTdata_in,Fs,bistatic,sizeX,sizeY,sizeZ,SOSmap,depth_range);
                
                % Depth Slice to Save
                Image2Save(:,:,num) = Reconstructed_Image(ind,:,:);
                num = num+1;
                
                if xi==1 && yi==1
                    stop_time = toc(ti);
                    fprintf('Estimated Image Reconstruction Time: %.2f hours \n', stop_time*2700/3600);
                end
                
            end
            
        end
        
        %%% Figure Display %%%
        im2Display = rescale((abs(squeeze(sum((Image2Save),3))))');
        figure; 
        imagesc(im2Display(1:120,:).^1.4)
        axis equal tight
        colormap(copper)
        cc = colorbar;
        cc.YTick = [0.2 0.4 0.6 0.8 1];
        axis equal tight
        %yticks([1 33 67 100 133 167 200])
        %yticklabels([-12 -8 -4 0 4 8 12])
        yticks([7 33 60 87 113])
        yticklabels([-8 -4 0 4 8])
        xticks([8 41 75 108 141 175 208])
        xticklabels([-12 -8 -4 0 4 8 12])
        set(gca,'fontsize',18)
        ylabel('Y (cm)','fontsize',20)
        xlabel('X (cm)','fontsize',20)
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 4
        dataFile = 'small_waves_sim.mat';
        load(dataFile);
        
        %%% File Contents %%%
        % CMUTdata: matched filtered acoustic measurement data
        % Fs: sampling frequency of measurement data
        % bistatic: [1, hyd. offset in x, hyd. offset in y, cmut height]
        % sizeX: size of reconstruction grid in X
        % sizeZ: size of reconstruction grid in Z
        % SOSmap: spatial distribution of the speed of sound
        % depth_range: depth range to reconstruct
        % cmap: colormap used to display simulation environment
        % sim_env: depiction of the simulation environment
        % rec_mask: location of the CMUT receivers
        
        % Pick the Array Size
        array_size = 0.3;
        Nx = size(SOSmap,2);
        dx = sizeX/Nx;
        array_size = min(array_size,sizeX);
        mask = round((Nx/2-array_size/2/dx)/4):round((Nx/2+array_size/2/dx)/4);
        rec_mask2 = rec_mask(mask);
        CMUTdata_in = zeros(size(CMUTdata,1),Nx);
        CMUTdata_in(:,rec_mask2) = CMUTdata(:,mask);
        
        %%% Image Reconstruction Algorithm %%%
        Reconstructed_Image = GPWSAR_2D(CMUTdata_in,Fs,bistatic,sizeX,sizeZ,SOSmap,depth_range);
        
        %%% Figure Display %%%
        figure; 
        subplot(1,2,1)
        imagesc(sim_env)
        colormap(gca,cmap)
        set(gca,'fontsize',14)
        xlabel('X (cm)','fontsize',16)
        ylabel('Z (cm)','fontsize',16)
        xticks([57 138 219 300 381 462 543])
        xticklabels([-30 -20 -10 0 10 20 30])
        yticks([1 81 162 243 324 405])
        yticklabels([0 10 20 30 40 50])
        axis image;
        subplot(1,2,2)
        im2Display = abs(Reconstructed_Image).^1.4;
        im2Display = im2Display./(max(im2Display(270:end,:),[],'all'));
        imagesc(im2Display)
        set(gca,'fontsize',14)
        axis image;
        xlabel('X (cm)','fontsize',16)
        ylabel('Z (cm)','fontsize',16)
        xticks([57 138 219 300 381 462 543])
        xticklabels([-30 -20 -10 0 10 20 30])
        yticks([1 81 162 243 324 405])
        yticklabels([0 10 20 30 40 50])
        caxis([0 1])
        pos = get(gca,'Position');
        cbh = colorbar('location','SouthOutside','Position', [pos(1) pos(2)+0.1 pos(3) 0.03]);
        cbh.Ticks = linspace(0,1,6);  
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 5
        dataFile = 'large_waves_sim.mat';
        load(dataFile);
        
        %%% File Contents %%%
        % CMUTdata: matched filtered acoustic measurement data
        % Fs: sampling frequency of measurement data
        % bistatic: [1, hyd. offset in x, hyd. offset in y, cmut height]
        % sizeX: size of reconstruction grid in X
        % sizeZ: size of reconstruction grid in Z
        % SOSmap: spatial distribution of the speed of sound
        % depth_range: depth range to reconstruct
        % cmap: colormap used to display simulation environment
        % sim_env: depiction of the simulation environment
        % rec_mask: location of the CMUT receivers
        
        % Pick the Array Size
        array_size = 0.3;
        Nx = size(SOSmap,2);
        dx = sizeX/Nx;
        array_size = min(array_size,sizeX);
        mask = round((Nx/2-array_size/2/dx)/4):round((Nx/2+array_size/2/dx)/4);
        rec_mask2 = rec_mask(mask);
        CMUTdata_in = zeros(size(CMUTdata,1),Nx);
        CMUTdata_in(:,rec_mask2) = CMUTdata(:,mask);
        
        %%% Image Reconstruction Algorithm %%%
        Reconstructed_Image = GPWSAR_2D(CMUTdata_in,Fs,bistatic,sizeX,sizeZ,SOSmap,depth_range);
        
        %%% Figure Display %%%
        figure; 
        subplot(1,2,1)
        imagesc(sim_env)
        colormap(gca,cmap)
        set(gca,'fontsize',14)
        xlabel('X (cm)','fontsize',16)
        ylabel('Z (cm)','fontsize',16)
        xticks([57 138 219 300 381 462 543])
        xticklabels([-30 -20 -10 0 10 20 30])
        yticks([1 81 162 243 324 405])
        yticklabels([0 10 20 30 40 50])
        axis image;
        subplot(1,2,2)
        im2Display = abs(Reconstructed_Image).^1.4;
        im2Display = im2Display./(max(im2Display(270:end,:),[],'all'));
        imagesc(im2Display)
        set(gca,'fontsize',14)
        axis image;
        xlabel('X (cm)','fontsize',16)
        ylabel('Z (cm)','fontsize',16)
        xticks([57 138 219 300 381 462 543])
        xticklabels([-30 -20 -10 0 10 20 30])
        yticks([1 81 162 243 324 405])
        yticklabels([0 10 20 30 40 50])
        caxis([0 1])
        pos = get(gca,'Position');
        cbh = colorbar('location','SouthOutside','Position', [pos(1) pos(2)+0.1 pos(3) 0.03]);
        cbh.Ticks = linspace(0,1,6);   
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



