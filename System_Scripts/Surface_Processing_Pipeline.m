%% Post-Processing of Raw Surface Mapping

clear;
clc;

% Add folder to path
addpath(genpath('Data'))

numScansX = 54;
numScansY = 50;

for xrec_num = 1:numScansX
        
    % Load Raw Depth Sensor Point Cloud
    file2load = sprintf('raw_point_cloud_%d.mat',xrec_num);
    load(file2load);

    % Define Computational Domain in X,Y
    sizeX = 0.26;
    sizeY = 0.24;
    dr = 1135*4.21875e-6/4;
    X = linspace(-sizeX/2,sizeX/2,round(sizeX/dr));
    Y = linspace(-sizeY/2,sizeY/2,round(sizeY/dr));
    
    pointsX = vertices(:,:,1);
    pointsY = vertices(:,:,2);
    pointsZ = vertices(:,:,3);
    clear vertices;
    surface_maps = zeros(size(pointsX,1),length(X),length(Y));
    count_maps = zeros(size(pointsX,1),length(X),length(Y));

    for yrec_num = 1:numScansY

        % What is offset in X,Y of camera?
        camera_offset_X = 0.038;
        camera_offset_Y = 0.02;
        [offset_map_Y, offset_map_X] = meshgrid(linspace(-sizeY/2,sizeY/2-sizeY/numScansY,numScansY),linspace(-sizeX/2,sizeX/2-sizeX/numScansX,numScansX));
        offset_X = offset_map_X(xrec_num,yrec_num)+camera_offset_X;
        offset_Y = offset_map_Y(xrec_num,yrec_num)+camera_offset_Y;

        % Convert Point Clouds to Surface Maps Defined over Domain
        fn = yrec_num;
        for i = 1:size(pointsX,2)
            if abs(pointsX(fn,i)+offset_Y)<=sizeY/2 && abs(pointsY(fn,i)+offset_X)<=sizeX/2
                [~,indY] = min(abs(Y-(pointsX(fn,i)+offset_Y)));
                [~,indX] = min(abs(X-(pointsY(fn,i)+offset_X)));
                surface_maps(fn,indX,indY) = surface_maps(fn,indX,indY)+pointsZ(fn,i);
                count_maps(fn,indX,indY) = count_maps(fn,indX,indY)+1;
            end
        end
    end
    count_maps(count_maps==0) = 1;
    surface_maps = surface_maps./count_maps;
    clear count_maps;

    % Surface Interpolation - removing holes and discontinuities
    window_size = round(length(X)/4);
    surf_fit = zeros(size(surface_maps));
    for fn = 1:size(surface_maps,1)
        surf_temp = filloutliers(surface_maps(fn,:,:),'nearest');
        surface = filloutliers(squeeze(surf_temp),'nearest','movmedian',window_size,'ThresholdFactor',1);
        surf_fit(fn,:,:) = filloutliers(surface,'nearest','movmedian',window_size,2,'ThresholdFactor',1);
        Check_outlier = isoutlier(squeeze(surf_fit(fn,:,:)),'movmedian',window_size,1,'ThresholdFactor',1.5);
        flag = any(Check_outlier,'all');
        if(flag)
            surf_fit(fn,:,:) = filloutliers(squeeze(surf_fit(fn,:,:)),'nearest','movmedian',window_size,1,'ThresholdFactor',1.5);
        end    
    end
    
    % Surface Filtering - removing high spatial frequencies
    surface_maps_smooth = zeros(size(surface_maps));
    for fn = 1:size(surface_maps,1)
        surface = squeeze(surf_fit(fn,:,:));

        [M, N] = size(surface);
        DimPadded = max(2^nextpow2(M), 2^nextpow2(N)); 
        MtoPad = round((DimPadded - M)/2);
        NtoPad = round((DimPadded - N)/2);
   
        surfacePadded = padarray(surface,[MtoPad NtoPad],'replicate','both');
        if size(surfacePadded,1)~=DimPadded(1) || size(surfacePadded,2)~=DimPadded(2)
            surfacePadded = surfacePadded(1:DimPadded,1:DimPadded);
        end

        fftSurfacePadded = fftshift(fft2(surfacePadded));
    
        % Low pass Butterworth Filter
        LamdaC = 0.02; % 2cm wavelength waves as cutoff
        % Formula : deltaF (resolution of frequency bin) = 1/(dr*DimPadded) and Fc = 1/LamdaC
        fcBins = dr*DimPadded/(LamdaC);
        X_grid = (0:DimPadded-1);
        Y_grid = (0:DimPadded-1);
        [X_grid, Y_grid] = meshgrid(X_grid,Y_grid);
        Cx = 0.5*DimPadded;
        Cy = 0.5*DimPadded;

        radius = sqrt((X_grid-Cx).^2 + (Y_grid-Cy).^2);
        filtOrder = 4;
        ButterworthLowPassFilter = 1./(1.0 + (radius./fcBins).^(2*filtOrder));
    
        filtered = fftSurfacePadded.*ButterworthLowPassFilter;
        smoothSurface = real(ifft2(ifftshift(filtered)));
    
        surface_maps_smooth(fn,:,:) = smoothSurface(MtoPad+1:MtoPad+M,NtoPad+1:NtoPad+N);
    end

    clear surface fitProfileCol fitProfileRow f x;

    file2save = sprintf('surface_maps_rowX%d',xrec_num);
    save(file2save,'surface_maps_smooth')
end

%% Display Post-Processed Surface Maps

for xrec_num = 1:numScansX
        
    % Load Raw Depth Sensor Point Cloud
    file2load = sprintf('surface_maps_rowX%d.mat',xrec_num);
    load(file2load);

    for fn = 1:size(surface_maps_smooth,1)
            surf(-squeeze(surface_maps_smooth(fn,1:end,1:end)));
            zlim([-0.61 -0.55])
            view([-156 45])
            colormap(flip(jet))
            yticks([9,34,59,84,109,134,159,184,209])
            yticklabels([12, 9, 6, 3, 0, -3, -6, -9, -12])
            xticks([26,51,76,101,126,151,176])
            xticklabels(flip([-9, -6, -3, 0, 3, 6, 9]))
            zticks(flip([-0.54 -0.56 -0.58 -0.60 -0.62]))
            zticklabels(flip([54 56 58 60 62]))
            xlim([1 200])
            ylim([1 216])
            set(gca,'fontsize',22)
            xlabel('Y (cm)','fontsize',24)
            ylabel('X (cm)','fontsize',24)
            zlabel('Z (cm)', 'fontsize',24)
            getframe;  
    end

end