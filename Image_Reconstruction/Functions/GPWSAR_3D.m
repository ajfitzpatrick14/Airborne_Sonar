function Im = GPWSAR_3D(CMUTdata, Fs, bistatic, sizeX, sizeY, sizeZ, SOSmap,depths) 
    
    % Time-gain
    dt = 1/Fs;
    t_end = size(CMUTdata,1)*dt;
    t = 0:dt:t_end-dt;
    timeGain = t./t_end;
    CMUTdata = CMUTdata(1:end,:,:).*(timeGain(1:size(CMUTdata,1)))';

    % Medium Properties
    Ca = min(SOSmap,[],'all');
    Cw = max(SOSmap,[],'all');
    zSpeedMax = max(SOSmap,[],[2,3]);
    maxWaveIndex = find(zSpeedMax==Cw,1,'first');
    zSpeedMin = min(SOSmap,[],[2,3]);
    minWaveIndex = find(zSpeedMin==Cw,1,'first');

    % CMUT Measurement Parameters
    [Nt,~,~] = size(CMUTdata);
    
    % Define Computational Domain
    [Nz,Nx,Ny] = size(SOSmap);
    X = linspace(-sizeX/2,sizeX/2,Nx);
    Y = linspace(-sizeY/2,sizeY/2,Ny);
    Z = linspace(0,sizeZ, Nz);
    dx = abs(X(2)-X(1));
    dy = abs(Y(2)-Y(1));
    dz = abs(Z(2)-Z(1));
    [Xgrid,Zgrid,Ygrid] = meshgrid(X,Z,Y);
   
    % Spatial Frequency Variables
    NKx = min(2*2^(nextpow2(Nx)),384);
    NKy = min(2*2^(nextpow2(Ny)),384);
    Kx = [linspace(0,pi/dx,NKx/2) linspace(-pi/dx,0-((pi/dx)/(NKx/2)),NKx/2)];
    Ky = [linspace(0,pi/dy,NKy/2) linspace(-pi/dy,0-((pi/dy)/(NKy/2)),NKy/2)];
    
    % Plane Wave Decomposition
    S_f_Kx_Ky = fftn(CMUTdata,[Nt NKx NKy]);
    freq = 0:Fs/Nt:Fs-Fs/Nt;
    freq = freq(freq<140e3);
    step = floor(750/(freq(2)-freq(1)));
    S_f_Kx_Ky = S_f_Kx_Ky(1:step:length(freq),:,:);
    freq = freq(1:step:end);
    
    % Set up Kza and Kzw tensors
    Ka = 2*pi*freq'/Ca;
    Kw = 2*pi*freq'/Cw;
    [mKx,mKy] = ndgrid(Kx,Ky);
    mKx_f = permute(repmat(mKx,1,1,length(freq)),[3 1 2]);
    mKy_f = permute(repmat(mKy,1,1,length(freq)),[3 1 2]);
    Kza = real(sqrt(Ka.^2-mKx_f.^2-mKy_f.^2));
    Kzw = real(sqrt(Kw.^2-mKx_f.^2-mKy_f.^2));
    
    % Prepare Bistatic Phase Map Correction
    if bistatic(1) == 1
        r = sqrt((Xgrid-bistatic(2)).^2+(Ygrid-bistatic(3)).^2+(Zgrid-bistatic(4)).^2);
    else
        r = zeros(size(SOSmap));
    end
    
    % Migrate to Above Water Surface
    numRows2migrate = maxWaveIndex-1;
    S_update = S_f_Kx_Ky.*exp(1j*Kza*dz*numRows2migrate);
    clear S_f_Kx_Ky
    
    % Migrate through Water Waves
    for di = maxWaveIndex:minWaveIndex-1
        S_air = S_update.*exp(1j*Kza*dz);
        %s_air = ifft2(S_air,[],[2,3]);
        s_air = ifft(ifft(S_air,[],2),[],3);
        S_water = S_update.*exp(1j*Kzw*dz);
        %s_water = ifft2(S_water,[],[2,3]);
        s_water = ifft(ifft(S_water,[],2),[],3);
        A_air = (Cw-SOSmap(di,:,:))/(Cw-Ca);
        A_water = (SOSmap(di,:,:)-Ca)/(Cw-Ca);
        s = s_air(:,1:Nx,1:Ny).*A_air + s_water(:,1:Nx,1:Ny).*A_water;
        %S_update = fft2(s,[NKx,NKy],[2,3]);
        S_update = fft(fft(s,NKx,2),NKy,3);
    end
    
    % Migrate to Water Depths
    Im = zeros(size(SOSmap));
    depthInds = [round(depths(1)/dz) round(depths(2)/dz)] + minWaveIndex;
    depthInds(2) = min(depthInds(2),Nz);
    for di = depthInds(1):depthInds(2)
        inc = di-minWaveIndex-1;
        S = S_update.*Kzw.*exp(1j*Kzw*dz*inc);
        %s = ifft2(S,[],[2,3]);
        s = ifft(ifft(S,[],2),[],3);
        s = s(:,1:Nx,1:Ny).*exp(1j*Kw.*r(di,:,:));
        Im(di,:,:) = sum(s,1);
    end
    
    
end