function Im = GPWSAR_2D(CMUTdata, Fs, bistatic, sizeX, sizeZ, SOSmap,depths) 

    % Time-gain
    dt = 1/Fs;
    t_end = size(CMUTdata,1)*dt;
    t = 0:dt:t_end-dt;
    timeGain = t./t_end;
    CMUTdata = CMUTdata(1:end,:).*sqrt(timeGain(1:size(CMUTdata,1)))';
    
    % Medium Properties
    Ca = min(SOSmap,[],'all');
    Cw = max(SOSmap,[],'all');
    zSpeedMax = max(SOSmap,[],2);
    maxWaveIndex = find(zSpeedMax==Cw,1,'first');
    zSpeedMin = min(SOSmap,[],2);
    minWaveIndex = find(zSpeedMin==Cw,1,'first');

    % CMUT Measurement Parameters
    [Nt,~] = size(CMUTdata);
    
    % Define Computational Domain
    [Nz,Nx] = size(SOSmap);
    X = linspace(-sizeX/2,sizeX/2,Nx);
    Z = linspace(0,sizeZ, Nz);
    dx = abs(X(2)-X(1));
    dz = abs(Z(2)-Z(1));
    [Xgrid,Zgrid] = meshgrid(X,Z);
   
    % Spatial Frequency Variables
    NKx = 2*2^(nextpow2(Nx));
    Kx = [linspace(0,pi/dx,NKx/2) linspace(-pi/dx,0-((pi/dx)/(NKx/2)),NKx/2)];
    
    % Plane Wave Decomposition
    S_f_Kx = fft2(CMUTdata,Nt,NKx);
    freq = 0:Fs/Nt:Fs-Fs/Nt;
    freq = freq(freq<240e3);
    step = 1;%round(800/(freq(2)-freq(1)));
    S_f_Kx = S_f_Kx(1:step:length(freq),:);
    freq = freq(1:step:end);
    
    % Set up Kza and Kzw tensors
    Ka = 2*pi*freq'/Ca;%repmat(2*pi*freq'/Ca,1,NKx);
    Kw = 2*pi*freq'/Cw;%repmat(2*pi*freq'/Cw,1,NKx);
    mKx = repmat(Kx,length(freq),1);
    Kza = real(sqrt(Ka.^2-mKx.^2));
    Kzw = real(sqrt(Kw.^2-mKx.^2));  
    
    % Prepare Bistatic Phase Map Correction
    if bistatic(1) == 1
        r = sqrt((Xgrid-bistatic(2)).^2+bistatic(3)^2+(Zgrid-bistatic(4)).^2);
    else
        r = zeros(size(SOSmap));
    end
    
    % Migrate to Above Water Surface
    numRows2migrate = maxWaveIndex-1;
    S_update = S_f_Kx.*exp(1j*Kza*dz*numRows2migrate);
    
    % Migrate through Water Waves
    for di = maxWaveIndex:minWaveIndex-1
        S_air = S_update.*exp(1j*Kza*dz);
        s_air = ifft(S_air,[],2);
        S_water = S_update.*exp(1j*Kzw*dz);
        s_water = ifft(S_water,[],2);
        A_air = (Cw-SOSmap(di,:))/(Cw-Ca);
        A_water = (SOSmap(di,:)-Ca)/(Cw-Ca);
        s = s_air(:,1:Nx).*A_air + s_water(:,1:Nx).*A_water;
        S_update = fft(s,NKx,2);
    end
    
    % Migrate to Water Depths
    Im = zeros(size(SOSmap));
    depthInds = [round(depths(1)/dz) round(depths(2)/dz)] + minWaveIndex;
    depthInds(2) = min(depthInds(2),Nz);
    for di = depthInds(1):depthInds(2)
        inc = di-minWaveIndex-1;
        S = S_update.*Kzw.*exp(1j*Kzw*dz*inc);
        s = ifft(S,[],2);
        s = s(:,1:Nx).*exp(1j*Kw.*r(di,:));
        Im(di,:) = sum(s,1);
    end
    
end