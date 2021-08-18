function data = preprocess(data,Fs_old,Fs_new,matched_filter)

    data_down = zeros(size(data,1)*Fs_new/Fs_old,size(data,2),size(data,3));
    for xi = 1:size(data,2)
        for yi = 1:size(data,3)
            % Matched Filter
            data_filtered = xcorr(data(:,xi,yi),matched_filter);
            data_filtered = data_filtered(round(length(data_filtered)/2):round(length(data_filtered)/2)+size(data,1)-1);
            % Downsample
            data_down(:,xi,yi) = resample(data_filtered,Fs_new,Fs_old);
        end
    end
    
    data = squeeze(data_down);
    
end