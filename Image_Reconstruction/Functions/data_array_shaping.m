function outdata = data_array_shaping(data,rec_spacing,domain,dr)

    outdata = zeros(size(data,1),size(domain,2),size(domain,3));

    X = linspace(-size(domain,2)*dr/2,size(domain,2)*dr/2,size(domain,2));
    Y = linspace(-size(domain,3)*dr/2,size(domain,3)*dr/2,size(domain,3));

    X_rec = linspace(-size(data,2)*rec_spacing/2,size(data,2)*rec_spacing/2,size(data,2));
    Y_rec = linspace(-size(data,3)*rec_spacing/2,size(data,3)*rec_spacing/2,size(data,3));

    for xi = 1:length(X_rec)
        for yi = 1:length(Y_rec)
            [~,indX] = min(abs(X-X_rec(xi)));
            [~,indY] = min(abs(Y-Y_rec(yi)));
            outdata(:,indX,indY) = data(:,xi,yi);
        end
    end
        
end