function absorption_coefficient = calculate_absorption(wavelength,wavenumber)
    
    % Read CSV file containing absorption coefficients vs. wavelength
    M = csvread('PureWaterRefractiveIndex.csv',1,0);
    M = M(1:115,:);
    lambdas = M(:,1)*1e-6;
    absorptions = 4*pi*M(:,2)./lambdas;
    % Interpolate the data to be valid for nanometer increments
    lambdas_interp = 0.23e-6:0.001e-6:11e-6;
    absorptions_interp = interp1(lambdas,absorptions,lambdas_interp);

    if strcmp(wavelength,'Optimum')
        absorption_coefficient = wavenumber;
        index = [];
        epsilon = 0.1;
        while isempty(index)
            index = find(absorptions_interp>wavenumber-epsilon & absorptions_interp<wavenumber+epsilon);
            epsilon = epsilon+0.1;
        end
            
        opt_wavelength = sprintf('Optimal Laser Wavelength: %d nm \n',round(lambdas_interp(index)*1e9));
        fprintf(opt_wavelength)   
    else 
        index = [];
        epsilon = 0.1e-9;
        while isempty(index)
            index = find(lambdas_interp>wavelength*1e-9-epsilon & lambdas_interp<wavelength*1e-9+epsilon);
            epsilon = epsilon+0.1e-9;
        end
        
        absorption_coefficient = absorptions_interp(index);  
    end
    
end