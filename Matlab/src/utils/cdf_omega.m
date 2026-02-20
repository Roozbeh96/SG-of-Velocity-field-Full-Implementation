function  cdf_omega(obj)
    % Find the distribution of mean omega in Shear Layers for the
    % matter of matching decile process.
    Lambda_ci_filtered = obj.Lambda_ci;
    omega_filtered = obj.omega;
    
    for r = 1:size(obj.Lambda_ci, 1)
            Lambda_ci_filtered(r, :, :) = obj.Lambda_ci(r, :, :) .* (abs(obj.Lambda_ci(r, :, :)) >= 0.35*obj.Lambda_cirms(r));
            omega_filtered(r,:,:) = obj.omega(r, :, :) .* (abs(obj.Lambda_ci(r, :, :)) >= 0.35*obj.Lambda_cirms(r));
    end
    
    ind = 1;
    mean_omega_y_norm =zeros();
    
    for S = 1: size(omega_filtered,3)
        Matrix = omega_filtered(:,:,S);
        
        binaryMatrix = Matrix ~= 0;
        
        [labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);
        for j = 1: numComponents
            [rows, cols] = find(labeledMatrix == j);
            values = Matrix(labeledMatrix == j);
    
            row_centroid = sum(rows .* values) / sum(values);
            col_centroid = sum(cols .* values) / sum(values);
    
            row_centroid = round(row_centroid);
            col_centroid = round(col_centroid);
            if row_centroid <= 0 || row_centroid >= size(omega_filtered,1)
                continue
            end
            if col_centroid <= 0 || col_centroid >= size(omega_filtered,2)
                continue
            end
            if Matrix(row_centroid,col_centroid) == 0
                continue
            end
            mean_omega_y_norm(1,ind) = mean(values)/obj.omega_rms(row_centroid,1);
    
            ind = ind+1;
        end
    end
    
    [cdf_, mean_omega_] = ecdf(abs(mean_omega_y_norm));
    obj.cdf_omega  = cdf_;
    obj.mean_omega = mean_omega_;

end