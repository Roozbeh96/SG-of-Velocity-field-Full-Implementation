% function  cdf_omega(obj)
%     % Find the distribution of mean omega in Shear Layers for the
%     % matter of matching decile process.
%     Lambda_ci_filtered = obj.Lambda_ci;
%     omega_filtered = obj.omega;
% 
%     for r = 1:size(obj.Lambda_ci, 1)
%             Lambda_ci_filtered(r, :, :) = obj.Lambda_ci(r, :, :) .* (abs(obj.Lambda_ci(r, :, :)) >= 0.35*obj.Lambda_cirms(r));
%             omega_filtered(r,:,:) = obj.omega(r, :, :) .* (abs(obj.Lambda_ci(r, :, :)) >= 0.35*obj.Lambda_cirms(r));
%     end
% 
%     ind = 1;
%     mean_omega_y_norm =zeros();
% 
%     for S = 1: size(omega_filtered,3)
%         Matrix = omega_filtered(:,:,S);
% 
%         binaryMatrix = Matrix ~= 0;
% 
%         [labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);
%         for j = 1: numComponents
%             [rows, cols] = find(labeledMatrix == j);
%             values = Matrix(labeledMatrix == j);
% 
%             row_centroid = sum(rows .* values) / sum(values);
%             col_centroid = sum(cols .* values) / sum(values);
% 
%             row_centroid = round(row_centroid);
%             col_centroid = round(col_centroid);
%             if row_centroid <= 0 || row_centroid >= size(omega_filtered,1)
%                 continue
%             end
%             if col_centroid <= 0 || col_centroid >= size(omega_filtered,2)
%                 continue
%             end
%             if Matrix(row_centroid,col_centroid) == 0
%                 continue
%             end
%             mean_omega_y_norm(1,ind) = mean(values)/obj.omega_rms(row_centroid,1);
% 
%             ind = ind+1;
%         end
%     end
% 
%     [cdf_, mean_omega_] = ecdf(abs(mean_omega_y_norm));
%     obj.cdf_omega  = cdf_;
%     obj.mean_omega = mean_omega_;
% 
% end


function cdf_omega(obj)
    tic
    % Find the distribution of mean omega in Shear Layers
    % 1. Vectorized Filtering (keep strong shear layers)
    fprintf('Finding percentile of the omega in the field in parallel...\n');
    gamma_filter = 0.35;
    threshold = gamma_filter * obj.Lambda_cirms; % Assuming this is a column vector per row
    mask = abs(obj.Lambda_ci) >= threshold;
    omega_filtered = obj.omega .* mask;

    % Prepare for parallel processing
    num_frames = size(omega_filtered, 3);
    omega_rms_local = obj.omega_rms;
    sz1 = size(omega_filtered, 1);
    sz2 = size(omega_filtered, 2);

    % Temporary cell array to store results from each worker
    temp_results = cell(num_frames, 1);

    % 2. Parallel Loop over frames
    parfor S = 1:num_frames
        Matrix = omega_filtered(:, :, S);
        binaryMatrix = Matrix ~= 0;
        [labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);

        % Local storage for this specific frame
        frame_values = []; 

        for j = 1:numComponents
            % Logical indexing is faster than 'find' for gathering values
            component_mask = (labeledMatrix == j);
            values = Matrix(component_mask);
            [rows, cols] = find(component_mask);

            sum_val = sum(values);
            if sum_val == 0, continue; end % Avoid division by zero

            row_centroid = round(sum(rows .* values) / sum_val);
            col_centroid = round(sum(cols .* values) / sum_val);

            % Bounds and validity check
            if row_centroid <= 0 || row_centroid > sz1 || ...
               col_centroid <= 0 || col_centroid > sz2 || ...
               Matrix(row_centroid, col_centroid) == 0
                continue
            end

            % Normalize and store locally
            val_norm = mean(values) / omega_rms_local(row_centroid, 1);
            frame_values = [frame_values, val_norm];
        end
        temp_results{S} = frame_values;
    end

    % 3. Combine results and compute CDF
    mean_omega_y_norm = [temp_results{:}];

    [cdf_, mean_omega_] = ecdf(abs(mean_omega_y_norm));
    obj.cdf_omega  = cdf_;
    obj.mean_omega = mean_omega_;
    toc
end