function str_func(obj)
    %{
    More info about how to compute the structure function in:
    https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/
    article/stochastic-modal-velocity-field-in-roughwall-turbulence/
    800D97E60609F4FE4D92847CCDB1C7A0
    %}
    uprime = obj.Gen_u_HRVF - mean(obj.Gen_u_HRVF,3);
    concat_uprime = reshape(uprime,size(uprime, 1), []);
    dx = obj.Gen_x_HRVF(2) - obj.Gen_x_HRVF(1);
    win15del = ceil(15*obj.delta/dx);
    Numb_frame = floor(size(concat_uprime,2)/win15del);
    
    % --- LOGARITHMIC SKIP SETUP ---
    % Define how many points you want in your log scale (e.g., 50 or 100 points)
    num_log_points = 50; 
    
    % Generate unique, rounded integer steps logarithmically from 0 to win12del
    % We do logspace from 0 to log10(win12del), round it, and add 0.
    r_indices = unique([0, round(logspace(0, log10(win15del), num_log_points))]);
    % Ensure it doesn't accidentally exceed bounds
    r_indices(r_indices > win15del) = []; 
    num_r = length(r_indices);
    % ------------------------------

    % Preallocate only for the skipped points to save memory and speed up max()
    obj.D11 = zeros(size(concat_uprime,1), num_r, Numb_frame);
    obj.r11 = dx * r_indices; % Vectorized calculation of physical distances
    
    obj.epsilon_str = zeros(size(concat_uprime,1),Numb_frame);
    obj.eta_str = zeros(size(concat_uprime,1),Numb_frame);
    
    progressbar('Computing Structure Function')
    
    [~, ind_i]   = find(obj.z/obj.delta >= 0.015, 1, 'first');
    [~, ind_ii]  = find(obj.z/obj.delta >= 0.05,   1, 'first');
    [~, ind_iii] = find(obj.z/obj.delta >= 0.0919, 1, 'first');
    [~, ind_iv]  = find(obj.z/obj.delta >= 0.14,   1, 'first');
    [~, ind_v]   = find(obj.z/obj.delta >= 0.24,   1, 'first');

    Z = [ind_i, ind_ii, ind_iii, ind_iv, ind_v];
    
    for S = 1:Numb_frame
        frame_start = (S-1)*win15del + 1;
        frame_end   = S*win15del + 1;
        
        for elev = Z
            % Loop through our log-spaced indices instead of every integer
            for ri = 1:num_r
                r = r_indices(ri); 
                
                obj.D11(elev, ri, S) = mean((concat_uprime(elev, frame_start : frame_end-r) - ...
                    concat_uprime(elev, frame_start+r : frame_end)).^2, 2);
            end
        end
        
        % Avoid division by zero at r=0 for the epsilon calculation
        % We ignore the first point (r=0) by indexing from 2 onwards
        [max_val] = max(obj.D11(:, 2:end, S) .* (obj.r11(2:end).^(-2/3)), [], 2);
        
        obj.epsilon_str(:,S) = (2./max_val).^(-3/2);
        obj.eta_str(:,S) = (obj.nu^3 ./ obj.epsilon_str(:,S)).^0.25;
        progressbar(S/Numb_frame)
    end
end



% function str_func(obj)
%     % --- Pre-processing (Same as before) ---
%     uprime = obj.Gen_u_HRVF - mean(obj.Gen_u_HRVF, 3);
%     concat_uprime = reshape(uprime, size(uprime, 1), []);
%     dx = obj.Gen_x_HRVF(2) - obj.Gen_x_HRVF(1);
%     win12del = ceil(12*obj.delta/dx);
%     Numb_frame = floor(size(concat_uprime, 2) / win12del);
% 
%     % --- Determine Elevations ---
%     [~, ind_i]   = find(obj.z/obj.delta >= 0.015, 1, 'first');
%     [~, ind_ii]  = find(obj.z/obj.delta >= 0.05,   1, 'first');
%     [~, ind_iii] = find(obj.z/obj.delta >= 0.0919, 1, 'first');
%     [~, ind_iv]  = find(obj.z/obj.delta >= 0.14,   1, 'first');
%     [~, ind_v]   = find(obj.z/obj.delta >= 0.24,   1, 'first');
% 
%     Z = [ind_i, ind_ii, ind_iii, ind_iv, ind_v];
% 
%     % --- Local Variable Assignment (CRITICAL for parfor) ---
%     % We create local arrays so parfor can "slice" them correctly.
%     r11_local = (0:win12del) * dx;
%     D11_local = zeros(size(concat_uprime, 1), win12del + 1, Numb_frame);
%     eps_local = zeros(size(concat_uprime, 1), Numb_frame);
%     eta_local = zeros(size(concat_uprime, 1), Numb_frame);
%     nu_val    = obj.nu; % Local copy of scalar property
% 
%     % Setup for progress monitoring (standard progressbar won't work inside parfor)
%     fprintf('Computing Structure Function in Parallel...\n');
% 
%     % --- Parallel Loop ---
%     parfor S = 1:Numb_frame
%         % Temporary storage for the current frame's D11 results
%         % This minimizes communication between workers
%         frame_D11 = zeros(length(Z), win12del + 1);
% 
%         for i = 1:length(Z)
%             elev = Z(i);
%             % Pre-calculate frame indices to keep the inner loop clean
%             start_idx = (S-1)*win12del + 1;
%             end_idx   = S*win12del + 1;
% 
%             for r = 0:win12del
%                 % Compute structure function for specific r and elevation
%                 u_segment_1 = concat_uprime(elev, start_idx : end_idx - r);
%                 u_segment_2 = concat_uprime(elev, start_idx + r : end_idx);
% 
%                 frame_D11(i, r+1) = mean((u_segment_1 - u_segment_2).^2, 2);
%             end
%         end
% 
%         % Assign frame results into the sliced local array
%         % We use indexing to help MATLAB understand how to distribute the data
%         D11_local(Z, :, S) = frame_D11;
% 
%         % Compute epsilon and eta for this frame
%         % Note: We use the local r11 and the frame_D11 we just computed
%         [max_val] = max(frame_D11 .* (r11_local.^(-2/3)), [], 2);
% 
%         % Map frame results back to the sliced output arrays
%         eps_local(Z, S) = (2 ./ max_val).^(-3/2);
%         eta_local(Z, S) = (nu_val^3 ./ eps_local(Z, S)).^0.25;
%     end
% 
%     % --- Post-processing: Assign back to object ---
%     obj.r11 = r11_local;
%     obj.D11 = D11_local;
%     obj.epsilon_str = eps_local;
%     obj.eta_str = eta_local;
% 
%     fprintf('Computation Complete.\n');
% end