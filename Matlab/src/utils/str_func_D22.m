function str_func_D22(obj)
    %{
    More info about how to compute the structure function in:
    https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/
    article/stochastic-modal-velocity-field-in-roughwall-turbulence/
    800D97E60609F4FE4D92847CCDB1C7A0
    %}
    wprime = obj.Gen_w_HRVF - mean(obj.Gen_w_HRVF,3);
    concat_wprime = reshape(wprime,size(wprime, 1), []);
    dx = obj.Gen_x_HRVF(2) - obj.Gen_x_HRVF(1);
    win15del = ceil(15*obj.delta/dx);
    Numb_frame = floor(size(concat_wprime,2)/win15del);
    
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
    obj.D22 = zeros(size(concat_wprime,1), num_r, Numb_frame);
    obj.r11 = dx * r_indices; % Vectorized calculation of physical distances
    
    progressbar('Computing Structure Function D22')
    
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
                
                obj.D22(elev, ri, S) = mean((concat_wprime(elev, frame_start : frame_end-r) - ...
                    concat_wprime(elev, frame_start+r : frame_end)).^2, 2);
            end
        end
        
        progressbar(S/Numb_frame)
    end
end