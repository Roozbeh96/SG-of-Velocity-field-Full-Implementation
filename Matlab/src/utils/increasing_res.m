% function increasing_res(obj, L0, factor)
% 
%     L = ceil(L0*obj.delta/obj.Delta_x); %Number of intervals in one delta
%     Numrator = floor((size(obj.Gen_u_LRVF,2)-1)/L);
% 
%     obj.Kr = factor;
%     obj.Gen_u_HRVF = zeros(size(obj.Gen_u_LRVF,1), obj.Kr*(L) +1, Numrator);
%     obj.Gen_w_HRVF = zeros(size(obj.Gen_w_LRVF,1), obj.Kr*(L) +1, Numrator);
%     obj.Gen_x_HRVF = 0:obj.Delta_x/obj.Kr: (L)*obj.Delta_x;
% 
%     obj.Gen_x_LRVF = 0:obj.Delta_x:(L)*obj.Delta_x;
% 
%     progressbar('Increasing Resolution')
% 
%     k_start = 36;  
%     k_mid   = 236;   
%     k_end   = 490;    
%     min_val = 0.70;
% 
%     N = length(obj.Gen_x_HRVF);
%     dx = obj.Gen_x_HRVF(2) - obj.Gen_x_HRVF(1); 
%     k = (2*pi/(N*dx)) * [0:floor(N/2) -floor((N-1)/2):-1];
%     k_abs = abs(k);
% 
%     control_k = [0, k_start, k_mid, k_end, max(k_abs)];
%     control_v = [1, 1, min_val, 1, 1];
%     mask = interp1(control_k, control_v, k_abs, 'linear');
% 
%     for S= 1:Numrator
% 
%         ucroped = obj.Gen_u_LRVF(:, (S-1)*L+1:(S)*L+1);
%         wcroped = obj.Gen_w_LRVF(:, (S-1)*L+1:(S)*L+1);
% 
%         for r = 1:size(obj.z,2)
%             % Interpolate the row
%             % obj.Gen_u_HRVF(r, :, S) = interp1(obj.Gen_x_LRVF, ucroped(r, :),...
%             %     obj.Gen_x_HRVF, 'makima');
% 
%             u = interp1(obj.Gen_x_LRVF, ucroped(r, :),...
%                 obj.Gen_x_HRVF, 'makima');
%             U_hat = fft(u);
%             U_filtered = U_hat .* mask;
%             obj.Gen_u_HRVF(r, :, S) = real(ifft(U_filtered));
% 
%             obj.Gen_w_HRVF(r, :, S) = interp1(obj.Gen_x_LRVF, wcroped(r, :),...
%                 obj.Gen_x_HRVF, 'makima');
%         end
% 
%         progressbar((S)/Numrator)
%     end
% 
% end


function increasing_res(obj, L0, factor)
    L = ceil(L0*obj.delta/obj.Delta_x); 
    Numrator = floor((size(obj.Gen_u_LRVF,2)-1)/L);
    obj.Kr = factor;

    % Pre-calculate coordinate vectors
    Gen_x_HRVF = 0:obj.Delta_x/obj.Kr: (L)*obj.Delta_x;
    Gen_x_LRVF = 0:obj.Delta_x:(L)*obj.Delta_x;
    obj.Gen_x_HRVF = Gen_x_HRVF;
    obj.Gen_x_LRVF = Gen_x_LRVF;

    % Create temporary local variables to bypass object slicing restrictions
    % This is necessary for parfor to work efficiently
    LRVF_u = obj.Gen_u_LRVF;
    LRVF_w = obj.Gen_w_LRVF;
    num_rows = size(obj.z, 2);
    num_cols_HR = length(Gen_x_HRVF);

    % Pre-allocate temporary output arrays
    HRVF_u_temp = zeros(num_rows, num_cols_HR, Numrator);
    HRVF_w_temp = zeros(num_rows, num_cols_HR, Numrator);

    fprintf('Increasing Resolution in Parallel...\n');

    k_start = 90;%36;  
    k_mid   = 200;   
    k_end   = 418;    
    min_val = 0.80;

    N = length(obj.Gen_x_HRVF);
    dx = obj.Gen_x_HRVF(2) - obj.Gen_x_HRVF(1); 
    k = (2*pi/(N*dx)) * [0:floor(N/2) -floor((N-1)/2):-1];
    k_abs = abs(k);

    control_k = [0, k_start, k_mid, k_end, max(k_abs)];
    control_v = [1, 1, min_val, 1, 1];
    mask = interp1(control_k, control_v, k_abs, 'linear');

    % figure
    % plot(k_abs*0.055e-3,mask,'k-')
    % ylabel('$\Phi(k_{1})$','VerticalAlignment','bottom','HorizontalAlignment','center',...
    % 'Interpreter','latex');
    % xlabel('$k_{1}\eta$','Interpreter','latex');
    % set(gca,'TickLabelInterpreter','latex','XGrid','on','XMinorTick','on',...
    % 'XScale','log','YGrid','on');

    % Parallel loop over the frames
    parfor S = 1:Numrator
        % Extract local slices for this worker
        u_slice = LRVF_u(:, (S-1)*L+1:(S)*L+1);
        w_slice = LRVF_w(:, (S-1)*L+1:(S)*L+1);

        % Local storage for the current frame to avoid overhead
        frame_u = zeros(num_rows, num_cols_HR);
        frame_w = zeros(num_rows, num_cols_HR);

        for r = 1:num_rows
            % Interpolate each row
            frame_u(r, :) = interp1(Gen_x_LRVF, u_slice(r, :), Gen_x_HRVF, 'makima');
            % frame_u(r, :) = smoothdata(interp1(Gen_x_LRVF, u_slice(r, :), Gen_x_HRVF, 'makima'),'gaussian',50);

            % U_hat = fft(interp1(Gen_x_LRVF, u_slice(r, :), Gen_x_HRVF, 'makima'));
            % U_filtered = U_hat .* mask;
            % frame_u(r, :) = real(ifft(U_filtered));

            frame_w(r, :) = interp1(Gen_x_LRVF, w_slice(r, :), Gen_x_HRVF, 'makima');
            
            % W_hat = fft(interp1(Gen_x_LRVF, w_slice(r, :), Gen_x_HRVF, 'makima'));
            % W_filtered = W_hat .* mask;
            % frame_w(r, :) = real(ifft(W_filtered));
        end

        % Assign the processed frame back to the temp 3D array
        HRVF_u_temp(:, :, S) = frame_u;
        HRVF_w_temp(:, :, S) = frame_w;
    end

    % Assign back to the object
    obj.Gen_u_HRVF = HRVF_u_temp;
    obj.Gen_w_HRVF = HRVF_w_temp;
    obj.Gen_u_HRVF(:,:,1:149) = [];
    obj.Gen_w_HRVF(:,:,1:149) = [];

    fprintf('Done.\n');
end