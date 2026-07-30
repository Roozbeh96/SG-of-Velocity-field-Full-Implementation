% function Lambda_ci(obj)
% 
%     obj.Lambda_ci = zeros(size(obj.Gen_u_HRVF));
%     DelxGen = obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1);
%     DelzGen = obj.z(2)-obj.z(1);
%     progressbar('Computing Lambda_{ci} Field')
%     for S= 1:size(obj.Lambda_ci,3)
% 
%         [dudx, dudz] = gradient(obj.Gen_u_HRVF(:,:, S),DelxGen,DelzGen);
%         [dwdx, dwdz] = gradient(obj.Gen_w_HRVF(:,:, S),DelxGen,DelzGen);
% 
%         omega = dwdx - dudz;
%             for r = 1:size(obj.Lambda_ci,1)
% 
%                 for c = 1:size(obj.Lambda_ci,2)
% 
%                     temp = [dudx(r,c) dudz(r,c); dwdx(r,c) dwdz(r,c)];
%                     obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
%                         *sign(omega(r,c));
% 
%                     if  unique(abs(imag(eig(temp)))) ~= 0
%                         obj.omega(r,c,S) = omega(r,c);
%                     end
% 
%                 end
%             end
%         progressbar((S)/size(obj.Lambda_ci,3))
% 
%     end
%     Lambda_cirms_ = zeros(size(obj.Lambda_ci,1),1);
%     omega_rms_ = zeros(size(obj.Lambda_ci,1),1);
% 
%     for z_ = 1:size(obj.Lambda_ci,1)
%         temp_ = obj.Lambda_ci(z_,:,:);
%         temp_ = temp_(temp_~=0);
%         Lambda_cirms_(z_,1) = std(temp_,0,1);
% 
%         temp__ = obj.omega(z_,:,:);
%         temp__ = temp__(temp__~=0);
%         omega_rms_(z_,1)= std(temp__,0,1);
% 
%     end
% 
% 
%     obj.Lambda_cirms = Lambda_cirms_;
%     obj.omega_rms = omega_rms_;
% end

function Lambda_ci(obj)
    tic
    % 1. Extract data to local variables (Required for parfor)
    Gen_u_HRVF = obj.Gen_u_HRVF;
    Gen_w_HRVF = obj.Gen_w_HRVF;
    Gen_x_HRVF = obj.Gen_x_HRVF;
    z = obj.z;

    [nRows, nCols, nFrames] = size(Gen_u_HRVF);

    % Preallocate local arrays
    Lambda_ci_local = zeros(nRows, nCols, nFrames);
    omega_local_out = zeros(nRows, nCols, nFrames);

    DelxGen = Gen_x_HRVF(2) - Gen_x_HRVF(1);
    DelzGen = z(2) - z(1);

    fprintf('Computing Lambda_{ci} Field in parallel...\n');

    % 2. Parallel loop over frames
    parfor S = 1:nFrames
        % Compute gradients for the current frame
        [dudx, dudz] = gradient(Gen_u_HRVF(:,:,S), DelxGen, DelzGen);
        [dwdx, dwdz] = gradient(Gen_w_HRVF(:,:,S), DelxGen, DelzGen);

        omega = dwdx - dudz;

        % Vectorized eigenvalue calculation for 2x2 matrix
        % The imaginary part of eigenvalues for [a b; c d] is 
        % sqrt(abs(min(0, (a+d)^2 - 4*(a*d - b*c)))) / 2
        tr = dudx + dwdz;
        detJ = (dudx .* dwdz) - (dudz .* dwdx);
        discriminant = (tr.^2) - 4.*detJ;

        % Calculate swirling strength
        swirl = zeros(nRows, nCols);
        mask = discriminant < 0;
        swirl(mask) = sqrt(abs(discriminant(mask))) / 2;

        % Assign to local frame results
        Lambda_ci_local(:,:,S) = swirl .* sign(omega);

        % Store omega where swirling exists
        temp_omega = zeros(nRows, nCols);
        temp_omega(mask) = omega(mask);
        omega_local_out(:,:,S) = temp_omega;
    end

    % 3. Statistics Calculation (RMS)
    Lambda_cirms_ = zeros(nRows, 1);
    omega_rms_ = zeros(nRows, 1);

    for z_ = 1:nRows
        temp_ = Lambda_ci_local(z_,:,:);
        temp_ = temp_(temp_ ~= 0);
        if ~isempty(temp_)
            Lambda_cirms_(z_,1) = std(temp_, 0, 1);
        end

        temp__ = omega_local_out(z_,:,:);
        temp__ = temp__(temp__ ~= 0);
        if ~isempty(temp__)
            omega_rms_(z_,1) = std(temp__, 0, 1);
        end
    end

    % 4. Assign results back to the object
    obj.Lambda_ci = Lambda_ci_local;
    obj.omega = omega_local_out;
    obj.Lambda_cirms = Lambda_cirms_;
    obj.omega_rms = omega_rms_;
    toc
end