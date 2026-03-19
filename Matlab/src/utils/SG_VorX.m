function SG_VorX(obj, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
    u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,...
    omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,...
    omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles)
    %#ok<*AGROW>
        %{
            The Pearson correlation coefficient between radius and the
            maximum azimuthal velocity of the generated vortices at
            different wall normal elevations.
            Near : z/k_s < 2
            Far : z/k_s >= 2
        %}
        % Find the omega of the generated vortices
        % omega_near = u_wOu_tau_near./r_wOlambda_T_near*obj.u_tau/obj.lambda;
        % omega_near_norm= omega_near/std(omega_near,0,1);
        % 
        % omega_far = u_wOu_tau_far./r_wOlambda_T_far*obj.u_tau/obj.lambda;
        % omega_far_norm= omega_far/std(omega_far,0,1);
        % 
        % % Dividing the generated vortices into deciles for the matter of
        % % matching
        % 
        % omega_near_deciles = cell(1, 10); % Create a container for 10 groups
        % u_wOu_tau_near_deciles = cell(1, 10);
        % r_wOlambda_T_near_deciles = cell(1, 10);
        % 
        % for i = 1:10
        %     L = (i-1) * 10;
        %     U = i * 10;
        % 
        %     valL = prctile(omega_near_norm, L);
        %     valU = prctile(omega_near_norm, U);
        % 
        %     % Store the samples for this specific bracket
        %     omega_near_deciles{i} = omega_near_norm(omega_near_norm >= valL & omega_near_norm <= valU);
        %     u_wOu_tau_near_deciles{i} = u_wOu_tau_near(omega_near_norm >= valL & omega_near_norm <= valU);
        %     r_wOlambda_T_near_deciles{i} = r_wOlambda_T_near(omega_near_norm >= valL & omega_near_norm <= valU);
        % end
        % 
        % omega_far_deciles = cell(1, 10); % Create a container for 10 groups
        % u_wOu_tau_far_deciles = cell(1, 10);
        % r_wOlambda_T_far_deciles = cell(1, 10);
        % for i = 1:10
        %     L = (i-1) * 10;
        %     U = i * 10;
        % 
        %     valL = prctile(omega_far_norm, L);
        %     valU = prctile(omega_far_norm, U);
        % 
        %     % Store the samples for this specific bracket
        %     omega_far_deciles{i} = omega_far_norm(omega_far_norm >= valL & omega_far_norm <= valU);
        %     u_wOu_tau_far_deciles{i} = u_wOu_tau_far(omega_far_norm >= valL & omega_far_norm <= valU);
        %     r_wOlambda_T_far_deciles{i} = r_wOlambda_T_far(omega_far_norm >= valL & omega_far_norm <= valU);
        % end

        % some counter parameters
        ii = 1; % counter for r_omega and u_omega of z/ks<2 vortices
        %(could be near-wall or lambda_ci vortices)
        iii = 1;% counter for z/ks>=2 vortices (lambda_ci vortices)
        progressbar('Adding VorX to the Velocity Field')
        % Filter low intense lambda_ci regions.
        Lambda_ci_filtered = obj.Lambda_ci;
        omega_filtered = obj.omega;

        for r = 1:size(obj.Lambda_ci, 1)
                Lambda_ci_filtered(r, :, :) = obj.Lambda_ci(r, :, :) .* (abs(obj.Lambda_ci(r, :, :)) >= 0.35*obj.Lambda_cirms(r));
                omega_filtered(r,:,:) = obj.omega(r, :, :) .* (abs(obj.Lambda_ci(r, :, :)) >= 0.35*obj.Lambda_cirms(r));
        end

        % Find the distribution of mean omega in Shear Layers for the
        % matter of matching decile process.
        % ind = 1;
        % mean_omega_y_norm =zeros();
        % 
        % for S = 1: size(omega_filtered,3)
        %     Matrix = omega_filtered(:,:,S);
        % 
        %     binaryMatrix = Matrix ~= 0;
        % 
        %     [labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);
        %     for j = 1: numComponents
        %         [rows, cols] = find(labeledMatrix == j);
        %         values = Matrix(labeledMatrix == j);
        % 
        %         row_centroid = sum(rows .* values) / sum(values);
        %         col_centroid = sum(cols .* values) / sum(values);
        % 
        %         row_centroid = round(row_centroid);
        %         col_centroid = round(col_centroid);
        %         if row_centroid <= 0 || row_centroid >= size(omega_filtered,1)
        %             continue
        %         end
        %         if col_centroid <= 0 || col_centroid >= size(omega_filtered,2)
        %             continue
        %         end
        %         if Matrix(row_centroid,col_centroid) == 0
        %             continue
        %         end
        %         mean_omega_y_norm(1,ind) = mean(values)/obj.omega_rms(row_centroid,1);
        % 
        %         ind = ind+1;
        %     end
        % end
        % 
        % [cdf_, mean_omega_] = ecdf(abs(mean_omega_y_norm));

        cdf_omega_ = obj.cdf_omega;
        mean_omega_ = obj.mean_omega;


        uprime = obj.Gen_u_HRVF - mean(obj.Gen_u_HRVF,3); % fluctuating velocity
        % uprime = obj.Gen_u_HRVF - 1/0.39*log(obj.z'./obj.z_0); % fluctuating velocity

        for S = 1:size(obj.Gen_u_HRVF,3) %S is iterating over the Velocity field with a length of \delta
            x_r_margin_near_wall = 0;
            temp_near_wall_VorX_counter = 1;
            %{
                Adding near-wall vortices(square shape). As long as we are
                in the velocity field region, we keep add near-wall
                vortices
            %}
            while x_r_margin_near_wall < max(obj.Gen_x_HRVF)

                N = 1;
                r_omega = r_wOlambda_T_near(ii)*obj.lambda;
                u_omega = u_wOu_tau_near(ii)*obj.u_tau;
                TVorX = [1, rho_uw_VorX_near(ii); 0, sqrt(1 - rho_uw_VorX_near(ii)^2)];
                Gama=2*pi*u_omega*r_omega/(1-exp(-1)); 
                ii = ii+1;
                % Find center position of the generated vortex.
                zc_near_wall = N*r_omega + obj.z(1);
                xc_near_wall = N*r_omega + x_r_margin_near_wall;
                x_r_margin_near_wall = xc_near_wall + N*r_omega;
                [indzmid] = find(obj.z>=zc_near_wall,1);
                [indxmid] = find(obj.Gen_x_HRVF>=xc_near_wall,1);

                %{ 
                    Two conditions for putting near-wall vortices:
                    1. uprime at center of the vortex should be positive.
                    2. randomly generated binary number should be one. 
                %}
                if uprime(indzmid,indxmid,S)<0
                    continue
                end
                if randi([0, 1])==0
                    continue
                end
                % Find distance between the near-wall vortices. 
                if temp_near_wall_VorX_counter > 1
                    near_wall_VorX_data{S}{temp_near_wall_VorX_counter} = ...
                        dictionary(["xc_VorX[m]","zc_VorX[m]","u_omega[m/s]","r_omega[m]",...
                        "Dist_to_prev_VorX[m]"],[xc_near_wall, zc_near_wall,u_omega,...
                        r_omega, xc_near_wall-prev_pos_xc_near_wall_VorX]); 
                    prev_pos_xc_near_wall_VorX = xc_near_wall;
                    temp_near_wall_VorX_counter = temp_near_wall_VorX_counter +1;

                else 
                    near_wall_VorX_data{S}{temp_near_wall_VorX_counter} = ...
                        dictionary(["xc_VorX","zc_VorX","u_omega","r_omega",...
                        "Dist_to_prev_VorX"],[xc_near_wall, zc_near_wall,u_omega,...
                        r_omega,0]); 
                    prev_pos_xc_near_wall_VorX = xc_near_wall;
                    temp_near_wall_VorX_counter = temp_near_wall_VorX_counter +1;
                end


                % Generating vortex field. Explained in methodology of the paper. 
                [indz] = find(obj.z>=zc_near_wall-N*r_omega &...
                    obj.z<=zc_near_wall+N*r_omega);
                [indx] = find(obj.Gen_x_HRVF>=xc_near_wall-N*r_omega &...
                    obj.Gen_x_HRVF<=xc_near_wall+N*r_omega);
                [X,Z] = meshgrid(obj.Gen_x_HRVF(indx),obj.z(indz));
                % square shape vortex
                r = sqrt((X-xc_near_wall).^2+(Z-zc_near_wall).^2);

                % Generating the azimuthal velocity using equation 3.1.
                uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
                nan_logical_array = isnan(uazi);
                [nan_indices] = find(nan_logical_array);
                uazi(nan_indices) = 0;

                Theta = atan2((Z-zc_near_wall),(X-xc_near_wall));
                nan_logical_array = isnan(Theta);
                [nan_indices] = find(nan_logical_array);
                Theta(nan_indices) = 0;

                u_VorX = uazi.*sin(Theta);
                %Check whether generated vortex can be resolved in our mesh
                %domain or not.
                if isempty(u_VorX)
                    continue
                end
                % Subtracting the center velocity from u_VorX
                u_VorX(uazi~=0)=u_VorX(uazi~=0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                w_VorX = -uazi.*cos(Theta);

                A = [reshape(u_VorX,[],1) reshape(w_VorX,[],1)];
                B = A*TVorX;
                u_TVorX = reshape(B(:,1), size(u_VorX));
                w_TVorX = reshape(B(:,2), size(w_VorX));


                % Adding to the velocity field.
                obj.Gen_u_HRVF(indz,indx,S) = obj.Gen_u_HRVF(indz,indx,S)+ u_TVorX;
                obj.Gen_w_HRVF(indz,indx,S) = obj.Gen_w_HRVF(indz,indx,S)+ w_TVorX;


            end

            % Lambda_ci vortices (Primary vortices)
            Matrix_lambdaci = Lambda_ci_filtered(:,:,S);
            Matrix_omega = omega_filtered(:,:,S);
            binaryMatrix = Matrix_lambdaci ~= 0;

            % Label connected components
            [labeledMatrix, numb_clusters] = bwlabel(binaryMatrix, 4);  % 4-connectivity

            % Initialize variables to store centroids
            centroids = zeros(1,6);
            ix = 1;
            % Calculate the mass center (centroid) for each connected component
            for x = 1:numb_clusters
                % Extract the indices of the current component
                [rows, cols] = find(labeledMatrix == x);
                values = Matrix_lambdaci(labeledMatrix == x);
                sign_region_mean = sign(mean(values));
                omega_region_mean = mean(Matrix_omega(labeledMatrix == x));
                % Calculate the centroid
                row_centroid = sum(rows .* values) / sum(values);
                col_centroid = sum(cols .* values) / sum(values);

                % Round the centroids to the nearest integer
                row_centroid = round(row_centroid);
                col_centroid = round(col_centroid);
                if row_centroid <= 0 || row_centroid > size(obj.z,2)
                    continue
                end
                if col_centroid <= 0 || col_centroid > size(obj.Gen_x_HRVF,2)
                    continue
                end
                % Store the centroid
                centroids(ix, :) = [row_centroid, col_centroid, sign_region_mean,...
                    min(cols), max(cols), omega_region_mean/obj.omega_rms(row_centroid,1)];
                ix = ix + 1;
            end

            zcs = obj.z(centroids(:, 1));
            xcs = obj.Gen_x_HRVF(centroids(:, 2));
            rot = centroids(:, 3);
            omega_ = centroids(:, 6);
            Numb_vortices  = size(centroids,1);
            pro_counter = 1; %counter for lambda_ci prograde vortices
            retro_counter = 1; %counter for lambda_ci retrograde vortices
            N = 2;
            for i=1:Numb_vortices               
                if zcs(i)< 2*obj.ks % close to the wall
                    %Primary vortex
                    ind_Pri_VorX = ii;
                    [pro_counter, retro_counter,r_wOlambda_T_near,u_wOu_tau_near] = ...
                    Lambdaci_VorX(obj,i,S,xcs,zcs,rot,r_wOlambda_T_near,...
                    u_wOu_tau_near, rho_uw_VorX_near(ii),ii, pro_counter, retro_counter, N,...
                    cdf_omega_, mean_omega_,omega_(i),omega_near_deciles, r_wOlambda_T_near_deciles, u_wOu_tau_near_deciles);
                    x_ = xcs(i)+N*r_wOlambda_T_near(ind_Pri_VorX)*obj.lambda;
                    ii = ii+1;
                    %Filling the right side till reaching to right margin(Scondary vortices)
                    while  x_+N*r_wOlambda_T_near(ii)*obj.lambda<...
                            obj.Gen_x_HRVF(centroids(i, 5))
                        [pro_counter, retro_counter,r_wOlambda_T_near,u_wOu_tau_near] = ...
                            Lambdaci_VorX(obj,i,S,xcs,zcs,rot,r_wOlambda_T_near,...
                            u_wOu_tau_near, rho_uw_VorX_near(ii), ii, pro_counter, retro_counter, N,...
                            cdf_omega_, mean_omega_,omega_(i),omega_near_deciles, r_wOlambda_T_near_deciles, u_wOu_tau_near_deciles);
                        x_ = x_+2*N*r_wOlambda_T_near(ii)*obj.lambda;
                        ii = ii+1;
                    end
                    x_ = xcs(i)-N*r_wOlambda_T_near(ind_Pri_VorX)*obj.lambda;
                    ii = ii +1;
                    %Filling the left side till reaching to left margin(Scondary vortices)
                    while x_-N*r_wOlambda_T_near(ii)*obj.lambda>...
                            obj.Gen_x_HRVF(centroids(i, 4))
                        [pro_counter, retro_counter, r_wOlambda_T_near,u_wOu_tau_near] = ...
                            Lambdaci_VorX(obj,i,S,xcs,zcs,rot,r_wOlambda_T_near,...
                            u_wOu_tau_near, rho_uw_VorX_near(ii), ii, pro_counter, retro_counter, N,...
                            cdf_omega_, mean_omega_,omega_(i),omega_near_deciles, r_wOlambda_T_near_deciles, u_wOu_tau_near_deciles);
                        x_ = x_-2*N*r_wOlambda_T_near(ii)*obj.lambda;
                        ii = ii+1;
                    end
                    ii = ii +1;
                else %Far from the wall(z/ks>=2)
                    %Primary vortex
                    ind_Pri_VorX = iii;
                    [pro_counter, retro_counter, r_wOlambda_T_far, u_wOu_tau_far] = ...
                    Lambdaci_VorX(obj,i,S,xcs,zcs,rot,r_wOlambda_T_far,...
                    u_wOu_tau_far, rho_uw_VorX_far(iii), iii, pro_counter, retro_counter, N,...
                    cdf_omega_, mean_omega_,omega_(i),omega_far_deciles, r_wOlambda_T_far_deciles, u_wOu_tau_far_deciles);
                    x_ = xcs(i)+N*r_wOlambda_T_far(ind_Pri_VorX)*obj.lambda;
                    iii = iii + 1;
                    %Filling the right side till reaching to right margin(Scondary vortices)
                    while  x_+N*r_wOlambda_T_far(iii)*obj.lambda<...
                            obj.Gen_x_HRVF(centroids(i, 5))
                        [pro_counter, retro_counter, r_wOlambda_T_far, u_wOu_tau_far] = ...
                            Lambdaci_VorX(obj,i,S,xcs,zcs,rot,r_wOlambda_T_far,...
                            u_wOu_tau_far, rho_uw_VorX_far(iii), iii, pro_counter, retro_counter, N,...
                            cdf_omega_, mean_omega_,omega_(i),omega_far_deciles, r_wOlambda_T_far_deciles, u_wOu_tau_far_deciles);
                        x_ = x_+2*N*r_wOlambda_T_far(iii)*obj.lambda;
                        iii = iii + 1;
                    end
                    x_ = xcs(i)-N*r_wOlambda_T_far(ind_Pri_VorX)*obj.lambda;
                    iii = iii +1;
                    %Filling the left side till reaching to left margin(Scondary vortices)
                    while x_-N*r_wOlambda_T_far(iii)*obj.lambda>...
                            obj.Gen_x_HRVF(centroids(i, 4))
                        [pro_counter, retro_counter, r_wOlambda_T_far, u_wOu_tau_far] = ...
                            Lambdaci_VorX(obj,i,S,xcs,zcs,rot,r_wOlambda_T_far,...
                            u_wOu_tau_far, rho_uw_VorX_far(iii), iii, pro_counter, retro_counter, N,...
                            cdf_omega_, mean_omega_,omega_(i),omega_far_deciles, r_wOlambda_T_far_deciles, u_wOu_tau_far_deciles);
                        x_ = x_-2*N*r_wOlambda_T_far(iii)*obj.lambda;
                        iii = iii + 1;
                    end
                    iii = iii +1;
                end
            end
            progressbar((S)/size(obj.Gen_u_HRVF,3))
        end

end
