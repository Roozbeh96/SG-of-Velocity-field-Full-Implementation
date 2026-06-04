% function [Gen_u_prof, Gen_w_prof, log_data] = SG_VelProf(obj)
% %SGVF generate step like velocity profiles. Detail of line will explain in
% %the folllowing lines.
% %#ok<*AGROW>
%     N_rand = 1e6;
%     repro_numb = 137; %for reproducibility, we generate one set of random numbers
%     [rand_u,rand_w] = genGaussCop(obj.ro_uw, N_rand, repro_numb);
%     Gen_u_prof = zeros(size(obj.z,2),obj.N_prof);
%     Gen_w_prof = zeros(size(obj.z,2),obj.N_prof);
%     kappa = 0.39;
%     rng(37)
%     rand_h = rand(N_rand,1);
%     counter_1 = 1;
%     %Townsend-Perry constants
%     A = 1.26;
%     B = 1.00;
%     progressbar('Profile Generation')
%         for prof_num = 1:obj.N_prof
%             z_i = obj.z(1);
%             z_ind = 1;
%             counter_2 = 1;
%             while z_i < obj.z(end)
%                 %{
%                  Fitted model is avaliable in section 5.2:
%                  https://www.cambridge.org/
%                  core/journals/journal-of-fluid-mechanics/article/
%                  stochastic-modelling-of-the-instantaneous-velocity
%                  -profile-in-roughwall-turbulent-boundary-layers
%                  /492F3CD03C8C3E7ED306E9117B848B5E
%                 %}
%                 hm_i = z_i*exp(-3.59*(z_i/obj.delta)^0.91+1*sqrt(2)*erfinv(2*rand_h(counter_1)-1));
% 
%                 % In case if you want to remove one-step profiles
%                 % if z_i + hm_i > obj.z(end) & counter_2 == 1
%                 %     counter_1 = counter_1+1;
%                 %     continue
%                 % end
% 
%                 row = find(obj.z > z_i + hm_i,1,'first');
%                 if row
%                     % std_uu_plus = 2;
%                     std_uu_plus = sqrt(-A*log((z_i+hm_i/2)/obj.delta)+B);
%                     mean_u_plus = 1/kappa*log((z_i+hm_i/2)/obj.z_0);
%                 else
%                     std_uu_plus = sqrt(-A*log(((z_i+obj.z(end))/2)/obj.delta)+B);
%                     mean_u_plus = 1/kappa*log(((z_i+obj.z(end))/2)/obj.z_0);
%                     % std_uu_plus = sqrt(-A*log(obj.z(end)/obj.delta)+B);
%                     % mean_u_plus = 1/kappa*log(obj.z(end)/obj.z_0);
%                 end
% 
%                 um_i = obj.u_tau*(mean_u_plus+std_uu_plus*sqrt(2)*erfinv(2*rand_u(counter_1)-1));
%                 wm_i = obj.u_tau*(0+0.85*sqrt(2)*erfinv(2*rand_w(counter_1)-1));
% 
%                 if row
%                     log_data{prof_num}{counter_2} = ...
%                         dictionary(["hm_i[m]","zm_i[m]","um_i[m/s]","wm_i[m/s]"],...
%                         [hm_i, z_i+hm_i/2, um_i, wm_i]);
%                     Gen_u_prof(z_ind:row-1, prof_num) = um_i;
%                     Gen_w_prof(z_ind:row-1, prof_num) = wm_i;
%                     z_i = z_i + hm_i;
%                 else
%                     log_data{prof_num}{counter_2} = ...
%                         dictionary(["hm_i[m]","zm_i[m]","um_i[m/s]","wm_i[m/s]"],...
%                         [obj.z(end)-z_i, 0.5*(z_i+obj.z(end)), um_i, wm_i]);
%                     Gen_u_prof(z_ind:end, prof_num) = um_i;
%                     Gen_w_prof(z_ind:end, prof_num) = wm_i;
%                     z_i = obj.z(end);
%                 end
%                 z_ind = row;
%                 counter_2 = counter_2 + 1;
%                 counter_1 = counter_1 + 1;
%                 if counter_1 > numel(rand_u)
%                     error(['Not enough pre-generated random numbers for profile #%d. ' ...
%                         'Please increase the number of generated random samples. ' ...
%                         'Current limit: %d, required index: %d.'], ...
%                         prof_num, numel(rand_u), counter_1);
%                 end
%             end
%             progressbar((prof_num)/obj.N_prof)
% 
%         end
% end


function [Gen_u_prof, Gen_w_prof, log_data] = SG_VelProf(obj)
%SGVF generate step like velocity profiles. Detail of line will explain in
%the folllowing lines.
%#ok<*AGROW>
    N_rand = 1e6;
    repro_numb = 137; %for reproducibility, we generate one set of random numbers
    [rand_u,rand_w] = genGaussCop(obj.ro_uw, N_rand, repro_numb);
    Gen_u_prof = zeros(size(obj.z,2),obj.N_prof);
    Gen_w_prof = zeros(size(obj.z,2),obj.N_prof);
    kappa = 0.39;
    rng(37)
    rand_h = rand(N_rand,1);
    counter_1 = 1;
    %Townsend-Perry constants
    A = 1.26;
    B = 1.00;
    progressbar('Profile Generation')
        for prof_num = 1:obj.N_prof
            z_i = obj.z(1);
            z_ind = 1;
            counter_2 = 1;
            while z_i < obj.z(end)
                %{
                 Fitted model is avaliable in section 5.2:
                 https://www.cambridge.org/
                 core/journals/journal-of-fluid-mechanics/article/
                 stochastic-modelling-of-the-instantaneous-velocity
                 -profile-in-roughwall-turbulent-boundary-layers
                 /492F3CD03C8C3E7ED306E9117B848B5E
                %}
                hm_i = z_i*exp(-3.59*(z_i/obj.delta)^0.91+1*sqrt(2)*erfinv(2*rand_h(counter_1)-1));

                % In case if you want to remove one-step profiles
                % if z_i + hm_i > obj.z(end) & counter_2 == 1
                %     counter_1 = counter_1+1;
                %     continue
                % end
                
                row = find(obj.z > z_i + hm_i,1,'first');
                std_uu_plus = sqrt(-A*log((z_i+hm_i/2)/obj.delta)+B);
                mean_u_plus = 1/kappa*log((z_i+hm_i/2)/obj.z_0);

                if std_uu_plus^2 <= 0
                    z_i = obj.z(1);
                    z_ind = 1;
                    counter_2 = 1;
                    log_data{prof_num} = {};
                    counter_1 = counter_1 + 1;
                    if counter_1 > numel(rand_u)
                        error(['Not enough pre-generated random numbers for profile #%d ' ...
                            'after restarting. Please increase N_rand.'], prof_num);
                    end
                    continue
                end
                
                um_i = obj.u_tau*(mean_u_plus+std_uu_plus*sqrt(2)*erfinv(2*rand_u(counter_1)-1));
                wm_i = obj.u_tau*(0+0.85*sqrt(2)*erfinv(2*rand_w(counter_1)-1));
                
                if row
                    log_data{prof_num}{counter_2} = ...
                        dictionary(["hm_i[m]","zm_i[m]","um_i[m/s]","wm_i[m/s]"],...
                        [hm_i, z_i+hm_i/2, um_i, wm_i]);
                    Gen_u_prof(z_ind:row-1, prof_num) = um_i;
                    Gen_w_prof(z_ind:row-1, prof_num) = wm_i;
                    z_i = z_i + hm_i;
                else
                    log_data{prof_num}{counter_2} = ...
                        dictionary(["hm_i[m]","zm_i[m]","um_i[m/s]","wm_i[m/s]"],...
                        [obj.z(end)-z_i, 0.5*(z_i+obj.z(end)), um_i, wm_i]);
                    Gen_u_prof(z_ind:end, prof_num) = um_i;
                    Gen_w_prof(z_ind:end, prof_num) = wm_i;
                    z_i = obj.z(end);
                end
                z_ind = row;
                counter_2 = counter_2 + 1;
                counter_1 = counter_1 + 1;
                if counter_1 > numel(rand_u)
                    error(['Not enough pre-generated random numbers for profile #%d. ' ...
                        'Please increase the number of generated random samples. ' ...
                        'Current limit: %d, required index: %d.'], ...
                        prof_num, numel(rand_u), counter_1);
                end
            end
            progressbar((prof_num)/obj.N_prof)

        end
end