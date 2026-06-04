classdef stochastic_generation < dynamicprops
    %In this class, we define methods for SGVP, SGVF, SG of Vortex cores
    %and addition. The properties of the object is included as well.
    
    properties
        u_tau, z_0, ks, delta, lambda_T, nu, ro_uw, N_prof, z, Delta_x,...
            Gen_u_prof, Gen_w_prof, log_data_SGVP,...
            Gen_u_LRVF, Gen_w_LRVF, Gen_x_LRVF, hist_corr,...
            Kr, Gen_u_HRVF, Gen_w_HRVF, Gen_x_HRVF,...
            Lambda_ci, Lambda_cirms,omega, omega_rms,...
            ro_r_u_omega_near_wall, ro_r_u_omega_far_wall,...
            prograde_VorX_data, retrograde_VorX_data,...
            E11, k1, E22, E12, epsilon_E11, eta_E11,...
            D11, r11, epsilon_D11, eta_D11, D22, cdf_omega, mean_omega,...
            production, eta_prod
    end
    
    methods
        function obj = stochastic_generation(u_tau, z_0, ks,delta, lambda_T, nu, ...
                z_min, z_max, res_z)
            %STOCHASTIC_GENERATION Construct an instance of this class
            obj.u_tau = u_tau;
            obj.z_0 = z_0;
            obj.ks = ks;
            obj.delta = delta;
            obj.lambda_T = lambda_T;
            obj.nu = nu;
            obj.z = z_min:res_z:z_max;
        end
        
        function SGVP(obj, ro_uw, N)
            %In this function, we generate N profiles stochastically.
            %   Detailed explanation goes here
            obj.ro_uw = ro_uw;
            obj.N_prof = N;
            [obj.Gen_u_prof, obj.Gen_w_prof, obj.log_data_SGVP] =...
                SG_VelProf(obj);

            
        end
        function SGVF(obj, Delta_x, M)
            %In this function, we try to reorganize the profiles generated
            %in SGVP, to have correlated velocity field.
            %   Detailed explanation goes here
            obj.Delta_x = Delta_x;
            [obj.Gen_u_LRVF, obj.Gen_w_LRVF, obj.hist_corr] =...
                SG_VelField(obj, M);

            
        end
        function SGVorX(obj, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
            u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,...
            omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,...
            omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles)
            % In this function, we try to add vortices to the SGVF.
            %   Detailed explanation goes here
            SG_VorX(obj, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
                u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,...
                omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,...
                omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles);
            
        end
        function Structure_function_D11(obj)
            str_func_D11(obj)
            
        end
        function Structure_function_D22(obj)
            str_func_D22(obj)
            
        end
        function Spectral_analysis_E11(obj)
            spec_analysis_E11(obj) 

        end
        function Spectral_analysis_E22(obj)
            spec_analysis_E22(obj) 

        end
        function Spectral_analysis_E12(obj)
            spec_analysis_E12(obj) 

        end
        function Gauss_kernel(obj,kern_xsize, kern_ysize, Sigma)
            % Applying Multivariate Gaussian Kernel
            MVG_kernel(obj, kern_xsize, kern_ysize, Sigma)
            
        end
    end
end

