%{
To Stochatically Generate Velocity Profiles (SGVP), reorganize, and add vortex
cores, you need to provide initial physical parameters such as
u_{tau}[m/s](friction velocity), z_{0}[m](aerodynamic roughness length),
delta[m](boundary layer thickness), lambda[m](Taylor micro-scale). The details 
about the rest of the parameter will be providing in the following lines.
If you need details about the statistics, refer to the following papers:

**For SG of Velocity Profiles:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
stochastic-modelling-of-the-instantaneous-velocity-profile-in-roughwall-
turbulent-boundary-layers/492F3CD03C8C3E7ED306E9117B848B5E

**For SG of Velocity field:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
stochastic-modal-velocity-field-in-roughwall-turbulence/
800D97E60609F4FE4D92847CCDB1C7A0

**For SG of Vortex core:
coming soon
%}


root = fileparts(mfilename('fullpath'));  % project_root
addpath(fullfile(root, 'utils')); %addpath for using functions in utils.

%{
The physical properties of the Saddoughi High-velocity(50[m/s]) are used. 
You can find details about this dataset in the papers.
%}

u_tau = 2.325; %[m/s]
z_0 = 3.8e-4; %[m];
kappa = 0.39;
AFR = 8.5;
ks = z_0/exp(-kappa*AFR);
delta = 1.09; %[m]
lambda_T = 0.007; %[m]
nu = 1.57e-5; %[m^2/s]


%{
Define the wall-normal extent of the generated profiles which covers
logarithmic layer. The generated profiles start from z = 1.5*ks,
continues till z/delta = 0.25. The values here are conservatively considered.
The wall-normal resolution of profiles is selected 0.4*lambda/6. Since
average thickness of shear layer is 0.4*lambda, I divided it into 6 to
resolve the shear layer.
%}
z_min = 1.5*ks; %[m]
z_max = 0.25*delta; %[m]
res_z = 0.4*lambda_T/6; %[m]

% Let us generate an object. 
Gen_sample = stochastic_generation(u_tau, z_0, ks,delta, lambda_T, nu,z_min, z_max,...
    res_z);

%% Stochastic Generation of Velocity Profiles (SGVP)
%{
For the stochastic generation of velocity profiles, we need to start from
z_min, then stochastocally generate thickness (step), streamwise (or modal)
velocity, and vertical velocity. Repeat the process of step generation 
untill we reach to z_max, and then stop. We can independently generate 
profiles as much as we want. The minimum number of the profiles for having 
convergence in statistics (stationary dataset) is ~100 profiles. However, 
since we want to generate long velocity field in the next step, we need 
more profiles.
The parameter we need here is the Pearson correlation coefficient between
streamwise velocity and vertical velocity components. This parameter is
necessary in Reynolds shear stress (covariance) reproduction. This value in
the wind tunnel (m1) dataset is -0.4. You can tweak it w.r.t. your dataset.
N_prof is number of profiles we want to generate.
The parameters for the Townsend-Perry constant (A, B) can be set inside the
function.
%}

ro_uw = -0.33;
% Can also use Priyadarshana and Klevicki (2004) equation for the correlation coefficient:
U_inf = 50;
n = 4;
Theta = delta*n/((n+1)*(n+2));
Re_tau = u_tau*delta/nu;
% ro_uw = -0.03*log(U_inf*Theta/(u_tau*delta)*Re_tau)+0.63;

N_prof = 15e4;

Gen_sample.SGVP(ro_uw, N_prof)

%% Stochastic Generation of Velocity Field (SGVF)

%{ 
For the generated velocity field, you need to prescribe the distance
between the velocity profiles which here is lambda_T. The storage size M
is formulated through a Re dependent equation.
In the function, first, the shear layers between the adjacent UMZs are smoothed using
harmonic functions.
%}
Delta_x = lambda_T;
M = ceil(150*sqrt(Gen_sample.u_tau*Gen_sample.delta/Gen_sample.nu));
Gen_sample.SGVF(Delta_x, M)

%% Increasing the resolution of the SGVF and compute the lambda_ci field
%{
At the first step, we need to divide the VF into 1 * delta (L = 1)
intervals, and increase the resolution. This will be helpful in parallel
execution of this step.
For increasing the resolution, the factor that should be increased is 
another input parameter. Here Kr is the factor parameter which is 
consistent with the parameter in the paper.  All these steps are done 
increasing_res function.
The band-stop filter is implemented inside the function.
%}
L = 1; 
Kr = 10;
increasing_res(Gen_sample, L, Kr);
%% Lambda_ci(swirling strength) field
% In the next step, we need to find the Lambda_ci field (swirling strength)
% to insert the vortices in the shear layer. This code is implemented in
% parallel.
Lambda_ci(Gen_sample)

%% Generating radius, maximum azimuthal velocity, and corrcoeff(u,w)
%{
At this step, we genrate radius(r_{\omega}), maximum azimuthal velocity
(u_{\omega}), and the correlation coefficient between u and w component
of the generated vortex \rho^{VorX}_{u,w}. Generated values for r_omega
are normalized by lambda_T, for u_omega are normalized by u_tau.
For the decile matching algorithm, the vorticity of the generated vortices
is also computed.
%}

[u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
    u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,...
    omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,...
    omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles] = VorX_generator(Gen_sample);


%% Find the c.d.f. of the omega distribution
%{
At this step, we find the c.d.f. of the omega distribution in the
non_zero lambda_ci regions. The output of this step is used in decile
matching process for adaptive selection of vortices. This step is also done
in parallel.
%}
cdf_omega(Gen_sample)

%% Adding Stochastically Generated Vortex Cores to the HRVF (both near wall and shear layers). 

Gen_sample.SGVorX(u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
    u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,...
    omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,...
    omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles)

%% computing std of the small-scale 2-D Gaussian filter

w0 = 1e-5;
li = 0.4*Gen_sample.lambda_T;%0.4
f = @(s) 2*(-2*s.*log(2*w0*pi*s))^0.5 - li;

% Solve using fzero with an initial guess
s0 = 1e-7;   % initial guess(Change it if it is necessary)
sol = fzero(f, s0);

% Display solution
disp(sol)

%The next line is just a check to verify the correctness of computed std.
% li should be equal to 0.4*lambda
% li = 2*(-2*sol*log(2*pi*w0*sol))^0.5;

Sigma_small = [sol 0;...
    0 sol];

%% Convolving the small-scale 2-D Gaussian filter
kern_xsize = 0.4;
kern_ysize = 0.4;
Gen_sample.Gauss_kernel(kern_xsize, kern_ysize, Sigma_small)

%% Spectral and structure function analysis
%{
The spatial and spectral analysis can be done on different stages of the
generated velocity field.
LRVF, HRVF(w/ or w/o band-stop filter), HRVF(w/ or w/o inserting VorX),
HRVF(w/ or w/o passing 2-D filter).
If you want to compute it on LRVF, run increasing_res(Gen_sample, L, Kr)
with Kr = 1 firts, then computed the spectra or structure function.
%}
Gen_sample.Spectral_analysis_E11
Gen_sample.Spectral_analysis_E22
Gen_sample.Spectral_analysis_E12
Gen_sample.Structure_function_D11
Gen_sample.Structure_function_D22

%% Computing production

production_comp(Gen_sample)
%% Saving

save('C:\Users\ehsan010\Desktop\Gen_sample_largerstorage80_TPmodel.mat','Gen_sample')
% % save('Gen_sample.mat','Gen_sample')

%% Selective saving

% x = Gen_sample.Gen_x_HRVF;
% z =Gen_sample.z;
% u = Gen_sample.Gen_u_HRVF;
% w = Gen_sample.Gen_w_HRVF;
% delta = Gen_sample.delta;
% u_tau = Gen_sample.u_tau;
% z_0 = Gen_sample.z_0;
% lambda_T = Gen_sample.lambda_T;
% nu = Gen_sample.nu;
% 
% save('C:\Users\ehsan010\Desktop\SGHRVFVorX_Saddoughi(Corrected_U).mat', ...
%     'x', ...
%     'z', ...
%     'u', ...
%     'w', ...
%     'delta', ...
%     'u_tau', ...
%     'z_0', ...
%     'lambda_T', ...
%     'nu')
