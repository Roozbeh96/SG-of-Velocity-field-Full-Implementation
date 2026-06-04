%% Load Generated data

load("C:\Users\ehsan010\Desktop\Gen_sample_largerstorage10_TPmodel(u_Fourierfilter).mat")


%% Read Saddoughi data(Both 10 and 50 [m/s])
jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\mid-layer-u-spectra(Ue_10).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eta_10_Sadd = str2double({jsonData.x});
energy_10_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\inner-layer-compen-u-spectra(Ue_10).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eta_compen_inner_10_Sadd = str2double({jsonData.x});
compen_energy_inner_10_Sadd = str2double({jsonData.y});


jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\mid-layer-compen-u-spectra(Ue_10).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eta_compen_mid_10_Sadd = str2double({jsonData.x});
compen_energy_mid_10_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\profile_meanvelocity.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_U_o_Ue_10_Sadd = str2double({jsonData.x});
U_o_Ue_10_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\profile_uvar.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_uu_o_Ue2_10_Sadd = str2double({jsonData.x});
uu_o_Ue2_10_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\profile_wvar.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_ww_o_Ue2_10_Sadd = str2double({jsonData.x});
ww_o_Ue2_10_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\profile_uwcovar.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_uw_o_Ue2_10_Sadd = str2double({jsonData.x});
uw_o_Ue2_10_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=10[ms]\mid-layer-D11(Ue_10).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
r_o_eta_10_Sadd = str2double({jsonData.x});
compen_D11_10_Sadd = str2double({jsonData.y});

Ue_10_Sadd = 10;
u_tau_10_Sadd = 0.465;
z_0_10_Sadd = 3.8e-4;
z_mid_o_delta_10_Sadd = 0.515/1.09;
z_inner_o_delta_10_Sadd = 0.1/1.09;
delta_10_Sadd = 1.09;
%------------------------------
jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\mid-layer-u-spectra(Ue_50).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eta_50_Sadd = str2double({jsonData.x});
energy_50_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\inner-layer-uw-cospectra(Ue_50).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eps_S_50_Sadd = str2double({jsonData.x});
E12_50_Sadd = str2double({jsonData.y});


jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\inner-layer-compend-u-spectra(Ue_50).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eta_compen_inner_50_Sadd = str2double({jsonData.x});
compen_energy_inner_50_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\inner-layer-compend-w-spectra(Ue_50).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eta_E22_compen_inner_50_Sadd = str2double({jsonData.x});
compen_E22_inner_50_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\mid-layer-compend-u-spectra(Ue_50).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
k1_eta_compen_mid_50_Sadd = str2double({jsonData.x});
compen_energy_mid_50_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\profile_meanvelocity.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_U_o_Ue_50_Sadd = str2double({jsonData.x});
U_o_Ue_50_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\profile_uvar.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_uu_o_Ue2_50_Sadd = str2double({jsonData.x});
uu_o_Ue2_50_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\profile_wvar.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_ww_o_Ue2_50_Sadd = str2double({jsonData.x});
ww_o_Ue2_50_Sadd = str2double({jsonData.y});

jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\profile_uwcovar.json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
z_uw_o_Ue2_50_Sadd = str2double({jsonData.x});
uw_o_Ue2_50_Sadd = str2double({jsonData.y});


jsonFileName = "G:\My Drive\Research\SG-of-Velocity-field-Full-Implementation-\Matlab\src\Data\Saddoughi\Ue=50[ms]\mid-layer-D11(Ue_50).json"; 
jsonStr = fileread(jsonFileName);
jsonData = jsondecode(jsonStr);
r_o_eta_50_Sadd = str2double({jsonData.x});
compen_D11_50_Sadd = str2double({jsonData.y});

Ue_50_Sadd = 50;
u_tau_50_Sadd = 2.325;
z_0_50_Sadd = 3.8e-4;
z_mid_o_delta_50_Sadd = 0.4/1.09;
z_inner_o_delta_50_Sadd = 0.1/1.09;
delta_50_Sadd = 1.09;
