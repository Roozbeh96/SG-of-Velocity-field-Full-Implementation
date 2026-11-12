% Figures 5, 6 and 7 of the paper: the generated velocity field (Gen_sample)
% compared with the Saddoughi data (U_inf = 50 [m/s]). Load Gen_sample and
% run read_data.m first; then run each section to draw its figure.

%% Mean velocity and second moments of the HRVF vs Saddoughi (Figure 5 of the paper)
% (a) mean velocity, (b) streamwise variance versus z/delta, (c) vertical
% variance, (d) Reynolds shear stress. Saddoughi data for U_inf = 50 [m/s].
kappa = 0.39;
saddFace = [1.00,0.81,0.81]; % light face colour of the Saddoughi markers
gold = [0.93,0.69,0.13];
numPoints = 30; % Adjust as needed
logIndices = unique(round(logspace(0, log10(size(Gen_sample.z,2)), numPoints)));
logIndices_b = unique(round(logspace(0, log10(size(Gen_sample.z,2)-50), numPoints)));
s2 = Ue_50_Sadd^2/u_tau_50_Sadd^2;

U_Gen = mean(mean(Gen_sample.Gen_u_HRVF,3),2);
uu_Gen = mean(var(Gen_sample.Gen_u_HRVF,0,3),2);
ww_Gen = mean(var(Gen_sample.Gen_w_HRVF,0,3),2);
uw_Gen = mean(-sum((Gen_sample.Gen_u_HRVF - mean(Gen_sample.Gen_u_HRVF,3)) ...
    .*(Gen_sample.Gen_w_HRVF - mean(Gen_sample.Gen_w_HRVF,3)),3) ...
    / (size(Gen_sample.Gen_u_HRVF,3) - 1),2);

saddStyle = {'LineStyle','none','LineWidth',1.5,'Marker','o','MarkerSize',8,...
    'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',saddFace};
genStyle = {'Color',[1,0,0],'LineWidth',2,'Marker','x','MarkerSize',8,...
    'MarkerEdgeColor',[0,0,0]};
axStyle = {'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on',"XMinorGrid","on","YMinorGrid","on"};

figure
set(gcf,'Position',[885,209,798,902])

% (a) Mean velocity
axes('Position',[0.0794,0.7761,0.8997,0.2139])
semilogy(1/kappa*log(30:3000),30:3000,'Color',gold,'LineWidth',2.5)
hold on
semilogy(U_o_Ue_50_Sadd*Ue_50_Sadd/u_tau_50_Sadd, z_U_o_Ue_50_Sadd/z_0_50_Sadd, saddStyle{:})
semilogy(U_Gen(logIndices)/Gen_sample.u_tau, Gen_sample.z(logIndices)/Gen_sample.z_0, genStyle{:})
xlim([8 22])
ylim([5 2e4])
set(gca,axStyle{:},'YTick',[1e1 1e2 1e3 1e4])
xlabel('$\overline{U}/u_{\tau}$','Interpreter','latex','FontSize',15)
ylabel('$z/z_{0}$','Interpreter','latex','FontSize',15)
legend('','Saddoughi($U_{\infty}=50$[m/s])','Gen($U_{\infty}=50$[m/s])',...
    'Interpreter','latex','FontSize',12,'Position',[0.1621,0.9268,0.2818,0.0524],...
    'NumColumns',1,'Orientation','vertical','Color',[1,1,1]);

% (b) Streamwise variance versus z/delta
axes('Position',[0.0802,0.4992,0.8990,0.2095])
semilogx(logspace(-2, 0, 20),-1.26*log(logspace(-2, 0, 20))+1.0,'Color',gold,'LineWidth',2.5)
hold on
semilogx(z_uu_o_Ue2_50_Sadd/delta_50_Sadd, uu_o_Ue2_50_Sadd*s2, saddStyle{:})
semilogx(Gen_sample.z(logIndices_b)/Gen_sample.delta, uu_Gen(logIndices_b)/Gen_sample.u_tau^2, genStyle{:})
xlim([1e-3 1.335])
ylim([0 8])
set(gca,axStyle{:},'XTick',[1e-3 1e-2 1e-1 1e0])
xlabel('$z/\delta$','Interpreter','latex','FontSize',15)
ylabel('$\overline{u^{\prime^{2}}}/u_{\tau}^{2}$','Interpreter','latex','FontSize',15)

% (c) Vertical variance
axes('Position',[0.0789,0.2794,0.9002,0.1682])
semilogx(z_ww_o_Ue2_50_Sadd/z_0_50_Sadd, ww_o_Ue2_50_Sadd*s2, saddStyle{:})
hold on
semilogx(Gen_sample.z(logIndices)/Gen_sample.z_0, ww_Gen(logIndices)/Gen_sample.u_tau^2, genStyle{:})
xlim([10 5000])
ylim([0 2])
set(gca,axStyle{:},'XTick',[1e1 1e2 1e3])
xlabel('$z/z_{0}$','Interpreter','latex','FontSize',15)
ylabel('$\overline{w^{\prime^{2}}}/u_{\tau}^{2}$','Interpreter','latex','FontSize',15)

% (d) Reynolds shear stress
axes('Position',[0.0815,0.0401,0.8977,0.1927])
semilogx(z_uw_o_Ue2_50_Sadd/z_0_50_Sadd, -uw_o_Ue2_50_Sadd*s2, saddStyle{:})
hold on
semilogx(Gen_sample.z(logIndices)/Gen_sample.z_0, uw_Gen(logIndices)/Gen_sample.u_tau^2, genStyle{:})
xlim([10 5000])
ylim([-0.5 1.6])
set(gca,axStyle{:},'XTick',[1e1 1e2 1e3],'YTick',[-0.5 0 0.5 1 1.5])
xlabel('$z/z_{0}$','Interpreter','latex','FontSize',15)
ylabel('$-\overline{u^{\prime}w^{\prime}}/u_{\tau}^{2}$','Interpreter','latex','FontSize',15)

% Panel labels
labelPos = [0.0847,0.9481; 0.0759,0.6665; 0.0772,0.4049; 0.0847,0.1898];
labels = {'(a)','(b)','(c)','(d)'};
for p = 1:4
    annotation(gcf,'textbox',[labelPos(p,:),0.0522,0.0458],'String',labels{p},...
        'Interpreter','latex','FontSize',15,'FitBoxToText','off','EdgeColor','none');
end

%% Spectra and structure function of the HRVF + VorX (Figure 6 of the paper)
% (a) compensated E11 at z/delta = 0.09, (b) E11 in Kolmogorov scaling at
% z/delta = 0.14, (c) D11 at z/delta = 0.05, 0.09 and 0.14 with the
% premultiplied E11 at the same elevations as an inset. The lines are
% sampled at numPoints log-spaced wavenumbers. The noisy spectra are shown as
% dots (at numPoints_dots wavenumbers in (a) and numPoints_dots_inset in the
% inset) with their moving average as a line: movmean(.,[2,2]) in (a) and
% movmean(.,[1,1]) in the inset; (b) is not smoothed.
[~, ind_05] = find(Gen_sample.z/Gen_sample.delta>= 0.05, 1, 'first');
[~, ind_09] = find(Gen_sample.z/Gen_sample.delta>= 0.0919, 1, 'first');
[~, ind_14] = find(Gen_sample.z/Gen_sample.delta>= 0.14, 1, 'first');
kappa = 0.39;
epsilon_Sadd = Gen_sample.u_tau^3./(kappa*Gen_sample.z');
eta_Sadd = (Gen_sample.nu^3./epsilon_Sadd).^0.25;
saddFace = [1.00,0.81,0.81]; % light face colour of the Saddoughi markers
k1 = Gen_sample.k1(:);
numPoints = 50; % Adjust as needed
logIndices = unique(round(logspace(0, log10(size(k1,1)), numPoints)));
numPoints_dots = 1000; % Adjust as needed (dots of figure (a))
logIndices_dots = unique(round(logspace(0, log10(size(k1,1)), numPoints_dots)));
numPoints_dots_inset = 500; % Adjust as needed (dots of the inset of figure (c))
logIndices_dots_inset = unique(round(logspace(0, log10(size(k1,1)), numPoints_dots_inset)));
genLabel = @(ind) sprintf('Gen($U_{\\infty}=50$[m/s]) HRVF + VorX $z/\\delta$ = %.2f',...
    Gen_sample.z(ind)/Gen_sample.delta);
saddLabel = @(zod) sprintf('Saddoughi($\\mathrm{U}_{\\mathrm{\\infty}}=50$[m/s]) z$/\\delta$ = %.2f', zod);
axStyle = {'TickLabelInterpreter','latex','FontSize',15,'Box','on',...
    'XGrid','on','YGrid','on',"XMinorGrid","on","YMinorGrid","on"};
legStyle = {'Interpreter','latex','FontSize',12,'NumColumns',1,'Orientation','vertical',...
    'Color',[1,1,1]};

figure
set(gcf,'Position',[848,362,984,620])

% (a) Compensated E11 at z/delta = 0.09
axes('Position',[0.0681,0.6226,0.5504,0.3538])
semilogx(linspace(1e-5,1e5,2),1.5*18/55*ones(1,2),'k-.','LineWidth',2)
hold on
semilogx(k1_eta_compen_inner_50_Sadd,compen_energy_inner_50_Sadd,'LineStyle','none',...
    'LineWidth',1.5,'Marker','v','MarkerSize',9,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',saddFace)
semilogx(k1(logIndices_dots)*0.055e-3,...
    Gen_sample.E11(ind_09,logIndices_dots).*k1(logIndices_dots)'.^(5/3)*342^(-2/3),...
    'LineStyle','none','Marker','.','MarkerSize',4,'Color',[1,0,0])
y = Gen_sample.E11(ind_09,logIndices).*k1(logIndices)'.^(5/3)*342^(-2/3);
semilogx(k1(logIndices)*0.055e-3, movmean(y,[2,2]),'-','LineWidth',2,'Color',[1,0,0],...
    'Marker','x','MarkerSize',8,'MarkerEdgeColor',[0,0,0])
xlim([1e-5 1])
ylim([0 1.5])
set(gca,axStyle{:},'YTick',[0 0.5 1 1.5])
xlabel('$k_{1}\eta$','Interpreter','latex','FontSize',15)
ylabel('$\varepsilon^{-2/3}k_{1}^{5/3}E_{11}(k_{1})$','Interpreter','latex','FontSize',15)
legend('',saddLabel(z_inner_o_delta_50_Sadd),'',genLabel(ind_09),legStyle{:},...
    'Position',[0.0785,0.8694,0.3059,0.0984]);

% (b) E11 in Kolmogorov scaling at z/delta = 0.14
axes('Position',[0.6694,0.0978,0.3203,0.8892])
loglog((3e-4:1e-5:1),(3e-4:1e-5:1).^(-5/3),'LineWidth',2,'Color',[1.00,0.74,0.17])
hold on
loglog(linspace(5e-7,5,2),[1 1],'k--','LineWidth',1.5)
loglog([1 1],[1e-8 10],'k--','LineWidth',1.5)
loglog(k1_eta_50_Sadd,energy_50_Sadd,'LineStyle','none','LineWidth',1.5,'Marker','o',...
    'MarkerSize',9,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',saddFace)
loglog(k1(logIndices)*eta_Sadd(ind_14,1),...
    Gen_sample.E11(ind_14,logIndices)/(epsilon_Sadd(ind_14,1)*Gen_sample.nu^5)^(1/4),'-',...
    'LineWidth',2,'Color',[0,0,1],'Marker','x','MarkerSize',8,'MarkerEdgeColor',[0,0,0])
ylim([1e-4 3.2e6])
set(gca,axStyle{:},'YTick',[1e0 1e6],'XTick',[1e-4 1e0])
xlabel('$k_{1}\eta$','Interpreter','latex','FontSize',15)
ylabel('$E_{11}(\mathrm{k}_{1})/(\varepsilon \nu^5)^{1/4}$','Interpreter','latex','FontSize',15)
annotation(gcf,'textbox',[0.7993,0.8618,0.0847,0.0802],'String','$(k_{1}\eta)^{-5/3}$',...
    'Interpreter','latex','FontSize',15,'FitBoxToText','off','EdgeColor','none');
legend('','','',saddLabel(z_mid_o_delta_50_Sadd),genLabel(ind_14),legStyle{:},...
    'Position',[0.6774,0.1150,0.3059,0.0984]);

% (c) D11 at z/delta = 0.05, 0.09 and 0.14
axes('Position',[0.0681,0.0926,0.5504,0.4283])
semilogx(Gen_sample.r11/Gen_sample.z(ind_05), mean(Gen_sample.D11(ind_05,:,:),3)/Gen_sample.u_tau^2,...
    ':','LineWidth',2.5,'Color',[0,0,0])
hold on
semilogx(Gen_sample.r11/Gen_sample.z(ind_09), mean(Gen_sample.D11(ind_09,:,:),3)/Gen_sample.u_tau^2,...
    '-','LineWidth',2.5,'Color',[1,0,0])
semilogx(Gen_sample.r11/Gen_sample.z(ind_14), mean(Gen_sample.D11(ind_14,:,:),3)/Gen_sample.u_tau^2,...
    '--','LineWidth',2.5,'Color',[0,0,1])
xx = logspace(0.05, 0.5, 25);
plot(xx, 2*1.26*log(xx)+4.2,'-','LineWidth',2.5,'Color',[0.93,0.69,0.13])
xlim([0.1 20])
ylim([0 14])
set(gca,axStyle{:},'YTick',0:2:14)
xlabel('$r/z$','Interpreter','latex','FontSize',15)
ylabel('$D_{11}(r)/u_{\tau}^2$','Interpreter','latex','FontSize',15)
legend(genLabel(ind_05),genLabel(ind_09),genLabel(ind_14),'',legStyle{:},...
    'Position',[0.3074,0.1083,0.3059,0.0888]);
annotation(gcf,'textarrow',[0.5697 0.5700],[0.4565 0.3000],'String','Increasing $z$',...
    'Interpreter','latex','FontSize',15);

% Inset of (c): premultiplied E11 at the same elevations
axes('Position',[0.1011,0.3199,0.2905,0.1630])
inds = [ind_05, ind_09, ind_14];
cols = {[0,0,0],[1,0,0],[0,0,1]};
lstyles = {':','-','-.'};
for j = 1:3
    % dots at numPoints_dots_inset wavenumbers, smoothed line at numPoints
    semilogx(k1(logIndices_dots_inset)*Gen_sample.z(inds(j)),...
        Gen_sample.E11(inds(j),logIndices_dots_inset).*k1(logIndices_dots_inset)'/Gen_sample.u_tau^2,...
        'LineStyle','none','Marker','.','MarkerSize',4,'Color',cols{j})
    hold on
    y = Gen_sample.E11(inds(j),logIndices).*k1(logIndices)'/Gen_sample.u_tau^2;
    semilogx(k1(logIndices)*Gen_sample.z(inds(j)), movmean(y,[1,1]),lstyles{j},...
        'LineWidth',2,'Color',cols{j})
end
xlim([1e-2 1e2])
ylim([0 2])
set(gca,axStyle{:},'FontSize',13,'XAxisLocation','top','YAxisLocation','right','YTick',[0 1 2])
xlabel('$k_{1}z$','Interpreter','latex','FontSize',13,'Units','normalized','Position',[0.5,-0.12,0])
ylabel('$k_{1}E_{11}(k_{1})/u_{\tau}^{2}$','Interpreter','latex','FontSize',13,...
    'Units','normalized','Position',[-0.07,0.5,0])
annotation(gcf,'textarrow',[0.1363 0.1831],[0.4200 0.3452],'String','Increasing $z$',...
    'Interpreter','latex','FontSize',13);

% Panel labels
annotation(gcf,'textbox',[0.0032,0.9403,0.0414,0.0630],'String','(a)','Interpreter','latex',...
    'FontSize',15,'FitBoxToText','off','EdgeColor','none');
annotation(gcf,'textbox',[0.9341,0.9131,0.0414,0.0630],'String','(b)','Interpreter','latex',...
    'FontSize',15,'FitBoxToText','off','EdgeColor','none');
annotation(gcf,'textbox',[0.0080,0.5225,0.0414,0.0630],'String','(c)','Interpreter','latex',...
    'FontSize',15,'FitBoxToText','off','EdgeColor','none');

%% Co-spectrum E12 of the HRVF + VorX at z/delta = 0.09 (Figure 7 of the paper)
% Normalized with epsilon = 342 and S = 83.5 at z/delta = 0.09 (Saddoughi).
% The raw values are shown as dots at numPoints_dots wavenumbers and their
% moving average movmean(.,[2,2]) at numPoints wavenumbers as a line.
[~, ind_ii] = find(Gen_sample.z/Gen_sample.delta>= 0.0919, 1, 'first');
saddFace = [1.00,0.81,0.81]; % light face colour of the Saddoughi markers
numPoints = 50; % Adjust as needed
logIndices = unique(round(logspace(0.3, log10(size(Gen_sample.k1, 1)), numPoints)));
numPoints_dots = 1000; % Adjust as needed
logIndices_dots = unique(round(logspace(0.3, log10(size(Gen_sample.k1, 1)), numPoints_dots)));

figure
set(gcf,'Position',[1136,320,577,386])
axes('Position',[0.1263,0.1528,0.8344,0.8057])
loglog((1e-1:1e1:1e3),(1e-1:1e1:1e3).^(-7/3),'LineWidth',2,'Color',[1.00,0.74,0.17],'LineStyle','-')
hold on
loglog(k1_eps_S_50_Sadd,E12_50_Sadd,'LineStyle','none','LineWidth',1.5,'Marker','v',...
    'MarkerSize',9,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',saddFace)
loglog(Gen_sample.k1(logIndices_dots,1)*342^(0.5)*83.5^(-3/2),...
    abs(Gen_sample.E12(ind_ii,logIndices_dots))/(342^3*83.5^(-5))^(1/2),'LineStyle','none',...
    'Marker','.','MarkerSize',4,'Color',[1.0, 0.0, 0.0])
loglog(Gen_sample.k1(logIndices,1)*342^(0.5)*83.5^(-3/2),...
    movmean(abs(Gen_sample.E12(ind_ii,logIndices))/(342^3*83.5^(-5))^(1/2),[2,2]),'LineStyle','-',...
    'LineWidth',2,'Color',[1.0, 0.0, 0.0],'Marker','x','MarkerSize',8,'MarkerEdgeColor',[0,0,0])
xlim([1e-3 1e2])
ylim([1e-5 1e2])
set(gca,'TickLabelInterpreter','latex','FontSize',15,'Box','on','XGrid','on','YGrid','on',...
    "XMinorGrid","on","YMinorGrid","on",'YTick',[1e-4 1e-2 1e0 1e2])
xlabel('$k_{1}\epsilon^{1/2}S^{-3/2}$','Interpreter','latex','FontSize',15)
ylabel('$-E_{12}(k_{1})/(\varepsilon^{3} S^{-5})^{1/2}$','Interpreter','latex','FontSize',15)
annotation(gcf,'textbox',[0.5775,0.8616,0.2405,0.0824],...
    'String',{'$(k_{1}\epsilon^{1/2}S^{-3/2})^{-7/3}$'},'Interpreter','latex',...
    'FontSize',15,'FitBoxToText','off','EdgeColor','none');
legend('',sprintf('Saddoughi($U_{\\infty}=50$[m/s]) $z/\\delta$ = %.2f',z_inner_o_delta_50_Sadd),...
    '',sprintf('Gen($U_{\\infty}=50$[m/s]) $z/\\delta$ = %.2f',Gen_sample.z(ind_ii)/Gen_sample.delta),...
    'Interpreter','latex','FontSize',12,'Position',[0.1613,0.1885,0.5107,0.1225],...
    'NumColumns',1,'Orientation','vertical','Color',[1,1,1]);
