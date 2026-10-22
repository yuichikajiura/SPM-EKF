clc;
clear;
close all;

%% Parameters
% LiCo02.

run params_setting % load params

% For LiCo02, with voltage lim [2.5, 4.05], n_Li_s = 2.50, then capacity is
% 29.601016543579075 Ah/m^2
OneC = 29.601016543579075; %[Ah/m^2]

%%% Finite Difference Model for particle
p.Nr = 15;
Nr = p.Nr;
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;
r_vec = (0:p.delta_r_n:1)';
r_vecx = r_vec(2:end-1);

%%% Finite Difference points along x-coordinate
p.Nxn = 10;
p.Nxs = 5;
p.Nxp = 10;
p.Nx = p.Nxn + p.Nxs + p.Nxp;  
Nx = p.Nx - 3;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

% Model setup
[An,Ap,Bn,Bp,Cn,Cp,Dn,Dp,MN_mats] = spm_plant_obs_mats_YK_modified(p);
sys_n = ss(An,Bn,Cn,Dn);
sys_p = ss(Ap,Bp,Cp,Dp);

% Data setup
data = readtable('data/SPM_FDM_nr100_simulation_US06_Extended.csv');
data.Properties.VariableNames = {'t', 'I', 'V_spm', 'css_n', ...
    'cs_ave_n', 'css_p', 'cs_ave_p'};
I = data.I;
v_true = data.V_spm;
t = data.t;

noise = 0.03;  % the magnitude of white noise to standard deviation of data
std_current = 38;  % standard deviation of current data in Ampere
std_voltage = 0.16;  % standard deviation of voltage data in Volt
I = I + normrnd(0, std_current * noise, size(I));
v_true = v_true + normrnd(0, std_voltage * noise, size(v_true));

input_data = [t, I, v_true];

% Set nLi as identified (true value: 2.5)
p.n_Li_s = 2.4655;

% Initial Condition
[cs_n0,cs_p0] = init_cs(p, 4.0);
c_p0 = cs_p0 * ones(Nr-1, 1);
c_n0 = cs_n0 * ones(Nr-1, 1);
Sig_0 = eye(2*(Nr-1), 2*(Nr-1));
states0 = [c_p0; c_n0; reshape(Sig_0, [2*(Nr-1)*2*(Nr-1), 1])];

% calculate the size of the unit current noise with respect to the surface
% concentration
coef_I_n = - (2*p.delta_r_n * p.R_s_n)/(p.D_s_n * p.Faraday * p.a_s_n * p.Area * p.L_n); %positive noise is discharge, making Cs,n smaller than true
coef_I_p = (2*p.delta_r_p * p.R_s_p)/(p.D_s_p * p.Faraday * p.a_s_p * p.Area * p.L_p);

% Electrochemical Model Parameters
alpha_n = p.D_s_n / (p.R_s_n * p.delta_r_n)^2;
alpha_p = p.D_s_p / (p.R_s_p * p.delta_r_p)^2;

% Noise Covariances (calculated based on d(d(Cs,n)/dt)/dI)
Lp = zeros(Nr-1, 1);
Lp(Nr-1) = ((Nr-2)/(Nr-1)-2/4+Nr/(Nr-1)*(1/3)) * alpha_p * coef_I_p * std_current * noise;
Lp(Nr-2) = (-2-(Nr-1)/(Nr-2)*(1/4))* alpha_p * coef_I_p * std_current * noise;
Lp(Nr-3) = (Nr-2)/(Nr-3) * alpha_p * coef_I_p * std_current * noise;
Ln = zeros(Nr-1, 1);
Ln(Nr-1) = ((Nr-2)/(Nr-1)-2/4+Nr/(Nr-1)*(1/3)) * alpha_n * coef_I_n * std_current * noise;
Ln(Nr-2) = (-2-(Nr-1)/(Nr-2)*(1/4)) * alpha_n * coef_I_n * std_current * noise;
Ln(Nr-3) = (Nr-2)/(Nr-3) * alpha_n * coef_I_n * std_current * noise;
L = [Lp; Ln];
W = L*L';
N = [(std_voltage * noise)^2, 0; 0, 0.05^2];

% Simulate Extended Kalman Filter
[tsim,states] = ode45(@(t,x) ode_ekf_spmfdm(t,x,Ap, An, Bp, Bn, input_data, W, N, p),t,states0);

% Parse States
c_p_start = 1;
c_p_end = c_p_start + Nr-2;
c_n_start = c_p_end + 1;
c_n_end = c_n_start + Nr-2;
c_p_hat = states(:, c_p_start:c_p_end);
c_n_hat = states(:, c_n_start:c_n_end);
csp_hat = 4/3 * c_p_hat(:, Nr-1) - 1/3 * c_p_hat(:, Nr-2) ...
    +(2*p.delta_r_p * p.R_s_p)/(3*p.D_s_p * p.Faraday * p.a_s_p * p.Area * p.L_p)*I;
csn_hat = 4/3 * c_n_hat(:, Nr-1) - 1/3 * c_n_hat(:, Nr-2) ...
    -(2*p.delta_r_n * p.R_s_n)/(3*p.D_s_n * p.Faraday * p.a_s_n * p.Area * p.L_n)*I;
c_p0_hat = 4/3 * c_p_hat(:, 1) - 1/3 * c_p_hat(:, 2);
c_n0_hat = 4/3 * c_n_hat(:, 1) - 1/3 * c_n_hat(:, 2);
cs_ave_n_hat = calc_c_ave([c_n0_hat, c_n_hat, csn_hat]);
cs_ave_p_hat = calc_c_ave([c_p0_hat, c_p_hat, csp_hat]);


NT = size(t);
V_hat = zeros(NT);
nLi_hat = zeros(NT);
for k = 1:NT
    V_hat(k) = nonlinearSPMOutputVoltage_Scott(p,csn_hat(k),csp_hat(k),p.c_e,p.c_e,p.c_e,I(k));
    nLi_n_hat = cs_ave_n_hat(k) * p.Area * p.L_n * p.epsilon_s_n;
    nLi_p_hat = cs_ave_p_hat(k) * p.Area * p.L_p * p.epsilon_s_p;
    nLi_hat(k) = nLi_n_hat + nLi_p_hat;
end

% Plot results
fs = 15;    % Font Size for plots
figure(1); clf;
set(gcf,'Position',[99,29,300,300]);

subplot(321);
plot(t,csp_hat,'b-',t,data.css_p,'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{s,p}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])

subplot(322);
plot(t,csn_hat,'b-',t,data.css_n,'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{s,n}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])

subplot(323);
plot(t,cs_ave_p_hat,'b-',t,data.cs_ave_p,'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{ave,p}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])

subplot(324);
plot(t,cs_ave_n_hat,'b-',t,data.cs_ave_n,'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{ave,n}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])

subplot(325);
plot(t,V_hat,'b-',t,data.V_spm,'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('V','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])

subplot(326);
plot(t,nLi_hat,'b-',t, ones(NT)*p.n_Li_s,'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('nLi','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])

% Plot results
fs = 15;    % Font Size for plots
figure(2); clf;
set(gcf,'Position',[99,29,300,300]);

subplot(321);
plot(t(end-4250:end),csp_hat(end-4250:end),'b-',t(end-4250:end),data.css_p(end-4250:end),'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{s,p}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([t(end-4250), t(end)])

subplot(322);
plot(t(end-4250:end),csn_hat(end-4250:end),'b-',t(end-4250:end),data.css_n(end-4250:end),'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{s,n}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([t(end-4250), t(end)])

subplot(323);
plot(t(end-4250:end),cs_ave_p_hat(end-4250:end),'b-',t(end-4250:end),data.cs_ave_p(end-4250:end),'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{ave,p}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([t(end-4250), t(end)])

subplot(324);
plot(t(end-4250:end),cs_ave_n_hat(end-4250:end),'b-',t(end-4250:end),data.cs_ave_n(end-4250:end),'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('C_{ave,n}','FontSize',fs);
set(gca,'FontSize',fs)
xlim([t(end-4250), t(end)])

subplot(325);
plot(t(end-4250:end),V_hat(end-4250:end),'b-',t(end-4250:end),data.V_spm(end-4250:end),'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('V','FontSize',fs);
set(gca,'FontSize',fs)
xlim([t(end-4250), t(end)])

subplot(326);
plot(t(end-4250:end),nLi_hat(end-4250:end),'b-',t(end-4250:end), ones(4251,1)*p.n_Li_s,'r-','LineWidth',1);
leg_css = {'Estimated';'True'};
legend(leg_css,'FontSize',fs,'Interpreter','latex', 'Location','southwest')
ylabel('nLi','FontSize',fs);
set(gca,'FontSize',fs)
xlim([t(end-4250), t(end)])


