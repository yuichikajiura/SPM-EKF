%% Params for Electrochemical Model
%   Created Mar 22, 2016 by Saehong Park
%   For SPM

%% Geometric Params
% Thickness of each layer
p.L_n = 100e-6;     % Thickness of negative electrode [m]
p.L_s = 25e-6;     % Thickness of separator [m]
p.L_p = 100e-6;     % Thickness of positive electrode [m]

L_ccn = 25e-6;    % Thickness of negative current collector [m]
L_ccp = 25e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 10e-6;   % Radius of solid particles in negative electrode [m]
p.R_s_p = 10e-6;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.6;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.5;      % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 0.3;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 1.0;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.3;   % Volume fraction in electrolyte for pos. electrode

epsilon_f_n = 0.1;  % Volume fraction of filler in neg. electrode
epsilon_f_p = 0.2;  % Volume fraction of filler in pos. electrode

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]


%% Transport Params

% Diffusion coefficient in solid
p.D_s_n = 3.9e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p = 1e-13;  % Diffusion coeff for solid in pos. electrode, [m^2/s]

% Diffusion coefficient in electrolyte
p.D_e = 5.34e-10*exp(-0.65);    % Diffusion coeff for electrolyte, [m^2/s]


p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 10;    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% p.sig_eff_n = p.sig_n * p.epsilon_s_n^p.brug;    % Eff. conductivity in neg. electrode, [1/Ohms*m]
% p.sig_eff_p = p.sig_p * p.epsilon_s_p^p.brug;    % Eff. conductivity in pos. electrode, [1/Ohms*m]

% Conductivity of electrolyte

% Miscellaneous
p.t_plus = 0.4;       % Transference number
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 1;           % Electrode current collector area [m^2]

%% Kinetic Params

p.R = 8.314472;       % Gas constant, [J/mol-K]
p.alph = 0.5;         % Charge transfer coefficients
% p.R_SEI = 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_n = 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;       % Resistivity of SEI layer, [Ohms*m^2]

% Reaction rates
% p.k_n = 1e-5;  % Reaction rate in neg. electrode, [(A/m^2)*(m^3/mol)^(1+alpha)]
% p.k_p = 3e-7; % Reaction rate in pos. electrode, [(A/m^2)*(m^3/mol)^(1+alpha)]

% Reaction rate constant [m2.5 mol-0.5 s-1]*[C/mol] = [(A/m^2)*(m^3/mol)^(1.5)] from https://github.com/davidhowey/Spectral_li-ion_SPM/blob/master/model_parameters/get_modelData.m
p.k_n = 1.764e-11 * p.Faraday;
p.k_p = 6.667e-11 * p.Faraday;

%% Thermodynamic Params

% Thermal dynamics
p.C_p = 2000;   % Heat capacity, [J/kg-K]
p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0

% Ambient Temperature
p.T_amp = 298; % [K]


%% Concentrations
% % Maxima based on 2.3Ah cell
% p.c_s_n_max = 1.4990e4;   % Max concentration in anode, [mol/m^3]
% p.c_s_p_max = 2.5609e4;    % Max concentration in cathode, [mol/m^3]

% Maxima based on DUALFOIL
p.c_s_n_max = 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]
p.c_s_p_max = 3.6e3 * 247 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.n_Li_s = 2.50;        % Total moles of lithium in solid phase [mol]
p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]

