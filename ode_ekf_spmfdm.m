function [xdot] = ode_ekf_spmfdm(t,x, Ap, An, Bp, Bn, input_data, W, N, p)

Nr = p.Nr;

% Parse States
c_p_start = 1;
c_p_end = c_p_start + Nr-2;
c_n_start = c_p_end + 1;
c_n_end = c_n_start + Nr-2;
Sig_start = c_n_end + 1;
Sig_end = Sig_start + (2*(Nr-1))^2 - 1;
c_p_hat = x(c_p_start:c_p_end);
c_n_hat = x(c_n_start:c_n_end);
Sig = reshape(x(Sig_start:Sig_end), [2*(Nr-1), 2*(Nr-1)]);

% Parse and interpolate input signal data
%   This is a coding trick for time-varying input data.
it = input_data(:,1);
iI = input_data(:,2);
iy_m = input_data(:,3);
I = interp1(it,iI,t);
y_m = interp1(it,iy_m,t);

%% Compute H(t)
H_p = zeros(2, Nr-1);
H_n = zeros(2, Nr-1);

% Differentiate overpotential by concentration
csp = 4/3 * c_p_hat(Nr-1) - 1/3 * c_p_hat(Nr-2) ...
    +(2*p.delta_r_p * p.R_s_p)/(3*p.D_s_p * p.Faraday * p.a_s_p * p.Area * p.L_p)*I;
csn = 4/3 * c_n_hat(Nr-1) - 1/3 * c_n_hat(Nr-2) ...
    -(2*p.delta_r_n * p.R_s_n)/(3*p.D_s_n * p.Faraday * p.a_s_n * p.Area * p.L_n)*I;
csp = min(csp, p.c_s_p_max-1);
csn = min(csn, p.c_s_n_max-1);
csp = max(1, csp);
csn = max(1, csn);
y_hat = nonlinearSPMOutputVoltage_EKF(p,csn, csp, p.c_e, I);

[i0n,i0p,di0n_dcsn, di0p_dcsp] = exch_cur_dens_for_ode(p, csn,csp,p.c_e);

RTaF=(p.R*p.T_amp)/(p.alph*p.Faraday);
Xn = I/(2*p.a_s_n*p.Area*p.L_n*i0n);
Xp = -I/(2*p.a_s_p*p.Area*p.L_p*i0p);

detap_di0p = RTaF*I/(2*p.a_s_p*p.Area*p.L_p*(Xp^2+1)^(1/2));
detan_di0n = - RTaF*I/(2*p.a_s_n*p.Area*p.L_n*(Xn^2+1)^(1/2));

detap_dcsp = detap_di0p * di0p_dcsp;
detan_dcsn = detan_di0n * di0n_dcsn;

% Differentiate open circuit potential by concentration
theta_p = csp / p.c_s_p_max;
theta_n = csn / p.c_s_n_max;

[~, dUp_dcsp] = refPotentialCathode(p, theta_p);
[~, dUn_dcsn] = refPotentialAnode(p, theta_n);

% Differentiate voltage equation by concentration
dg_dcsp = detap_dcsp + dUp_dcsp;
dg_dcsn = - detan_dcsn - dUn_dcsn;

H_p(1,Nr-2) = -1/3 * dg_dcsp;
H_n(1,Nr-2) = -1/3 * dg_dcsn;
H_p(1,Nr-1) = 4/3 * dg_dcsp;
H_n(1,Nr-1) = 4/3 * dg_dcsn;

%Compute d(nLi)/d(c_k) at each radical discretization point k for all k={1:Nr-1} 
Vol_n = p.Area * p.L_n * p.epsilon_s_n;
Vol_p = p.Area * p.L_p * p.epsilon_s_p;
for k = 1:Nr-1
    H_p(2,k) = (1/Nr)^3*(3*k^2+1/4); 
    H_n(2,k) = (1/Nr)^3*(3*k^2+1/4);
end
H_p(2,1) = H_p(2,1) + (4/3)*(2/Nr)^3;
H_n(2,1) = H_n(2,1) + (4/3)*(2/Nr)^3;
H_p(2,2) = H_p(2,2) - (1/3)*(2/Nr)^3;
H_n(2,2) = H_n(2,2) - (1/3)*(2/Nr)^3;
H_p(2,Nr-2) = H_p(2,Nr-2) - (1/3)*(1/Nr)^3*(3*Nr^2+1/4);
H_n(2,Nr-2) = H_n(2,Nr-2) - (1/3)*(1/Nr)^3*(3*Nr^2+1/4);
H_p(2,Nr-1) = H_p(2,Nr-2) + (4/3)*(1/Nr)^3*(3*Nr^2+1/4);
H_n(2,Nr-1) = H_n(2,Nr-2) + (4/3)*(1/Nr)^3*(3*Nr^2+1/4);
H_p(2, :) = H_p(2, :) * Vol_p;
H_n(2, :) = H_n(2, :) * Vol_n;
H = [H_p, H_n];

n_Li_s_p_hat = 0;
n_Li_s_n_hat = 0;
for k = 1:Nr-1
    n_Li_s_p_hat = n_Li_s_p_hat + c_p_hat(k) * H_p(2,k); 
    n_Li_s_n_hat = n_Li_s_n_hat + c_n_hat(k) * H_n(2,k); 
end
n_Li_s_hat = n_Li_s_p_hat + n_Li_s_n_hat;

% Compute Kalman Gain (Look at Chapter 3, Section 4)
L = Sig * H' * inv(N);
L_p = L(1:Nr-1, :);
L_n = L(Nr:2*(Nr-1), :);

% Kalman Filter equations
c_p_hat_dot = Ap * c_p_hat + Bp * I + L_p * [y_m - y_hat; p.n_Li_s - n_Li_s_hat];
c_n_hat_dot = An * c_n_hat + Bn * I + L_n * [y_m - y_hat; p.n_Li_s - n_Li_s_hat];

A = [Ap, zeros(Nr-1, Nr-1); zeros(Nr-1, Nr-1), An];

% Riccati Equation for Sigma
Sig_dot = Sig * A' + A * Sig + W - Sig * H' * inv(N) * H *Sig;

% Concatenate LHS
xdot = [c_p_hat_dot; c_n_hat_dot; reshape(Sig_dot, [2*(Nr-1)*2*(Nr-1), 1])];