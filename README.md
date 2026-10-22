# SPM-EKF

Extended Kalman filter (EKF) for state estimation with a single particle model (SPM) of a lithium-ion cell. The filter estimates the radial lithium concentration profiles in both electrode particles from measured current and voltage. It also uses the total cyclable lithium as a pseudo-measurement.

This is a baseline for physics-informed neural network (PINN) state estimation. It uses the SPM with the LiCoO2 parameter set from Scott Moura's SPMeT, not the LG M50 SPMe used in PINN-SPMe.

## Method

- **Model:** SPM with finite-difference radial discretization (`p.Nr = 15`, so 14 interior nodes per particle). The electrolyte concentration is held constant. The solid diffusion matrices come from `spm_plant_obs_mats_YK_modified.m`, and the voltage comes from `nonlinearSPMOutputVoltage_EKF.m` (OCP, Butler-Volmer overpotentials and film resistance).
- **Filter:** continuous-time EKF (Kalman-Bucy form) in `ode_ekf_spmfdm.m`, integrated with `ode45` over the measured time vector. The filter states are `[c_p (Nr-1); c_n (Nr-1); vec(Sigma)]`, so the Riccati equation is integrated together with the state estimate.
- **Measurements:**
  1. Terminal voltage. The Jacobian `H` is built analytically from the OCP slope and the Butler-Volmer derivative with respect to surface concentration (`exch_cur_dens_for_ode.m`, `refPotential*.m`).
  2. Total solid lithium `n_Li_s` as a pseudo-measurement, using a known or previously identified value. `Run_EKF.m` uses `p.n_Li_s = 2.4655`, the value identified by PSO in the companion parameter-identification repo SPM-OP (true value 2.50).
- **Noise model:** process noise `W = L*L'`, derived from 3% current noise propagated into the surface-node dynamics. Measurement noise is `N = diag((0.03*std_V)^2, 0.05^2)`.
- **Initialization:** both particles start uniform at the concentrations from `init_cs(p, 4.0)` (a 4.0 V guess), with `Sigma_0 = I`.

## Usage

Open MATLAB in this folder and run:

```matlab
Run_EKF
```

`Run_EKF.m` loads `data/SPM_FDM_nr100_simulation_US06_Extended.csv` and adds 3% Gaussian noise to current and voltage. It then runs the filter and plots estimated against true surface concentrations, average concentrations, voltage and total lithium, over the full profile (figure 1) and the last 4,250 s (figure 2).

To use another profile, change the `readtable` path in `Run_EKF.m`. The noise level (`noise`), the identified `n_Li_s` and the initial voltage passed to `init_cs` are also set at the top of `Run_EKF.m`.

## Files

| File | Purpose |
|---|---|
| `Run_EKF.m` | Entry point: parameters, data, noise, initialization, filter run, plots |
| `ode_ekf_spmfdm.m` | EKF right-hand side: state estimate, output Jacobian, Kalman gain, Riccati equation |
| `params_setting.m` | SPM parameters (LiCoO2 cell) |
| `params_setting_2.m` | The same parameter set with selected parameters scaled by 0.8 or 1.2, representing a "different" (degraded) cell |
| `spm_plant_obs_mats_YK_modified.m` | Solid-diffusion state-space matrices; returns both boundary concentrations (`spm_plant_obs_mats.m` is the original, which returns the surface value only) |
| `nonlinearSPMOutputVoltage_EKF.m` | Voltage for a constant electrolyte (used inside the filter) |
| `nonlinearSPMOutputVoltage_Scott.m` | Original SPMeT voltage function (used for plotting) |
| `exch_cur_dens.m`, `exch_cur_dens_for_ode.m` | Exchange current density; the `_for_ode` version also returns its derivative |
| `refPotentialAnode.m`, `refPotentialCathode.m` | Open-circuit potentials and their slopes |
| `init_cs.m` | Initial particle concentrations from a voltage |
| `calc_c_ave.m` | Volume-averaged particle concentration |
| `SPM Derivation.PDF`, `SPMFDM Kalman Filter.pdf` | Derivations of the SPM finite-difference model and of this EKF |
| `theta_calc.xlsx` | Hand calculation of stoichiometry windows and capacity for the parameter set |

## Data

- `data/SPM_FDM_nr{N}_simulation_{profile}.csv`: synthetic SPM simulations with radial resolution `N`. Columns (no header): time [s], current [A], voltage [V], surface and average anode concentration, surface and average cathode concentration [mol/m^3]. The `nr100` files serve as ground truth for the `Nr = 15` filter. Positive current is discharge.
- `data/11_*_SP20-2_*_80SOC.*`: CALCE drive-cycle tests (DST, FUDS, US06, BJDST at 80% SoC) on a different cell, kept as current-profile sources. `Run_EKF.m` does not use them.
- `data/UDDSx2/UDDSx2_DFN_Co_1sec.mat`, `data/Sig_last.mat`: a DFN simulation and a saved covariance from earlier experiments.

## Known issues

Check these before reusing or porting the filter:

- In `ode_ekf_spmfdm.m`, `n_Li_s_n_hat` sums `c_n_hat(k) * H_p(2,k)`, where `H_n(2,k)` (anode volume) looks intended.
- The lines setting `H_p(2,Nr-1)` and `H_n(2,Nr-1)` add to the `Nr-2` entry rather than the `Nr-1` entry.

Both affect only the lithium-inventory pseudo-measurement, not the voltage channel.

## Related repositories

- [SPM-OP](https://github.com/yuichikajiura/SPM-OP): PSO/GA parameter identification for the same SPM.
- [PINN-SPMe](https://github.com/yuichikajiura/PINN-SPMe) and [LGM50](https://github.com/ESCL-OU/LGM50): PINN state and SEI estimation for an SPMe with degradation, and the MATLAB simulator that generates its data.

## Credits

The model functions (`spm_plant_obs_mats.m`, `nonlinearSPMOutputVoltage_Scott.m`, `init_cs.m`, `refPotential*.m`, `exch_cur_dens.m`) are adapted from Scott Moura's SPMeT, which is licensed under GPL-3.0.
