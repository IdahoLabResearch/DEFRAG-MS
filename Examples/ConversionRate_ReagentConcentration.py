# -*- coding: utf-8 -*-
"""
RC Conversion Analysis with Knudsen Diffusion and Active Site Estimation
Author: kristy-sci 
Copyright 2026, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from numpy import sqrt, linspace, interp, nan_to_num
from matplotlib.ticker import FuncFormatter
from matplotlib.pyplot import cm

# === USER CONFIGURATION ===
path = R"C:\Users\KRISST\OneDrive - Idaho National Laboratory\Desktop\Cedar\Ga_May\redo_for_paper\h2\code_git"  # Folder containing flux_*.csv files
cat_ratio = 0.5
r = 2 / 1000  # m
d = r * 2
nmol = 10
mol = 10 / (10**9)
eb = 0.4
L = 39.2 / 1000  # m
cross_sectional_area = np.pi * d**2
units_c = mol / (eb * cross_sectional_area * L)

# === FILE HANDLING ===
files = sorted([f for f in os.listdir(path) if f.startswith("flux_knud") and f.endswith(".csv")],
               key=lambda f: int(f.split("knud")[1].split(".")[0]))

colors = cm.plasma(np.linspace(0, 1, len(files)))
tick_labels = [0, 100000, 200000, 300000, 400000, 500000]
fit_parameters = []
output_data = pd.DataFrame()

fig, ax = plt.subplots(figsize=(10, 6))

for idx, file in enumerate(files):
    df = pd.read_csv(os.path.join(path, file))
    df['time'] = np.arange(0, 5.997, 0.001)

    # === Average AMU columns ===
    base_columns = ['2', '15', '18', '26', '28', '29', '30', '40', '41', '44']
    averaged = {col: df[[c for c in df.columns if c.startswith(col)]].mean(axis=1) for col in base_columns}
    averaged_df = pd.DataFrame(averaged)
    averaged_df['time'] = df['time']

    # === Background subtraction ===
    bkg = averaged_df.iloc[1:2].mean()
    averaged_df = averaged_df - bkg.values
    averaged_df['time'] = df['time']
    averaged_df = averaged_df.drop(index=[0, 1, 2, 3])

    # === Inert and Reactant Flux ===
    inert_flux = averaged_df['40'].to_numpy()
    flux = averaged_df['29'].to_numpy()
    times = averaged_df['time'].to_numpy()

    # === Graham's Law Correction ===
    grahams_constant = sqrt(40 / 44)
    stop_point = int(round(len(times) / grahams_constant))
    original_m01 = integrate.simps(inert_flux, times)

    if grahams_constant >= 1:
        tail = inert_flux[stop_point:]
        inert_flux = np.concatenate([
            interp(linspace(times[0], times[-1], stop_point), times, inert_flux),
            tail
        ])
    else:
        inert_flux = interp(times * grahams_constant, times, inert_flux)

    inert_flux = nan_to_num(inert_flux)
    inert_flux *= original_m01 / integrate.simps(inert_flux, times)
    flux_diff = inert_flux - flux

    # === Rate Calculation ===
    time_scalar = -(1 - cat_ratio**2) * 1.5
    time_multiplier = np.concatenate([[0], times[1:]**time_scalar])
    rate = flux_diff * time_multiplier
    rate[0] = 0
    rate *= (integrate.simps(flux_diff, times) / integrate.simps(rate, times)) * (-time_scalar)

    # === Concentration Calculation ===
    time_scalar_c = -(1 - cat_ratio**2) / 6
    time_multiplier_c = np.concatenate([[0], times[1:]**time_scalar_c])
    concentration = flux * time_multiplier_c
    concentration *= integrate.simps(flux, times) / integrate.simps(concentration, times)

    # === Unit Conversion ===
    tp = times[np.argmax(flux_diff)]
    D = (1 / 6) * (eb * L**2) / tp
    units_r = nmol * D / (eb * L**2)
    rate *= units_r
    concentration *= units_c * 100

    # === Fit and Plot ===
    fit_range = 11
    x_fit = concentration[:fit_range]
    y_fit = rate[:fit_range]
    slope, intercept = np.polyfit(x_fit, y_fit, 1)
    delta_x = np.max(x_fit) - np.min(x_fit)
    delta_y = np.max(y_fit) - np.min(y_fit)
    fit_parameters.append((idx + 1, slope, intercept, delta_x, delta_y))

    ax.plot(concentration, rate, marker='x', color=colors[idx])
    ax.plot(concentration, slope * concentration + intercept, linestyle='-', color='r')

    # === Save individual data ===
    data = pd.DataFrame({'rate': rate, 'concentration': concentration})
    data.to_csv(os.path.join(path, f"RC_knud{idx + 1}.csv"), index=False)
    output_data = pd.concat([output_data, data], axis=1)

# === Save combined output ===
output_data.to_csv(os.path.join(path, "concatenated_output.csv"), index=False)

# === Colorbar ===
sm = plt.cm.ScalarMappable(cmap=cm.plasma, norm=plt.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = plt.colorbar(sm, ticks=np.linspace(0, 1, 6))
cbar.ax.set_yticklabels(tick_labels)
cbar.ax.set_title('Propane Pulsed [nmols]')
cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{tick_labels[pos]:.1e}'))

# === Labels ===
plt.xlabel('Propane Concentration [$10^2$ mol/$m^3$]')
plt.ylabel('Propane Conversion Rate [nmol/s $g_{cat}$]')
plt.tight_layout()
plt.savefig(os.path.join(path, "combined_plot_with_colorbar.png"))
plt.show()

# === Print Fit Table ===
print("\nTable of Fit Parameters:")
print(f"{'Fit #':<8}{'Slope':<12}{'Intercept':<12}{'Δx':<12}{'Δy':<12}")
for fit in fit_parameters:
    fit_num, slope, intercept, dx, dy = fit
    print(f"{fit_num:<8}{slope:<12.4f}{intercept:<12.4f}{dx:<12.4f}{dy:<12.4f}")

# === Active Site Estimation ===
active_sites = []
active_sites_std = []

for i in range(len(fit_parameters)):
    rate_col = output_data.iloc[:, 2 * i] * 1e-9  # mol/s/g_cat
    conc_col = output_data.iloc[:, 2 * i + 1]     # mol/m^3
    valid = (rate_col > 0) & (conc_col > 0)

    if valid.sum() < 11:
        N = sigma_N = np.nan
    else:
        R_vals = rate_col[valid].iloc[:11]
        C_vals = conc_col[valid].iloc[:11]
        R0 = R_vals.mean()
        C0 = C_vals.mean()
        sigma_R0 = R_vals.std(ddof=1)
        sigma_C0 = C_vals.std(ddof=1)

        if i == 0:
            k0 = fit_parameters[0][1]
            sigma_k0 = 0.05 * k0  # assume 5% error
        N = R0 / (k0 * C0) if k0 != 0 else np.nan
        rel_error = (sigma_R0 / R0)**2 + (sigma_k0 / k0)**2 + (sigma_C0 / C0)**2
        sigma_N = N * np.sqrt(rel_error)

    active_sites.append(N)
    active_sites_std.append(sigma_N)

# === Print Active Site Table ===
print("\nUpdated Table with Active Site Estimates and Standard Deviations:")
print(f"{'Fit #':<6}{'Slope':<10}{'Intercept':<12}{'Δx':<10}{'Δy':<10}{'Active Sites (mol/g_cat)':<28}{'Std Dev':<10}")
for i, fit in enumerate(fit_parameters):
    slope, intercept, dx, dy = fit[1], fit[2], fit[3], fit[4]
    N = active_sites[i]
    sigma_N = active_sites_std[i]
    print(f"{i+1:<6}{slope:<10.4f}{intercept:<12.4f}{dx:<10.4f}{dy:<10.4f}{N:<28.4e}{sigma_N:<10.2e}")
