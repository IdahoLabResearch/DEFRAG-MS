# -*- coding: utf-8 -*-
"""
Knudsen Conversion Plot with Error Bars and Divisors
Author: kristys-sci
Copyright 2026, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""

import pandas as pd
import matplotlib.pyplot as plt

# === USER CONFIGURATION ===
im_files = ["IM_knud1.csv", "IM_knud2.csv", "IM_knud10.csv"]
std_files = ["IM_std_knud1.csv", "IM_std_knud2.csv", "IM_std_knud10.csv"]

pattern = {'nmols_per_pulse': 9.3 + 93.5, 'pulses_per_row': 430}

# Plot all gases or just a subset
plot_all = True  # Set to False to plot only selected gases

# === HARDCODED DIVISORS FOR PLOTTING ===
divisors = {
    '2': 1,     # H2
    '15': 3,    # CH4
    '18': 1,    # H2O
    '26': 1.5,  # C2H4
    '28': 3,    # CO
    '29': 1,    # Propane
    '30': 1.5,  # C2H6
    '40': 1,    # Ar
    '41': 1,    # Propylene
    '44': 3     # CO2
}

colors = {
    '2': 'blue',
    '15': 'violet',
    '18': 'black',
    '26': 'red',
    '28': 'gold',
    '29': 'gray',
    '30': 'orange',
    '40': 'lime',
    '41': 'green',
    '44': 'fuchsia'
}

markers = {
    '2': 'x',
    '15': 'x',
    '18': 'o',
    '26': 'x',
    '28': 'o',
    '29': 'x',
    '30': 'o',
    '40': 'x',
    '41': 'x',
    '44': 'o'
}

labels = {
    '2': 'Hydrogen',
    '15': 'Methane',
    '18': 'Water',
    '26': 'Ethylene',
    '28': 'CO',
    '29': 'Propane',
    '30': 'Ethane',
    '40': 'Argon',
    '41': 'Propylene',
    '44': 'CO₂'
}

# === LOAD AND MERGE FILES ===
dfs = []
for im_file, std_file in zip(im_files, std_files):
    im = pd.read_csv(im_file)
    std = pd.read_csv(std_file)
    df = pd.concat([im, std], axis=1)
    dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

# === BUILD X-AXIS ===
total_nmols = []
cumulative_nmols = 0
for _ in range(len(df)):
    nmols_for_row = pattern['nmols_per_pulse'] * pattern['pulses_per_row']
    cumulative_nmols += nmols_for_row
    total_nmols.append(cumulative_nmols)
df['total_nmols'] = total_nmols

# === CONVERT TO NMOLS ===
for col in divisors.keys():
    nmols = []
    nmols_std = []
    for i in range(len(df)):
        base = df[col][i]
        std = df[f"{col}_std"][i]
        internal = df['40'][i]
        factor = (pattern['nmols_per_pulse'] / internal) * pattern['pulses_per_row']
        nmols.append(base * factor)
        nmols_std.append(std * factor)
    df[f'nmols_{col}'] = nmols
    df[f'nmols_{col}_std'] = nmols_std

# === PLOTTING ===
fig, ax = plt.subplots(dpi=300)

gases_to_plot = list(divisors.keys()) if plot_all else ['2', '18']

for col in gases_to_plot:
    ax.errorbar(
        df['total_nmols'],
        df[f'nmols_{col}'] / divisors[col],
        yerr=df[f'nmols_{col}_std'] / divisors[col],
        fmt=markers.get(col, 'o'),
        linestyle='-',
        color=colors.get(col, 'gray'),
        label=labels.get(col, col),
        capsize=5
    )

ax.set_xlabel('Propane Pulsed (nmols)')
ax.set_ylabel('Product Out (nmols)')
# ax.set_ylim(0, 30000)
plt.legend()
plt.tight_layout()
plt.show()
