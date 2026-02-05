
Copyright 2026, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
# DEFRAG-MS
DEFRAG-MS is a preprocessing tool for TAP mass spectrometry data that corrects instrument gain and background and defragments overlapping ion signals to recover accurate, time-resolved gas fluxes. It supports kinetic analysis and catalyst evaluation for energy and industrial reaction systems.


# Mass Spec Transient Analysis Toolkit - DEFRAG-MS

This repository contains tools for processing and analyzing transient mass spectrometry (MS) data from TAP (Temporal Analysis of Products) experiments. It includes a robust defragmentation pipeline which is the main purpose of this site, but also includes rate/concentration analysis scripts based on example data sets. 

## Background

TAP experiments provide time-resolved MS data for catalytic reactions. However, overlapping fragmentation patterns make quantification difficult. This toolkit implements a defragmentation method using calibrated fragmentation matrices and non-negative least squares (NNLS) regression to recover true species fluxes.

## 📁 Repository Structure

mass-spec-transient-analysis/
├── preprocess_mass_spec.py         # Core preprocessing script
├── data/                           # Example input/output files
├── examples/                       # Analysis scripts
├── docs/                           # SI and manuscript

## 🚀 Getting Started

### 1. Clone the repository

git clone https://github.com/yourusername/mass-spec-transient-analysis.git
cd mass-spec-transient-analysis

### 2. Install dependencies

pip install -r requirements.txt

### Calibration Files

This repository includes example calibration files:

- `gc.csv`: Gain correction factors
- `defrag_key.csv`: Defragmentation matrix

These are specific to the example data provided. For your own experiments, you must generate your own calibration files as described in the manuscript and Supporting Information (see `docs/SI.pdf`).

### 3. Run preprocessing

I have scripted this code for use in Spyder, but this should work elsewhere. 

python preprocess_mass_spec.py --input data/example_flux/knud1.xlsx --mode average --bkg --plot --integrate

This will generate:
- flux_knud{label}.csv: Time-resolved fluxes
- IM_knud{label}.csv: Integrated moments
- IM_std_knud{label}.csv: Standard deviations
where {label} is a file identifier specific to your experiment. 

## 📊 Example Analyses

All examples are in the `examples/` folder:

- experiment_specific_moment_calculations.py: Converts moment/std (IM) to nmols with error bars.
- ConversionRate_ReagentConcentration.py: Calculates conversion rate vs. reagent concentration using G procedure. 
- FormationRate_ReagentConcentration.py: Calculates propylene formation rate vs reagent concentration, with Savitzky-Golay smoothing.

## 📄 Documentation

- docs/main_manuscript.pdf: Main manuscript describing the utility of the steps in this script. 
- docs/SI.pdf: Supporting Information with calibration, defragmentation, and analysis details

## 📚 References

This work is based on the methodology described in:

> Kristy, S. et al. *Transient Catalytic Reaction Analysis Through Signal Defragmentation*. (2026)

VTAP example data is provided via TAP simulations via the framework developed by Wang et. al. Wang, Shengguang, et al. "A Simulation Framework for Understanding Transport and Kinetics in Transient Reactor Experiments." (2025).

## 🛠 Requirements

- Python 3.8+
- pandas
- numpy
- scipy
- matplotlib
- openpyxl

Install with:

pip install -r requirements.txt

## 📜 License

MIT License — You are free to use, modify, and distribute this code.

## 🙏 Acknowledgments

This work was supported by the U.S. Department of Energy (DOE), Office of Energy Efficiency and Renewable Energy (EERE), Industrial Efficiency and Decarbonization Office (IEDO), under contract DE-AC07-05ID14517.

Special thanks to Ross Kunz for advising on defragmentation workflows. Thanks to Shengguang Wang for providing example VTAP data. 
