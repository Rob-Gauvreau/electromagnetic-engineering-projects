# ELG3106 — Design Project: Broadband Anti-Reflection Coatings for Solar Power (Fall 2025)

This folder contains my final design project for **ELG3106 (Electromagnetic Engineering)**. The project evaluates how **broadband anti-reflection (AR) coatings** on a solar cell impact overall **electrical power production**, comparing **double-layer** and **triple-layer** coating designs.

## What this project covers
- Modeling wavelength-dependent **reflection and transmission** of multilayer dielectric stacks
- Comparing **2-layer vs 3-layer** AR coating behavior and tradeoffs (bandwidth vs complexity)
- Using a spectrum-weighted approach to estimate **delivered power** to the solar cell
- Exploring how refractive indices and layer thickness choices affect performance across a wide wavelength range

## Project structure
The report is organized into four main parts:
- **Part 1:** Baseline behavior and validation at a design wavelength (single-wavelength focus)
- **Part 2:** Full-spectrum reflectivity/transmission behavior and spectrum-weighted power calculation
- **Part 3:** Extending the model to a **triple-layer** coating and evaluating improved bandwidth performance
- **Part 4:** Numerical sweeps/optimization to maximize transmitted power over selected wavelength ranges

## What’s included
- **Report (PDF):** Full design report with results, discussion, and conclusions
- **Flowcharts:** High-level workflows for the simulations and parameter sweeps (appendix)
- **Code:** MATLAB implementations used to compute spectra, sweep parameters, and compare designs (appendix)
- **Simulation outputs:** Plots and screenshots verifying reflectivity trends and power comparisons (appendix)
- **Results + discussion:** Interpretation of why single-wavelength optimization differs from spectrum-wide optimization 

## Key takeaways
- AR coatings can significantly reduce reflection at a target wavelength, but maximizing total power requires **broadband (spectrum-wide) optimization** rather than tuning at a single wavelength.
- Adding a third layer can further lower reflectivity around the design wavelength, but the best overall design depends on the wavelength range used for the power metric.

## Notes
This repository is meant to present a complete design workflow: **theory → implementation → verification → optimization → interpretation**, with emphasis on making design decisions based on performance across the full operating spectrum.
