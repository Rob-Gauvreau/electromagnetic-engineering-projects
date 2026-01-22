# ELG3106 — Assignment 3: Microstrip Line Impedance Modeling (Fall 2025)

This folder contains my work for **ELG3106 (Electromagnetic Engineering) — Assignment 3**, focused on modeling a **lossless microstrip transmission line** and evaluating how well common curve-fit formulas match simulation results for different geometries and substrate permittivities. :contentReference[oaicite:0]{index=0}

## What this assignment covers 
- Forward modeling of **characteristic impedance vs. geometry** by sweeping the width-to-height ratio (w/h)
- Comparing calculated impedance to a **virtual microstrip simulation tool** at specific geometry checkpoints
- Inverse design: estimating required **w/h from a target impedance** using two approximation regimes
- Quantifying model accuracy using **relative error plots** and identifying where each approximation is valid

## What’s included
- **Report (PDF):** Full write-up with workflow, results, and discussion
- **Flow diagrams:** Computation pipeline for the forward and inverse comparisons
- **Code:** Script(s) used to generate calculated values and relative-difference plots
- **Virtual simulation screenshots:** Verification points for multiple w/h values and permittivities
- **Results + discussion:** Interpretation of error trends (including where errors grow and why)
