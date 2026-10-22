# GAA-Baseline-Audit-Supplement

**Journal:** *Micromachines* (MDPI)  
**Article type:** Article

## Supplementary Materials: Data and Analysis Code

**Manuscript:** *Baseline Uncertainty Limits Mechanistic Interpretation of Power-Law Degradation Exponents in GAA Transistors*

**Authors:** Liron Cohen, Lidor Bazon, and Joseph B. Bernstein  
Department of Electrical and Electronics Engineering, Ariel University, Ariel 40700, Israel

---

## 1. Overview

This repository contains the digitized time-domain degradation datasets and the Base-MATLAB analysis script used to generate the statistical audits and tables reported in the manuscript.

## 2. Included Files

| File | Description |
|------|-------------|
| `P03_Primary_GAA_Dataset.csv` | Digitized time-series degradation data for the primary GAA nanosheet dataset (P03, Ref. [9]), comprising 534 measurement points across 24 constant-stress curves (4 bias groups × 6 temperatures). |
| `External_Nanowire_Datasets.csv` | Digitized time-series degradation data for the external nanowire audit datasets: EXT-01 (Ref. [12], 18 curves), EXT-02 (Ref. [11], 5 curves), EXT-03 (Ref. [10], 6 curves, 3 time points per curve), and EXT-04 (Ref. [13], 3 curves). |
| `Run_GAA_Baseline_Audit_Supplement.m` | Self-contained Base-MATLAB script (no toolboxes required). |

> **Note on data sources:** The bracketed reference numbers ([9]–[13]) correspond to the reference list of the manuscript. Full citations for the original source publications of each dataset are provided in **Table 1** (primary P03 dataset) and **Table 4** (external nanowire datasets) of the manuscript.

### The analysis script executes:

- **(a)** Full-curve zero-curvature baseline reconstruction and root-multiplicity auditing on the primary P03 dataset (Table 2, Sections 3.1–3.3).
- **(b)** Single-point leave-one-out (LOO) stability auditing (Section 3.2).
- **(c)** Synthetic constant-exponent power-law verification benchmark (Section 2.3).
- **(d)** Within-curve 4-block evolution, leave-one-block-out (LOBO) cross-fitting, and matched residual-phase-permutation null testing (Section 3.4).
- **(e)** Inverse thermal energy regressions, baseline-exponent coupling, and threshold-time perturbation / LOTO auditing (Table 3, Section 3.5).
- **(f)** External nanowire zero-curvature root multiplicity, log-log R² comparison, and LOO jackknife stability auditing (Table 4, Section 3.6).

## 3. Execution Instructions

1. Download all files and place them in the same working folder in MATLAB (R2020b or newer).
2. Run the script from the MATLAB Command Window:

```matlab
   >> Run_GAA_Baseline_Audit_Supplement
```

3. The script prints all manuscript summary tables directly to the Command Window and exports the processed tables to the subfolder `./Supplement_Audit_Outputs/`.

## 4. Requirements

- MATLAB R2020b or newer
- No additional toolboxes required

## 5. Citation

If you use these data or code, please cite the associated manuscript:

> Cohen, L.; Bazon, L.; Bernstein, J.B. Baseline Uncertainty Limits Mechanistic Interpretation of Power-Law Degradation Exponents in GAA Transistors. *Micromachines* **2026**, *Volume*, Article number. https://doi.org/xxxxx

*(Volume, article number and DOI will be updated upon publication.)*

When using the digitized datasets, please also cite the original source publications listed in Tables 1 and 4 of the manuscript.
