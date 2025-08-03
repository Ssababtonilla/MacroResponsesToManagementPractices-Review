# MacroResponsesToManagementPractices-Review

![R ≥ 4.3.2](https://img.shields.io/badge/R-%3E%3D4.3.2-blue) ![Codespaces Ready](https://img.shields.io/badge/Codespaces-ready-orange) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

A reproducible meta-analysis of stream macroinvertebrate responses to best management practices (BMPs) across regions, land uses, practice types, and assessment metrics in the Chesapeake Bay watershed. All data, code, and outputs supporting:

> **Sabat‑Bonilla, S., et al.** (2025). *Meta‑analysis of Stream Macroinvertebrate Responses to Management Practices in Selected Regions of the Chesapeake Bay Watershed Vary with Region, Land Use, Practice Type and Assessment Metrics.*

## Repository Structure

```
├── .devcontainer/             # Dev container config for GitHub Codespaces & local dev
│   └── devcontainer.json      # Defines R image, Git, LFS, and recommended VSCode extensions
├── data/                      # Raw and processed datasets
│   └── raw/
│       └── Regional BMP Database - Data.csv
│   └── processed/             # (auto-generated)
├── src/                       # R analysis scripts
│   └── sabat_bonilla_bmp_review.R
├── outputs/                   # Tables, figures, model objects (auto-generated)
├── .github/workflows/ci.yml   # (optional) automated checks
├── LICENSE                    # MIT license for code
├── LICENSE-data.md            # CC-BY-4.0 license for data
├── CITATION.cff               # Citation metadata for repository
└── README.md                  # This file
```

## Quick Start

1. **Clone & restore**

   ```bash
   git clone https://github.com/Ssabatbonilla/MacroResponsesToManagementPractices-Review.git
   cd MacroResponsesToManagementPractices-Review
   ```
2. **Launch in Codespaces**

   * Click **Code → Open with Codespaces → New codespace**
   * Wait for R, Git, and VSCode extensions to install automatically
3. **Local dev with Dev Container**

   ```bash
   code .
   # When prompted, “Reopen in Container” to spin up Docker-based R environment
   ```
4. **Run the analysis**

   * In an R terminal (`⇧⌘P` → `R: Create R Terminal`):

     ```r
     source("src/sabat_bonilla_bmp_review.R")
     ```
   * Or headless:

     ```bash
     Rscript src/sabat_bonilla_bmp_review.R
     ```

All generated CSVs and PNGs land in `outputs/`.

## Environment & Reproducibility

* **R packages** managed via `renv`; restore with:

  ```r
  renv::restore()
  ```
* **Git LFS** tracks large CSVs (ensure `git lfs install` locally)

## .gitignore

```gitignore
# R
.Rhistory
.RData
.Rproj.user/

# VSCode
.vscode/

env/
outputs/
```

## License & Citation

* **Code**: MIT License – see [LICENSE](LICENSE)
* **Data**: CC-BY‑4.0 – see [LICENSE-data.md](LICENSE-data.md)
* **To cite**: refer to [CITATION.cff](CITATION.cff)

---

Maintained by **Sergi Sabat‑Bonilla** ([ssabatbonilla@vt.edu](mailto:ssabatbonilla@vt.edu)). Feel free to open issues or pull requests!
