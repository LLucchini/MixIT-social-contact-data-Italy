# MixIT post-pandemic social contact data for Italy — preprocessing & contact matrices

Preprocessing and analysis code for the **MixIT_2022_2023** social contact data collection in Italy, with notebooks to:
- clean and harmonize respondent/contact records,
- fit descriptive/statistical models,
- derive **age-stratified contact matrices** (including handling of indirect contacts),
- assess sample representativeness and produce summary tables.

This repository is primarily **analysis code**; access to the underlying raw survey microdata may be restricted depending on the project’s data governance.

---
## Reference
This code was developed within the scope of the work titled "Post-pandemic social contacts in Italy: implications for social distancing measures on in-person school and work attendance" by Lucchini et al (2026). For more details on the analysis and results, please refer to and cite the ArXiv preprint: https://arxiv.org/abs/2412.18549.

---

## Context

Contact-diary / contact-survey data are commonly used to build **social mixing matrices** for infectious-disease transmission models. For a widely used open ecosystem and schema conventions, see **socialcontactdata.org** (data organization and standard “participant/contact/household/survey-day/time-use/dictionary” tables).  
Note that the data requirements to run most of this pipeline require a standard schema such as the one proposed by the ``socialcontactdata`` community: https://socialcontactdata.org/data/

---

## Repository structure

Top-level files and directories (from the repo root):
- `MixIT/` — project-specific materials (raw/derived data and/or configuration, depending on your local setup)
- `preprocessing/` — preprocessing utilities/scripts (project-specific)
- `tables/` — exported tables for main and Supplementary Information
- Notebooks (intended to be run in order):
  - `1-statistical_model.ipynb`
  - `2-statistical_model_results.ipynb`
  - `3-matrices_contact_data_formatting.ipynb`
  - `4-exploding_indirect_contacts.ipynb`
  - `5-matrices_computation.ipynb`
  - `6-data_representativeness.ipynb`
- `7-tables.R` — table generation (R)
- `utils.py` — shared plotting + matrix utilities (Python)
- `LICENSE` — [MIT](https://github.com/LLucchini/MixIT-social-contact-data-Italy)

---

## The pipeline output
- cleaned **respondent** and **contact** records,
- harmonized categorical variables (e.g., contact location mapping, age groups),
- **contact matrices** (by age group, setting, and attendance status),
- bootstrap uncertainty for matrices (if enabled),
- representativeness comparisons vs. national benchmarks,
- publication-ready tables.

The provided `utils.py` includes functions for matrix symmetrization and bootstrapping (commonly used when producing reciprocal contact rates), plus plotting helpers.

---

## Data

This repo does **not** expose a public raw dataset at the root. You will generally need:
- respondent-level survey data (demographics, weights if applicable, survey-day metadata),
- contact-level diary entries (contact age/age group, setting, physical vs. conversational, etc.),
- optional household and time-use components (if collected).

If you intend to align with the SocialContactData conventions, the target relational tables are described here: [https://socialcontactdata.org/data/](https://socialcontactdata.org/data/?)

### Notes on open data
Note that the statistical model requires data that are not publicly shared and, if needed, should be requested from the corresponding author or data managers.
Notebook 3: transform MixIT data into socialcontactdata formatting style. Notebook 4-6 can be run with standard structure data available on [Zenodo](https://zenodo.org/records/17579537).

---

## Environment & installation

### R and Python package versions
The code was developed using R version 4.5.1 (2025-06-13) and Python version 3.12.11.

