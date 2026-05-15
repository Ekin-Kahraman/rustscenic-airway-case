# Security policy

## Supported versions

Security and correctness fixes are made on the default branch.

## Reporting a vulnerability

This repository processes local AnnData, parquet, CSV and PNG outputs. If you find unsafe path handling, dependency supply-chain risk, or unintended file writes, do not open a public issue.

Email: evk23umu@uea.ac.uk
Subject: `rustscenic-airway-case security: <short description>`

Please include:
- Reproduction steps
- Python version and operating system
- Input file type and source
- Expected impact

I aim to respond within 7 days.

## Data handling

Do not attach private patient-level data or restricted clinical metadata to public issues. The source case study uses the public Ziegler et al. atlas; users should supply their own local copy via `ZIEGLER_H5AD`.
