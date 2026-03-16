# Security Policy

## Supported Versions

Only the current `main` branch and actively developed PR branches are maintained for security updates.

## Dealing with Patient Data
**CRITICAL SECURITY IMPLICATION**: This repository models cardiovascular physiology and may be used with real patient data. 
- You MUST NOT commit patient-specific clinical data. 
- You MUST NOT share un-anonymized `.mat`, `.csv`, `.xlsx`, or `.json` files.
- You MUST ensure all patient parameters are fully pseudo-anonymized before pushing any configurations up to GitHub.

## Reporting a Vulnerability

If you discover a vulnerability or potential patient-data leak, do NOT open a public issue. 
Please email the repository owners directly, or notify through the internal communication channel.
