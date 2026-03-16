# Changelog

All notable changes to the Unified VSD Lumped-Parameter Model project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Sequential L-BFGS-B optimization pipeline with Sobol masking.
- Strict 5% gate targets for clinical indices (QpQs, SAP_mean, LVEF).
- Export functions for Phase 0 vs Final GSA Sobol overlay.

## [1.0.0] - 2026-03-16
### Added
- Initial validated physiology codebase for cardiovascular VSD defect modeling.
- Core physics equations mapping right and left circuits, compliance, resistances, and time-varying elastances.
