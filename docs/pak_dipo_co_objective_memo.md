# Pak Dipo Revision Memo: CO Comparator and Objective Restructure

Date: 2026-05-08

## Why We Changed The Objective

Pak Dipo's concern is valid: a model can fit clinical outputs while hiding
non-physiological fitted parameters. In Reyna, there is an additional issue:
the systemic flow target, systemic pressure, right atrial pressure, and derived
SVR are not independent facts. They describe the same systemic load relation:

```text
SVR = (MAP - RAP) / Qs
```

The previous objective penalized `CO_Lmin` as an independent primary metric,
then penalized it again in a systemic consistency bundle, then again in a
clinical guard. That made the optimization overreact to one internally
inconsistent clinical block.

## CO Comparator Decision

For Reyna pre-surgery, catheter Fick CO is treated as systemic flow `Qs`.
The model comparator remains:

```text
metrics.CO_Lmin = metrics.Qs_Lmin
```

This is defensible because Fick CO estimates systemic oxygen-delivery flow.
In unrepaired VSD, LV stroke output includes recirculated shunt volume, so
`LVCO_Lmin` should align with `Qp_Lmin`, not with Fick `Qs`.

## Known Clinical Consistency Limitation

Reyna's Teichholz LV volumes imply `LVSV*HR` below the Fick/QpQs-implied
pulmonary flow. The model should not pretend both are exact point targets.
Therefore, `CO_Lmin` now has an explicit widened uncertainty:

```text
CO_comparator = Qs_Lmin
CO_uncertainty_Lmin = 0.50 L/min
```

## New Objective Behavior

The systemic quartet is now handled as one uncertainty-aware bundle:

```text
SAP_mean, RAP_mean, CO_Lmin/Qs, derived SVR
```

Those metrics are removed from the independent objective loop when the bundle
is active. This preserves the clinical relationship without double-counting.

## VSD Shunt Handling

Reyna pre-surgery now declares `VSD_mode = orifice_bidirectional`, so the shunt
is anchored to measured defect geometry, pressure gradient, and shunt flow.
A small shunt-flow consistency term checks agreement between:

```text
model Qvsd
model Qp - Qs
clinical Qs * (Qp/Qs - 1)
```

This does not add a new patient-specific hack; it makes the existing clinical
shunt accounting explicit.

## Preschool Scaling Framing

Reyna is in a preschool age range where the adult-baseline-to-pediatric-scaling
anchor is extrapolated. The run manifest now records this and applies a small,
targeted prior relaxation for selected vascular terms. This is reported for
reproducibility rather than hidden inside the optimizer.
