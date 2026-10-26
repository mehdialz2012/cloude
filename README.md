# Conditionally-activated lytic peptide actuator — prostate cancer (PSMA × pHe)

Computational design campaign for a lytic peptide actuator that is catalytically silent in
blood and normal tissue and becomes lytic only inside a tumour, via **conditional activation**
(AND-logic: extracellular acidity + PSMA surface-enzyme cleavage) rather than attenuated
intrinsic potency.

PI: Mehdi Alizadeh. Compute environment: Kaggle notebooks (free tier).

## Status

| # | Deliverable | Status |
|---|---|---|
| 1 | §4 adversarial audit verdict | **Done** — `docs/01_audit_verdict.md` |
| — | Target / warhead / architecture verdict | Folded into audit (§4, §6) |
| 2 | FTO / patent scan (US7456146 + family) | Next |
| 3 | Compute budget + Kaggle method redesign | Pending |
| 4 | Ranked candidates + 4-state truth table + notebooks | Pending (gated on WP1) |

## Headline from the audit

- **Proceed**, with three changes. The 2001/2009 PSMA-amoebapore precedent (positive, unrefuted
  in vivo data) was abandoned for **institutional** reasons (PI left in 2009), not scientific —
  the project's direct precedent worked.
- **Architectural correction:** the mask must block **assembly** while permitting pH-gated
  binding (C9 as written defeats C6's spatial-focusing mechanism).
- **The pH gate is promoted to Work Package 1** and made a hard architecture gate: its clean
  selectivity window decides AND-logic vs a single excellent gate.
- Novelty narrowed (Pramanik 2026, Montaner 2025); PSMA kept for patient-selectability.

See `docs/01_audit_verdict.md` for the full argument, claim-by-claim verdicts, the running
record of disagreements with the brief, and primary-source citations.

## Layout

```
docs/        written deliverables (audit, budget, FTO, risk register)
notebooks/   runnable Kaggle notebooks (checkpoint-restart, version-pinned installs)
data/        parameter/provenance files, intermediate-state manifests
```
