STEP 1 CHECKPOINT

Files:
- step01_births_clean.rds: Cleaned birth-level dataset (5-year restriction). Includes HRFB flags, socio-demographic factors, IR covariates (v312, v602), and PSU contextual share_sec_plus.
- step01_svy_design.rds: srvyr survey design object (v021 PSU, strata coalesce(v022,v023), weight v005/1e6).

Notes:
- HRFB: age<18, age>34, interval<24 months, parity>=4.
- Unit of analysis: births with b19<60 (or v008-b3 fallback).
- Contextual: share_sec_plus = proportion women with >= secondary education by PSU (IR).
- Sensitivity: refit primary binary model using one birth per woman (most recent birth) to assess within-woman clustering.

Saved on: 2026-03-28 20:26:14 GMT
