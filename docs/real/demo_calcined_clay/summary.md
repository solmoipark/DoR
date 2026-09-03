# DoRGems 시나리오 A 요약 — CC-Gwangju-2026

- 역할 `calcined_clay` → InverseGems 슬롯 `metakaolin` (materials override: `materials.dorgems_c9828bc8202c.yaml`)
- DoR 예측 출처: `bayes`; 재령(일): [1.0, 3.0, 7.0, 28.0, 90.0, 180.0, 365.0]
- 실행 모드: real xGEMS; 분위 실행 q05: ok, q50: ok, q95: ok

| 재령 d | α q05 % | α q50 % | α q95 % |
|---|---|---|---|
| 1 | 4.9 | 13.5 | 23.2 |
| 3 | 7.6 | 20.8 | 35.7 |
| 7 | 10.1 | 27.7 | 47.2 |
| 28 | 14.1 | 38.9 | 66.0 |
| 90 | 16.1 | 44.5 | 75.8 |
| 180 | 16.5 | 45.7 | 78.3 |
| 365 | 16.6 | 46.1 | 78.9 |

- OOD: flags=['sparse_role:calcined_clay:4_papers', 'mahalanobis_beyond_p99:1134.62'], score_pct=100.0, sparse_role=True
- 근거 유사배합: 5배합 / 1편 (prediction.json → evidence)

| 재령 d | 변수 | q05 | q50 | q95 |
|---|---|---|---|---|
| 1 | phase_mass__Portlandite | 0 | 0 | 0 |
| 1 | phase_volume__Portlandite | 0 | 0 | 0 |
| 1 | scalar__pH | 13.47 | 12.61 | 10.52 |
| 3 | phase_mass__Portlandite | 0.0007338 | 0 | 0 |
| 3 | phase_volume__Portlandite | 0 | 0 | 0 |
| 3 | scalar__pH | 13.51 | 13.33 | 12.48 |
| 7 | phase_mass__Portlandite | 0.003288 | 0 | 0 |
| 7 | phase_volume__Portlandite | 0 | 0 | 0 |
| 7 | scalar__pH | 13.48 | 13.45 | 13.11 |
| 28 | phase_mass__Portlandite | 0.008616 | 0 | 0 |
| 28 | phase_volume__Portlandite | 0 | 0 | 0 |
| 28 | scalar__pH | 13.44 | 13.45 | 13.41 |
| 90 | phase_mass__Portlandite | 0.01035 | 0 | 0 |
| 90 | phase_volume__Portlandite | 0 | 0 | 0 |
| 90 | scalar__pH | 13.43 | 13.46 | 13.37 |
| 180 | phase_mass__Portlandite | 0.0107 | 0 | 0 |
| 180 | phase_volume__Portlandite | 0 | 0 | 0 |
| 180 | scalar__pH | 13.43 | 13.46 | 13.35 |
| 365 | phase_mass__Portlandite | 0.01098 | 0 | 0 |
| 365 | phase_volume__Portlandite | 0 | 0 | 0 |
| 365 | scalar__pH | 13.43 | 13.46 | 13.34 |

모든 수치의 원본: `prediction.json`, `envelope.csv`, `runs/*/forward/time_series.csv`. 경고는 `manifest.json`.
