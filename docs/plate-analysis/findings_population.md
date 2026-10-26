# Population-level findings: counts, divisions, deaths, clustering, morphology

Data: `well_summary.csv`, `frame_series.csv`, `plate_map.csv` and per-well `states_hoct.parquet` / `states_hoct.json` (48 wells, 576 frames, 2.5 min/frame, FOV 2.36 mm², 0.5 µm/px, 1 px² = 0.25 µm²).
Replication unit = well (2 replicate wells per condition × density). Effects are reported as ratios of condition mean to the mean of the matched No-coating wells at the same density (C1/D1, E1/F1, G1/H1); drugs vs E1/F1. No per-cell statistics.
Rates are per 1000 cells per hour; a Jurkat doubling time of 22 h corresponds to a per-cell division rate of 31.5 /1000/h.

Figures: `fig_population_counts.png`, `fig_population_growth.png`, `fig_population_clustering.png`, `fig_population_morphology.png`, `fig_population_drugs.png`.

## Key findings

1. **The observed count increase (median 4.9-fold in 24 h, range 1.9–8.0; apparent doubling time median 10.8 h) is roughly twice what real Jurkat proliferation can produce, and the excess is not a settling transient.** Counts rise exponentially in both halves of the movie (mean fold 2.1 in 0–12 h and 2.4 in 12–24 h, all coating wells), so there is no plateau after an early settling phase. The count-change budget for 12–24 h (mean of 32 coating wells) attributes the measured slope of 66 /1000/h to confirmed divisions 56 (of which only 15 are divisions of isolated cells) + non-division interior appearances net of losses/merges 21 + border exchange −9. Confidence **high** that fold change overstates proliferation; **medium** for the exact split.

2. **Division confirmation over-calls inside clusters, and the bias grows with density.** Per clustered cell-hour, confirmed divisions occur at 290 /1000/h vs 17 /1000/h for isolated cells (17×; a clustered cell would "divide" every 2.4 h, which is impossible). In-cluster events are 42 % / 63 % / 83 % of all confirmed divisions at 5000 / 10000 / 20000 cells/well, the total division rate correlates with cluster fraction across wells (r = 0.71), and it rises within wells over time (41 → 57 /1000/h) exactly as cluster fraction rises (5 → 22 %). The total `division_rate_per_1000cells_h` is therefore not a usable proliferation measure, particularly for density or coating comparisons. Confidence **high**.

3. **Corrected proliferation estimate: doubling time ≈ 18–25 h (24-h fold ≈ 2–2.5), i.e. about half of the observed fold change is detection/segmentation gain.** Bracketing: (lower bound) isolated-cell confirmed division rate 17 /1000/h → Td 41 h (strict criteria and track fragmentation reject many real divisions; median track length is only 1.1 h); (upper bound) count slope in the quietest window (6–12 h, before clusters build up and after the early transient) 30 and 46 /1000/h in E1/F1 → Td 15–23 h. Confirmed divisions do look like divisions where they are accepted: mothers swell ~11 % (range 0.94–1.22) in the 2 h before division, daughters are 0.50–0.52 of mother area at birth and recover to 0.75 within 2 h (consistent in all 48 wells). Confidence **medium**.

4. **The non-division count gain is not focus drift; it is progressive splitting of large objects and unconfirmed in-cluster events.** Median cell intensity changes by only +11 % (ratio of 21–24 h to 0–3 h: median 1.11, 0.97–1.2× in 38/48 wells; the 10 wells above 1.2× are mostly the low-intensity PDL and F6 wells) while median area falls 27 % (ratio 0.73, IQR 0.68–0.85) and the share of objects >300 µm² halves (31 % in 0–6 h → 14 % in 18–24 h). Interior appearances (median first area 129 µm² vs 174 µm² overall; 80 % touching another object at appearance) and lost ends (75 % touching) are dominated by segmentation changes within touching groups. Confidence **medium**.

5. **Coatings have only modest effects on apparent growth; anti-CD3 does not increase proliferation within 24 h.** 12–24 h count growth rates of all coatings are within 0.9–1.2× of matched No-coating wells (Collagen highest, 1.08–1.19×; anti-CD3 0.89–0.96×). Isolated-cell division rates: Retronectin 1.44 / 1.58 / 1.84× control at 5000 / 10000 / 20000 (both replicates above control at every density), anti-CD3 1.05–1.09×, Collagen 0.88–1.06×, Cell-Tak 1.03–1.35×. Confidence **medium** (two wells per condition; isolated rate is itself a lower-bound estimator).

6. **Rapamycin (10 nM) is the only drug with a clear antiproliferative signal.** 24-h count fold 1.89 (1.92 / 1.86) vs 4.26 (3.95 / 4.57) for E1/F1 (0.44×); 12–24 h growth rate 0.39×; confirmed divisions 0.74× (37 / 22 vs 40); isolated division rate 0.76× (16.7 / 8.7 vs 16.6); count declines after 20 h in A1. 2DG 100 µM, TWS119 1 µM, Metformin 150 µM and 2HC 5 mM show no reduction on any growth metric (ratios 1.0–1.4). Tumor coculture: fold 2.67 (0.63×) but 12–24 h growth 0.96× and only ~240 objects at t = 0 (vs ~590 in controls). Confidence **high** for the rapamycin direction, **medium** for its magnitude and for the null results.

7. **Death calls are too sparse and too inconsistent between replicates to compare conditions.** 702 deaths plate-wide; median well rate 0.43 /1000 cells/h (IQR 0.21–0.87) ≈ 1 % of cells per day, below the ~3–5 % dead expected in a Jurkat culture, so the heuristic under-calls. Replicates disagree by >10× in several pairs (Cell-Tak 10000: 2 vs 66; anti-CD3 20000: 25 vs 2; Retronectin 10000: 9 vs 31) and ~90 % of calls sit at the minimum passing death score (3/5). Elevated counts in Cell-Tak 5000 (27 / 33), PDL 10000 (29 / 34), Rapamycin (13 / 29) and tumor coculture (12 / 27) are qualitative at best. Confidence **low**.

8. **Clustering is mostly a density effect, with Cell-Tak, Collagen and anti-CD3 above and Retronectin below the density trend.** On No-coating wells cluster fraction scales as density^1.30 (0.028 → 0.077 → 0.236 at 5000 / 10000 / 20000), and the 6-h-window time courses of all wells move along the same curve, so the rise in clustering over time is explained by the rise in count. Density-corrected (observed/predicted): Cell-Tak 1.30–1.50×, Collagen 1.26–1.57×, anti-CD3 1.16–1.57×, Retronectin 0.58–0.82× (both replicates on the same side of 1 at every density). Confidence **medium-high**.

9. **Segmented areas are 1.5–2.6× a single Jurkat (113 µm² for 12 µm diameter) and depend on density and time, so only within-density, same-time comparisons are meaningful.** No-coating median area 245 / 178 / 168 µm² at 5000 / 10000 / 20000 and falling ~27 % over 24 h in every well. Within density, anti-CD3 (+18 / +46 / +8 %) and Retronectin (+22 / +20 / +16 %) increase median area; all coatings lower median eccentricity slightly (0.62–0.65 vs 0.65–0.68 for No coating; Cell-Tak and PDL roundest); solidity is 0.95–0.97 everywhere (no spread/irregular shapes). Confidence **medium** for the area differences; **low** for interpreting them as spreading or activation-induced growth (halo/clump inclusion is a plausible alternative).

10. **PDL wells at 5000 and 10000 (C5, D5, E5, F5) are anomalous and were excluded from quantitative conclusions.** They hold 5–10× fewer objects than matched wells (36 and 138 cells/mm² vs 213 and 543) while PDL 20000 (G5, H5) is normal (1191 cells/mm²); their early median intensity is also among the lowest on the plate (88–100; only F6, G5 and H5 are similarly low, 98–105; most wells 110–170). Confidence **high** that they are not comparable; cause (cells displaced to the well periphery, coating-related contrast/detection failure, or plating) cannot be determined from these data.

## Sanity issue: what drives the 4–8× count increase?

![counts](fig_population_counts.png)

![growth](fig_population_growth.png)

**Time course (frame_series.csv).** Counts rise from frame 0 onward with no plateau. Log-linear slopes (per h): 0–3 h 0.058–0.085, 3–12 h 0.054–0.056, 12–24 h 0.068–0.069 (means per density). Growth is at least as fast in the second half as in the first; a settling/arrival process would saturate within a few hours (10-µm cells sediment through a 48-well medium column in <1 h). So (b) "cells settling into focus" can at most explain part of the 0–3 h excess (slope 0.06–0.085 vs 0.054 afterwards; in 20000 wells the 0–3 h slope is 0.085).

| density | slope 0–3 h | slope 3–12 h | slope 12–24 h | fold 0–12 h | fold 12–24 h | fold 0–24 h |
|---|---|---|---|---|---|---|
| 5000 | 0.060 | 0.054 | 0.069 | 2.03 | 2.35 | 4.70 |
| 10000 | 0.058 | 0.056 | 0.068 | 2.11 | 2.38 | 4.99 |
| 20000 | 0.085 | 0.054 | 0.068 | 2.26 | 2.36 | 5.34 |

**Budget (per-track start/end types in states_hoct.parquet).** Net count change = (interior appearances − lost − merged) + (entered − left) + confirmed divisions − deaths; the identity closes exactly per well. Control wells E1/F1:

| window | mean n | count slope | confirmed div | isolated div | isolated div per isolated cell | in-cluster div per clustered cell | detection gain | border | cluster % | objects >300 µm² % | median area µm² | median intensity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0–6 h | 709 | 69.9 | 32.1 | 21.7 | 22.1 | 597 | 58.7 | −11.5 | 1.8 | 15.6 | 194 | 144 |
| 6–12 h | 981 | 38.2 | 31.3 | 17.7 | 18.2 | 536 | 28.6 | −17.0 | 2.6 | 15.0 | 201 | 142 |
| 12–18 h | 1416 | 68.2 | 39.2 | 15.5 | 16.7 | 345 | 45.0 | −13.7 | 6.9 | 8.1 | 180 | 144 |
| 18–24 h | 2017 | 66.0 | 47.2 | 11.8 | 13.5 | 276 | 32.1 | −17.1 | 12.9 | 3.9 | 163 | 148 |
| all | 1281 | 58.7 | 39.8 | 15.3 | 16.6 | 318 | 38.7 | −15.3 | 7.7 | 10.5 | 179 | 145 |

Mean of 32 coating wells (PDL 5k/10k excluded): count slope 83 / 46 / 73 / 69 per window; confirmed divisions 41 / 41 / 54 / 57; isolated divisions 23 / 18 / 18 / 13; detection gain 53 / 18 / 31 / 14; border −1 / −9 / −7 / −10; cluster % 4.6 / 6.4 / 12.7 / 21.9; objects >300 µm² 31 / 29 / 20 / 15 %.

Assessment of the four hypotheses:

- (a) **Real proliferation**: present. Accepted divisions have a genuine mitotic signature (mother area +11 % in the last 2 h, daughters 0.50–0.52 of mother at birth, recovery to ~0.75 by 2 h, identical in isolated and in-cluster events and in all 48 wells; figure panel b). The isolated-cell division rate (17 /1000/h) is stable or slowly declining with time, as expected for an unperturbed culture.
- (b) **Settling/arrival**: contributes only to the first ~3 h (detection gain 53–59 /1000/h in 0–6 h vs 14–32 later). Border exchange is net negative everywhere (more tracks end at the border than start there), so inflow from outside the FOV is not the driver.
- (c) **Focus drift**: not supported. Median intensity drifts upward by a median of +11 % (not downward, and not in step with the count), area declines steadily rather than stepping, and the counts show no step changes (`outlier_frames` median 1 per well, 75th percentile 6; three wells have >60 flagged frames but their count curves are smooth on the 6-h scale). The area decline tracks the disappearance of large (>300 µm²) objects, i.e. clumps/defocused blobs being resolved into several smaller objects as cells rearrange, not a global focal change.
- (d) **Over-called divisions**: strongly supported for in-cluster events. Rates per clustered cell-hour are 150–770 /1000/h (figure panel c), 10–40× the isolated-cell rate, and the share of in-cluster acceptances rises from 42 % (5000) to 83 % (20000). In addition, 11,120 of 12,035 linker candidates in E1 were rejected for `mother_touching`; the rejected candidate's second object becomes an "appeared" track, so unconfirmed real divisions and cluster rearrangements both feed the detection-gain term.

**Corrected, conservative proliferation estimate.** Lower bound from isolated-cell divisions: 17 /1000/h (Td ≈ 41 h; median over 48 wells, 12–23 /1000/h across densities, lower at 20000 where fewer cells are isolated). Upper bound from the quiet 6–12 h count slope in E1/F1: 30–46 /1000/h (Td 15–23 h). Best estimate: Td ≈ 18–25 h, 24-h fold ≈ 2–2.5, with the remainder of the observed 4–5× being segmentation/detection gain. Total detected area (a crude biomass proxy that is insensitive to splitting) grows 3.2–3.5× in E1/F1, intermediate between the count fold (4.0–4.6) and the proliferation estimate, which is consistent with part of the gain being newly detected material (defocused or clumped cells resolving) and part being over-segmentation.

**Reliability of division counts.** Isolated-cell confirmed divisions: specific (mitotic signature, rate physically plausible) but insensitive (strict criteria, fragmented tracks; use only as a relative measure between wells at similar density). In-cluster confirmed divisions: not reliable; they are mostly cluster rearrangements. The total `division_rate_per_1000cells_h` is inflated 1.5–5× depending on cluster fraction and should not be compared across densities.

## Q1. Cell count dynamics and growth

![counts](fig_population_counts.png)

24-h count fold (n24/n0, 5-frame smoothed), condition mean and replicate values:

| condition | 5000 | 10000 | 20000 |
|---|---|---|---|
| No coating | 5.64 (6.02 / 5.26) | 4.26 (3.95 / 4.57) | 4.38 (4.30 / 4.46) |
| Cell-Tak | 4.75 (4.66 / 4.84) | 5.22 (5.14 / 5.30) | 5.33 (5.29 / 5.38) |
| Collagen | 4.40 (4.66 / 4.13) | 6.46 (5.65 / 7.26) | 7.06 (6.14 / 7.98) |
| Retronectin | 4.19 (4.01 / 4.38) | 4.85 (4.82 / 4.88) | 5.21 (5.23 / 5.19) |
| PDL* | 4.67 (4.19 / 5.15) | 4.92 (4.64 / 5.20) | 5.02 (4.66 / 5.37) |
| anti-CD3 | 4.58 (4.13 / 5.02) | 4.25 (4.14 / 4.37) | 5.05 (4.39 / 5.72) |

Apparent count growth rate 12–24 h (d ln n/dt, per h; 0.0315 = Td 22 h) and ratio to matched No-coating:

| condition | 5000 | 10000 | 20000 | ratio 5000 | ratio 10000 | ratio 20000 |
|---|---|---|---|---|---|---|
| No coating | 0.068 (0.073 / 0.063) | 0.061 (0.057 / 0.065) | 0.063 (0.064 / 0.063) | – | – | – |
| Cell-Tak | 0.073 (0.068 / 0.077) | 0.069 (0.068 / 0.070) | 0.068 (0.069 / 0.067) | 1.07 | 1.13 | 1.07 |
| Collagen | 0.061 (0.061 / 0.061) | 0.066 (0.062 / 0.071) | 0.075 (0.075 / 0.075) | 0.90 | 1.08 | 1.19 |
| Retronectin | 0.061 (0.060 / 0.062) | 0.068 (0.067 / 0.068) | 0.070 (0.069 / 0.070) | 0.90 | 1.11 | 1.10 |
| PDL* | 0.090 (0.101 / 0.079) | 0.088 (0.088 / 0.088) | 0.070 (0.068 / 0.073) | 1.32 | 1.43 | 1.11 |
| anti-CD3 | 0.065 (0.064 / 0.066) | 0.054 (0.057 / 0.052) | 0.060 (0.054 / 0.067) | 0.95 | 0.89 | 0.96 |

\* PDL 5000/10000 wells have 5–10× fewer objects than matched wells (see Key finding 10); their high relative growth is from a very low base (25–175 objects) and is not interpreted.

Settling vs growth phase: all wells show a 0–3 h excess (slope 0.06–0.09 /h) followed by steady exponential growth at 0.054–0.075 /h; there is no identifiable plateau. Corrected growth (see above): Td ≈ 18–25 h for No-coating controls; coatings change the apparent rate by at most ±20 %, which is within the plausible range of segmentation effects (Collagen and Cell-Tak wells have the highest cluster fractions and the largest area declines, both of which inflate apparent growth). Interpretation: no coating produces a growth difference that can be separated from segmentation artefacts with 2 wells per condition.

## Q2. Division rate per condition

![growth](fig_population_growth.png)

Total confirmed division rate (per 1000 cells/h; mean (rep1 / rep2)) – **inflated by in-cluster calls, do not compare across densities**:

| condition | 5000 | 10000 | 20000 |
|---|---|---|---|
| No coating | 37.1 (39.2 / 34.9) | 39.8 (39.5 / 40.2) | 57.9 (53.9 / 61.9) |
| Cell-Tak | 35.1 (42.8 / 27.3) | 49.8 (55.7 / 43.8) | 64.0 (66.9 / 61.1) |
| Collagen | 47.2 (47.3 / 47.0) | 55.9 (59.4 / 52.4) | 59.1 (54.1 / 64.1) |
| Retronectin | 44.5 (40.6 / 48.3) | 53.8 (55.2 / 52.4) | 61.3 (57.4 / 65.2) |
| PDL* | 24.3 (25.1 / 23.4) | 37.0 (41.3 / 32.6) | 58.9 (59.2 / 58.5) |
| anti-CD3 | 40.9 (28.6 / 53.2) | 54.4 (57.8 / 51.0) | 60.8 (56.4 / 65.2) |

Isolated-cell confirmed division rate per isolated cell-hour (per 1000/h) – more specific, lower-bound estimator:

| condition | 5000 | 10000 | 20000 | ratio 5000 | ratio 10000 | ratio 20000 |
|---|---|---|---|---|---|---|
| No coating | 22.1 (23.0 / 21.3) | 16.6 (16.5 / 16.7) | 11.7 (12.1 / 11.3) | – | – | – |
| Cell-Tak | 22.8 (27.4 / 18.1) | 19.7 (22.4 / 17.0) | 15.9 (18.4 / 13.3) | 1.03 | 1.19 | 1.35 |
| Collagen | 23.5 (24.4 / 22.7) | 14.6 (14.7 / 14.5) | 10.4 (9.0 / 11.9) | 1.06 | 0.88 | 0.89 |
| Retronectin | 31.8 (28.3 / 35.4) | 26.1 (28.3 / 23.9) | 21.6 (17.9 / 25.3) | 1.44 | 1.58 | 1.84 |
| PDL* | 14.5 (16.0 / 13.1) | 15.2 (17.5 / 12.8) | 10.3 (9.4 / 11.2) | 0.66 | 0.91 | 0.88 |
| anti-CD3 | 23.4 (16.8 / 29.9) | 17.4 (19.8 / 15.1) | 12.8 (11.1 / 14.4) | 1.06 | 1.05 | 1.09 |

- The isolated rate falls with density (22 → 17 → 12 /1000/h for No coating) because at high density isolated cells are a shrinking, possibly non-random subset and daughters are more often rejected for touching; compare conditions within a density only.
- **anti-CD3**: isolated rate 1.05–1.09× control, total rate 1.05–1.37× (the 10000 excess is matched by a higher cluster fraction). No evidence of increased proliferation within 24 h, consistent with activation-induced blasting/division starting after 24–48 h (interpretation).
- **Retronectin** is the only coating with a consistent increase in isolated divisions (both replicates above control at every density, 1.44–1.84×). Interpretation uncertain: Retronectin wells also have the lowest cluster fractions (0.58–0.82× trend), so more divisions happen in isolation and pass the filter; this is at least partly a detectability effect.
- **Drugs** (vs E1/F1 = 39.8 total, 16.6 isolated): Rapamycin 29.3 (37.0 / 21.6) total, 12.7 (16.7 / 8.7) isolated → 0.74× / 0.76×, reduction in one replicate only but corroborated by the count data. 2DG 47.6 / 17.9, TWS119 50.2 / 20.2, Metformin 51.2 / 21.1, 2HC 46.1 / 19.2 → ratios 1.08–1.29, no reduction. Tumor coculture 29.8 / 16.2 → 0.75× / 0.98× (the lower total rate reflects the lower cluster-driven inflation at ~160 cells/mm²).

Reliability: isolated-division ratios between conditions at the same density are probably informative at the ±30 % level with 2 wells; total division rates are not.

## Q3. Death events

Death events (0–24 h, rep1 / rep2) and rates (per 1000 cells/h):

| condition | 5000 events | 10000 events | 20000 events | 5000 rate | 10000 rate | 20000 rate |
|---|---|---|---|---|---|---|
| No coating | 9 / 7 | 6 / 4 | 13 / 2 | 0.66 | 0.16 | 0.13 |
| Cell-Tak | 27 / 33 | 2 / 66 | 15 / 19 | 3.01 | 1.21 | 0.30 |
| Collagen | 9 / 23 | 17 / 3 | 5 / 18 | 0.89 | 0.26 | 0.23 |
| Retronectin | 16 / 9 | 9 / 31 | 5 / 23 | 0.67 | 0.54 | 0.25 |
| PDL* | 1 / 4 | 29 / 34 | 13 / 5 | 1.12 | 4.11 | 0.13 |
| anti-CD3 | 4 / 1 | 9 / 20 | 25 / 2 | 0.20 | 0.47 | 0.23 |
| Rapamycin | – | 13 / 29 | – | – | 1.15 | – |
| 2DG | – | 8 / 7 | – | – | 0.27 | – |
| TWS119 | – | 12 / 17 | – | – | 0.56 | – |
| Metformin | – | 8 / 15 | – | – | 0.44 | – |
| 2HC | – | 28 / 8 | – | – | 0.74 | – |
| Tumor coculture | – | 12 / 27 | – | – | 2.20 | – |

Credibility: (i) median rate 0.43 /1000/h ≈ 1 %/day, lower than the few-percent dead typical of Jurkat cultures → under-calling, as expected for a conservative heuristic; (ii) called deaths in high-count wells (F2, D2, B1) have plausible features – no touching neighbour, area 60–70 % of the well median, intensity 1.3–1.6× the well median (bright/condensed), but 90 % have the minimum passing score 3/5; (iii) replicate agreement is poor (E2 2 vs F2 66; G6 25 vs H6 2; E3 17 vs F3 3); (iv) deaths are ~0.1 % of all track ends (702 vs 454,236 "lost" ends), so dying cells that fragment or go out of focus are overwhelmingly booked as "lost". Conclusion: death rates can only be used qualitatively, and no condition effect is claimed. Consistently elevated pairs are Cell-Tak 5000 (27 / 33 vs 9 / 7), PDL 10000 (29 / 34, low-count wells), Rapamycin (13 / 29 vs 6 / 4) and tumor coculture (12 / 27 vs 6 / 4; may include tumor cells or Jurkat–tumor interactions); the rapamycin wells also show a count decline after 20 h in A1, which would be consistent with some loss, but 13 called deaths cannot account for a drop of ~200 objects.

## Q4. Clustering

![clustering](fig_population_clustering.png)

Cluster fraction 0–24 h (fraction of cell-frames in multi-cell clusters), mean (rep1 / rep2):

| condition | 5000 | 10000 | 20000 | cells/mm² 5000 | 10000 | 20000 |
|---|---|---|---|---|---|---|
| No coating | 0.028 (0.022 / 0.035) | 0.077 (0.081 / 0.073) | 0.236 (0.192 / 0.279) | 213 | 543 | 1108 |
| Cell-Tak | 0.036 (0.051 / 0.022) | 0.138 (0.167 / 0.109) | 0.266 (0.298 / 0.233) | 197 | 619 | 1082 |
| Collagen | 0.068 (0.069 / 0.067) | 0.153 (0.175 / 0.131) | 0.245 (0.277 / 0.212) | 316 | 657 | 1031 |
| Retronectin | 0.027 (0.023 / 0.031) | 0.085 (0.090 / 0.080) | 0.180 (0.203 / 0.157) | 331 | 664 | 1124 |
| PDL* | 0.039 (0.013 / 0.066) | 0.107 (0.137 / 0.078) | 0.343 (0.344 / 0.341) | 36 | 138 | 1191 |
| anti-CD3 | 0.049 (0.038 / 0.060) | 0.134 (0.153 / 0.116) | 0.225 (0.221 / 0.228) | 242 | 556 | 1007 |

Density control: a power law fitted to the six No-coating wells gives cluster fraction ∝ (cells/mm²)^1.30 (an exponent of 1 is expected for random encounters; the steeper slope indicates that clusters at high density are also larger/longer-lived). Observed / predicted cluster fraction:

| condition | 5000 | 10000 | 20000 | touching fraction ratio vs control (5000 / 10000 / 20000) |
|---|---|---|---|---|
| No coating | 1.08 (0.84 / 1.32) | 0.89 (0.88 / 0.90) | 1.07 (1.01 / 1.12) | – |
| Cell-Tak | 1.50 (1.56 / 1.43) | 1.35 (1.27 / 1.43) | 1.30 (1.09 / 1.51) | 1.07 / 1.38 / 1.02 |
| Collagen | 1.57 (1.63 / 1.51) | 1.37 (1.44 / 1.29) | 1.26 (1.10 / 1.42) | 2.08 / 1.78 / 1.05 |
| Retronectin | 0.58 (0.50 / 0.67) | 0.75 (0.78 / 0.73) | 0.82 (0.75 / 0.89) | 0.79 / 0.96 / 0.72 |
| PDL* | 13.3 (11.1 / 15.4) | 7.1 (8.2 / 6.1) | 1.43 (1.37 / 1.49) | 1.02 / 1.06 / 1.18 |
| anti-CD3 | 1.57 (1.50 / 1.65) | 1.49 (1.58 / 1.40) | 1.16 (1.13 / 1.20) | 1.41 / 1.57 / 0.96 |
| drugs (10000; control 0.89) | Rapamycin 1.05 (0.99 / 1.12), 2DG 1.15 (1.15 / 1.15), TWS119 1.15 (1.02 / 1.28), Metformin 1.25 (1.18 / 1.33), 2HC 1.03 (1.35 / 0.71), Tumor coculture 6.2 (6.9 / 5.6) | | | |

- Cell-Tak, Collagen and anti-CD3 sit 1.2–1.6× above the density trend at every density with both replicates on the same side; the excess shrinks at 20000 where crowding dominates. This matches the expectation that Cell-Tak and anti-CD3 promote adhesion/contact (interpretation), but Collagen shows it as strongly, so a general "coated surface" effect (or coating-related image contrast making touching cells more often segmented as clusters) is an alternative.
- Retronectin is consistently below the trend (0.58–0.82×), i.e. cells on Retronectin are more dispersed.
- PDL 5000/10000 values are extrapolations far below the fitted density range and are dominated by the well anomaly; PDL 20000 is 1.43× (both replicates), the highest cluster fraction on the plate (0.34).
- Tumor coculture at ~160 cells/mm² has a cluster fraction (0.11) that is 6× the Jurkat-only trend, consistent with Jurkat cells accumulating around tumor cells or tumor cells themselves forming clusters (interpretation; cell identities are not available).
- Time course (panel b): within every well the cluster fraction rises ~5-fold from 0–6 h to 18–24 h, but the per-window points follow the same density curve as the 24-h means, so clustering increases because density increases, not because of a time-dependent change in adhesion. anti-CD3 and Cell-Tak 10000 run parallel to and above the No-coating curve throughout; there is no sign of a delayed onset of clustering in anti-CD3 wells within 24 h.

## Q5. Morphology

![morphology](fig_population_morphology.png)

Median area (µm²), eccentricity and solidity, condition means (rep1 / rep2):

| condition | area 5000 | area 10000 | area 20000 | ecc 5000 | ecc 10000 | ecc 20000 | solidity 5000 / 10000 / 20000 |
|---|---|---|---|---|---|---|---|
| No coating | 245 (252 / 238) | 178 (174 / 183) | 168 (157 / 179) | 0.651 (0.643 / 0.658) | 0.675 (0.674 / 0.676) | 0.661 (0.663 / 0.660) | 0.966 / 0.958 / 0.956 |
| Cell-Tak | 210 (197 / 224) | 195 (164 / 225) | 170 (162 / 178) | 0.616 (0.624 / 0.607) | 0.630 (0.643 / 0.617) | 0.641 (0.641 / 0.642) | 0.965 / 0.961 / 0.957 |
| Collagen | 262 (270 / 255) | 198 (207 / 190) | 187 (162 / 211) | 0.641 (0.642 / 0.640) | 0.654 (0.645 / 0.662) | 0.655 (0.660 / 0.649) | 0.966 / 0.959 / 0.957 |
| Retronectin | 298 (299 / 297) | 214 (206 / 223) | 194 (178 / 210) | 0.629 (0.627 / 0.631) | 0.648 (0.658 / 0.637) | 0.638 (0.651 / 0.625) | 0.970 / 0.963 / 0.961 |
| PDL* | 202 (204 / 201) | 198 (224 / 172) | 160 (165 / 156) | 0.625 (0.627 / 0.623) | 0.619 (0.630 / 0.609) | 0.635 (0.635 / 0.635) | 0.965 / 0.962 / 0.954 |
| anti-CD3 | 289 (306 / 272) | 261 (239 / 283) | 182 (180 / 184) | 0.628 (0.616 / 0.639) | 0.634 (0.640 / 0.628) | 0.649 (0.648 / 0.651) | 0.969 / 0.964 / 0.957 |

Ratios vs matched No-coating – area: Cell-Tak 0.86 / 1.09 / 1.01, Collagen 1.07 / 1.11 / 1.11, Retronectin 1.22 / 1.20 / 1.16, PDL 0.83 / 1.11 / 0.96, anti-CD3 1.18 / 1.46 / 1.08; eccentricity: all coatings 0.92–0.99 (Cell-Tak 0.93–0.97, PDL 0.92–0.96); solidity 1.00–1.01 everywhere.

Area time trend (relative linear change over 24 h, `area_trend_rel`): No coating −0.13 / −0.24 / −0.34; Cell-Tak −0.04 / −0.34 / −0.33; Collagen −0.23 / −0.43 / −0.57; Retronectin −0.04 / −0.34 / −0.49; anti-CD3 −0.10 / −0.60 / −0.39 (10000 replicates −0.32 / −0.87). Every well's median area declines; the decline is steepest where cluster fraction grows most.

Interpretation and caveats:
- Absolute areas (160–300 µm²) are 1.5–2.6× the ~113 µm² expected for a 12-µm Jurkat and scale inversely with density (−39 % log-correlation across coating wells) and with time. This is a segmentation property (label-free halo and doublets/clumps segmented as single objects at low density; tighter masks and more splitting at high density), not cell size. Only within-density, same-time comparisons are meaningful, and the early share of objects >300 µm² (31 % in 0–6 h) means a large minority of early "cells" are multi-cell or defocused objects.
- Within density, anti-CD3 (+8 to +46 %) and Retronectin (+16 to +22 %) have larger objects in both replicates. An early activation-associated size increase under anti-CD3 is biologically plausible (interpretation), but anti-CD3 wells also have higher cluster fractions (1.2–1.6× trend) and the largest early areas (F6 0–3 h: 451 µm²), so clump inclusion is the more parsimonious explanation for part of the excess; Retronectin wells have *lower* clustering, so the Retronectin area excess is less confounded.
- Eccentricity 0.62–0.68 and solidity 0.95–0.97 in all conditions indicate mildly elongated, convex objects; coatings make objects slightly rounder (lower eccentricity) rather than more spread/irregular. There is no morphological evidence of spreading on adhesive coatings within 24 h (interpretation: Jurkat are weakly adherent; Cell-Tak/PDL rounding may reflect immobilised, less motile cells).

## Q6. Drug panel vs E1/F1

![drugs](fig_population_drugs.png)

Mean (rep A / rep B) and ratio to control mean:

| condition | count fold 24 h | growth 12–24 h (/h) | confirmed div (/1000/h) | isolated div (/1000/h) | death rate (/1000/h) | cluster fraction | density-corr. clustering | median area µm² | median ecc |
|---|---|---|---|---|---|---|---|---|---|
| Control E1/F1 | 4.26 (3.95 / 4.57) | 0.061 (0.057 / 0.065) | 39.8 (39.5 / 40.2) | 16.6 (16.5 / 16.7) | 0.16 (0.19 / 0.14) | 0.077 (0.081 / 0.073) | 0.89 (0.88 / 0.90) | 178 (174 / 183) | 0.675 (0.674 / 0.676) |
| Rapamycin 0.01 µM | 1.89 (1.92 / 1.86) **0.44×** | 0.024 (0.025 / 0.023) **0.39×** | 29.3 (37.0 / 21.6) 0.74× | 12.7 (16.7 / 8.7) 0.76× | 1.15 (0.55 / 1.75) 7.1× | 0.052 (0.061 / 0.044) 0.68× | 1.05 (0.99 / 1.12) 1.18× | 196 (173 / 218) 1.10× | 0.708 (0.720 / 0.696) 1.05× |
| 2DG 100 µM | 5.53 (5.78 / 5.28) 1.30× | 0.071 (0.074 / 0.068) 1.16× | 47.6 (48.2 / 47.1) 1.20× | 17.9 (19.4 / 16.4) 1.08× | 0.27 (0.33 / 0.21) 1.7× | 0.092 (0.073 / 0.112) 1.19× | 1.15 (1.15 / 1.15) 1.29× | 232 (262 / 203) 1.30× | 0.660 (0.652 / 0.668) 0.98× |
| TWS119 1 µM | 4.72 (4.99 / 4.44) 1.11× | 0.062 (0.068 / 0.056) 1.01× | 50.2 (47.6 / 52.7) 1.26× | 20.2 (19.6 / 20.9) 1.22× | 0.56 (0.49 / 0.64) 3.5× | 0.080 (0.067 / 0.093) 1.03× | 1.15 (1.02 / 1.28) 1.29× | 255 (244 / 266) 1.43× | 0.653 (0.657 / 0.649) 0.97× |
| Metformin 150 µM | 5.63 (5.63 / 5.62) 1.32× | 0.072 (0.072 / 0.073) 1.18× | 51.2 (48.9 / 53.5) 1.29× | 21.1 (21.8 / 20.4) 1.27× | 0.44 (0.29 / 0.58) 2.7× | 0.090 (0.087 / 0.093) 1.16× | 1.25 (1.18 / 1.33) 1.41× | 234 (211 / 258) 1.31× | 0.664 (0.677 / 0.650) 0.98× |
| 2HC 5 mM | 5.89 (7.25 / 4.54) 1.38× | 0.075 (0.083 / 0.067) 1.22× | 46.1 (47.4 / 44.9) 1.16× | 19.2 (16.4 / 22.0) 1.16× | 0.74 (1.19 / 0.28) 4.6× | 0.069 (0.083 / 0.055) 0.89× | 1.03 (1.35 / 0.71) 1.16× | 285 (341 / 229) 1.60× | 0.640 (0.630 / 0.651) 0.95× |
| Tumor coculture | 2.67 (2.64 / 2.69) 0.63× | 0.059 (0.059 / 0.059) 0.96× | 29.8 (29.3 / 30.4) 0.75× | 16.2 (16.0 / 16.4) 0.98× | 2.20 (1.29 / 3.11) 13.6× | 0.110 (0.127 / 0.094) 1.43× | 6.2 (6.9 / 5.6) 7.0× | 227 (240 / 214) 1.27× | 0.608 (0.608 / 0.609) 0.90× |

- **Rapamycin**: the only drug with a concordant reduction across all growth metrics in both replicates (fold 0.44×, 12–24 h growth 0.39×, isolated divisions 0.76×). Counts in A1 peak at 20 h (~1310) and fall to ~1080 by 24 h; the border term is the most negative on the plate (−28 and −19 /1000/h), so part of the low fold may be objects leaving the FOV or becoming undetectable rather than non-division. Objects are also the most elongated on the plate (ecc 0.708 vs 0.675) and clustering is lowest (0.68× control). Interpretation: consistent with mTOR inhibition slowing cell-cycle progression; the ~2-fold reduction in apparent growth is a plausible magnitude, but the 24-h fold of 1.9 is itself inflated by the same detection gain as the controls, so the real rapamycin growth may be close to zero.
- **2DG, TWS119, Metformin, 2HC**: no reduction in count growth, confirmed divisions or isolated divisions (all 1.0–1.4× control). These wells started at lower counts (290–550 vs 540–630 objects) and have larger median areas (232–285 vs 178 µm²; 2HC A5 341 µm² with 0–3 h median 387 µm²), which points to different initial clumping/focus rather than a biological effect. Interpretation: at these doses and within 24 h, glycolysis, GSK3, complex-I and ACLY inhibition do not slow Jurkat expansion measurably in this assay; the elevated death-rate ratios (1.7–4.6×) rest on 7–28 events per well with poor replicate agreement and are not interpretable.
- **Tumor coculture**: 24-h fold 2.67 (0.63×) but the 12–24 h growth rate (0.96×) and isolated division rate (0.98×) equal control; the initial object count (~240) is less than half the control's for the same total cell number, so either tumor cells are poorly detected or settle outside the FOV. Clustering is 6–7× the Jurkat density trend and eccentricity is the lowest (0.608). Without cell identities, nothing can be attributed to Jurkat vs tumor cells.

## Caveats

- Two wells per condition × density; all effect sizes are means of two values with the replicates shown. No statistical tests were performed; a ratio is called consistent only when both replicates lie on the same side of the control.
- `count_fold_change` and `apparent_doubling_time_h` include detection/segmentation gain of roughly the same magnitude as real proliferation (Key findings 1–4); use the 12–24 h slope, the isolated-cell division rate, or the corrected Td ≈ 18–25 h estimate, and treat condition differences below ~20 % as unresolved.
- `division_rate_per_1000cells_h` is dominated by in-cluster acceptances whose per-cell rate is physically impossible; it is biased upward by cluster fraction (and therefore by density, coating and time). Isolated-cell divisions are specific but under-count (strict criteria; median track length 1.1 h; `short_parent`/`short_daughter` rejections).
- Death calls are rare (median 0.43 /1000/h), mostly minimum-score, and replicate-inconsistent; only 0.15 % of track ends are deaths vs 454k "lost" ends.
- Tracking QC: `bidir_fraction` is 0.64–0.93 (median 0.85); the lowest values are the rapamycin wells (A1 0.64, B1 0.64) and the densest 20000 wells (G3 0.65, G6 0.68, H3 0.70, H1 0.71), so per-track metrics (including divisions) are least reliable there.
- Segmented area is a function of density, clustering and focus as well as cell size; eccentricity/solidity differences between conditions are ≤5 % and within what segmentation differences alone could produce.
- PDL 5000/10000 wells (C5, D5, E5, F5) hold 5–10× fewer objects than their matched wells and have among the lowest early intensities on the plate; they are flagged (*) and excluded from quantitative conclusions. PDL 20000 wells look normal.
- Tumor coculture wells contain two cell types that cannot be distinguished in these outputs.
- Border exchange is net negative in almost all wells (−1 to −28 /1000/h), which is implausible for a homogeneous well and suggests an asymmetry in how border starts/ends are classified; it was not used in the corrected growth estimate.
