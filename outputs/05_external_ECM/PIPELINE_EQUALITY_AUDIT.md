# PIPELINE_EQUALITY_AUDIT

- Real Reactome, real NABA, and matched random sets use cohort-internal mean-z scoring.
- Real and random models use fixed k=3 and full-data k-means nstart=500.
- Real and random bootstrap use the same 100 bootstrap index vectors saved in bootstrap_indices.rds.
- Bootstrap refit uses bootstrap samples with replacement, retaining duplicate samples.
- Bootstrap refit initializes from full-data centroids and uses the same Hungarian centroid matching.
- Collapse rule is identical: any missing cluster or minimum predicted cluster size <5 yields collapse=TRUE and ARI/NMI=NA.
- Reactome random cross-definition ARI is compared against fixed NABA cluster.
- NABA random cross-definition ARI is compared against fixed Reactome cluster.
- Projection formula, nearest-centroid margin, and low-confidence cutoff are identical for real and random sets.
- original32 is excluded from external specificity empirical P values and BH_FDR_external8.
