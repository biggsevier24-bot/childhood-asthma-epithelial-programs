source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 04_assign_nearest_centroid.R")

scores <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_program_scores.tsv"))
centroids <- projection_read_tsv(file.path(projection_out, "GSE152004_reference_centroids.tsv"))
state_order <- c("E1", "E2", "E3")
centroids <- centroids[match(state_order, centroids$state), ]
if (anyNA(centroids$state)) stop("Reference centroid state missing")

x <- as.matrix(scores[, c("T2_projection", "IFN_projection", "repair_projection")])
c_mat <- as.matrix(centroids[, c("T2", "IFN", "repair")])
dist_mat <- sapply(seq_len(nrow(c_mat)), function(i) {
  sqrt(rowSums((x - matrix(c_mat[i, ], nrow(x), ncol(x), byrow = TRUE))^2))
})
colnames(dist_mat) <- state_order
ord <- t(apply(dist_mat, 1, function(z) order(z, seq_along(z))))
nearest <- dist_mat[cbind(seq_len(nrow(dist_mat)), ord[, 1])]
second <- dist_mat[cbind(seq_len(nrow(dist_mat)), ord[, 2])]
margin <- second - nearest
threshold <- as.numeric(projection_config$projection$low_confidence_threshold)

distances <- data.frame(
  sample_id = scores$sample_id,
  tissue = scores$tissue,
  T2 = scores$T2_projection,
  IFN = scores$IFN_projection,
  repair = scores$repair_projection,
  distance_to_E1 = dist_mat[, "E1"],
  distance_to_E2 = dist_mat[, "E2"],
  distance_to_E3 = dist_mat[, "E3"],
  assigned_state = state_order[ord[, 1]],
  nearest_distance = nearest,
  second_nearest_distance = second,
  assignment_margin = margin,
  assignment_confidence = ifelse(margin <= threshold, "low", "standard")
)
projection_write_tsv(distances, file.path(projection_out, "GSE118761_projection_distances.tsv"))

assignments <- data.frame(
  sample_id = scores$sample_id,
  tissue = scores$tissue,
  assigned_state = distances$assigned_state,
  distance_E1 = distances$distance_to_E1,
  distance_E2 = distances$distance_to_E2,
  distance_E3 = distances$distance_to_E3,
  nearest_distance = nearest,
  second_nearest_distance = second,
  margin = margin,
  assignment_confidence = distances$assignment_confidence,
  atopy = scores$atopy,
  asthma = scores$asthma,
  wheeze = scores$wheeze,
  age = if ("age" %in% names(scores)) scores$age else NA,
  sex = if ("sex" %in% names(scores)) scores$sex else NA_character_
)
if (anyDuplicated(assignments$sample_id)) stop("Duplicate final sample assignment")
projection_write_tsv(assignments, file.path(projection_out, "GSE118761_projection_assignments_FINAL.tsv"))
projection_log("DISTANCE Euclidean; feature order T2, IFN, repair; tie order E1,E2,E3")
projection_log("END 04_assign_nearest_centroid.R")
