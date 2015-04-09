# tests getting statistics from different types of objects

context("getting statistics from objects")

# statistical results

c1 <- new("statistical_results",
          statistic_data = list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
                            counts = c(a1 = 5, a2 = 10, a3 = 1),
                            odds = c(a1 = 20, a2 = 100, a3 = 0)),
          annotation_id = c("a1", "a2", "a3"))

c2 <- data.frame(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
                 counts = c(a1 = 5, a2 = 10, a3 = 1),
                 odds = c(a1 = 20, a2 = 100, a3 = 0),
                 stringsAsFactors = FALSE)
rownames(c2) <- c("a1", "a2", "a3")

test_that("basic statistical_results works", {
  expect_equal(extract_statistics(c1), c2)
})

# combined_enrichment

stat1 <- new("statistical_results",
             annotation_id = c("a1", "a2", "a3"),
             statistic_data = list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
                               counts = c(a1 = 5, a2 = 10, a3 = 1),
                               odds = c(a1 = 20, a2 = 100, a3 = 0)))
stat2 <- new("statistical_results",
             annotation_id = c("a1", "a2", "a4"),
             statistic_data = list(pvalues = c(a1 = 0.01, a2 = 0.03, a4 = 0.0001),
                               counts = c(a1 = 5, a2 = 10, a4 = 1),
                               odds = c(a1 = 20, a2 = 100, a4 = 0)))

en1 <- new("enriched_result",
           features = letters,
           universe = letters,
           annotation = new("annotation"),
           statistics = stat1)

en2 <- new("enriched_result",
           features = letters,
           universe = letters,
           annotation = new("annotation"),
           statistics = stat2)

test_combined <- new("combined_enrichment",
                     enriched = list(en1 = en1, en2 = en2),
                     annotation = new("annotation"),
                     graph = new("graphNEL"))

out_stats_data <- matrix(NA, nrow = 4, ncol = 6)
rownames(out_stats_data) <- c("a1", "a2", "a3", "a4")
colnames(out_stats_data) <- c("en1.pvalues", "en1.counts", "en1.odds",
                                  "en2.pvalues", "en2.counts", "en2.odds")

en1_locs <- stat1@annotation_id
en2_locs <- stat2@annotation_id
out_stats_data[en1_locs, "en1.pvalues"] <- stat1@statistic_data$pvalues
out_stats_data[en1_locs, "en1.counts"] <- stat1@statistic_data$counts
out_stats_data[en1_locs, "en1.odds"] <- stat1@statistic_data$odds
out_stats_data[en2_locs, "en2.pvalues"] <- stat2@statistic_data$pvalues
out_stats_data[en2_locs, "en2.counts"] <- stat2@statistic_data$counts
out_stats_data[en2_locs, "en2.odds"] <- stat2@statistic_data$odds
out_stats_data <- as.data.frame(out_stats_data)

sig_matrix <- matrix(FALSE, 4, 2)
rownames(sig_matrix) <- rownames(out_stats_data)
colnames(sig_matrix) <- c("en1", "en2")

out_stats_combined <- new("combined_statistics",
                          statistic_data = out_stats_data,
                          annotation_id = rownames(out_stats_data),
                          which_enrichment = c("en1", "en1", "en1", "en2", "en2", "en2"),
                          which_statistic = c("pvalues", "counts", "odds", "pvalues", "counts", "odds"),
                          significant = new("significant_annotations",
                                           sig_calls = character(),
                                           significant = sig_matrix,
                                           measured = sig_matrix))

test_that("combined_enrichment works", {
  expect_equal(extract_statistics(test_combined), out_stats_combined)
})
