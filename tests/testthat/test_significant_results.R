# runs the tests for getting significant results from combined objects

context("basic logical indexing of list elements")

c1 <- list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
           counts = c(a1 = 5, a2 = 10, a3 = 1),
           odds = c(a1 = 20, a2 = 100, a3 = 0))

pvalue1 <- c(TRUE, FALSE, TRUE) # pvalues < 0.05
pvalue2 <- c(FALSE, FALSE, TRUE) # pvalues < 0.001

counts1 <- c(TRUE, TRUE, FALSE) # counts > 1
counts2 <- c(TRUE, FALSE, TRUE) # counts < 10

odds1 <- c(TRUE, TRUE, FALSE) # odds > 0

pvalue_counts1 <- c(TRUE, FALSE, TRUE) # pvalues < 0.5, counts >= 1

test_that("list indexing works properly", {
          expect_equal(pvalue1, multi_query_list(c1, pvalues < 0.05))
          expect_equal(pvalue2, multi_query_list(c1, pvalues < 0.001))
          expect_equal(counts1, multi_query_list(c1, counts > 1))
          expect_equal(counts2, multi_query_list(c1, counts < 10))
          expect_equal(odds1, multi_query_list(c1, odds > 0))
          expect_equal(pvalue_counts1, multi_query_list(c1, pvalues < 0.5, counts >= 1))
          expect_error(multi_query_list(c1, asdfghjkl))
          })

context("significant annotations from statistical_results")

test_stat <- new("statistical_results",
                  annotation_id = c("a1", "a2", "a3"),
                  statistic_data = list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
                    counts = c(a1 = 5, a2 = 10, a3 = 1),
                    odds = c(a1 = 20, a2 = 100, a3 = 0)))

test_that("get correct significant annotations", {
           expect_equal(c("a1", "a3"), get_significant_annotations(test_stat, pvalues < 0.05))
           expect_equal(c("a1", "a2"), get_significant_annotations(test_stat, odds > 10))
           expect_equal(c("a1", "a3"), get_significant_annotations(test_stat, pvalues < 0.05, counts >= 1))
})

# here we test what we get back from a combined_enrichment object
context("significant annotations from combined_enrichment")

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

combined_stats <- extract_statistics(test_combined)
test_combined@statistics <- combined_stats

# this uses pvalues < 0.05
meas_matrix <- matrix(c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE), nrow = 4, ncol = 2)
rownames(meas_matrix) <- c("a1", "a2", "a3", "a4")
colnames(meas_matrix) <- c("en1", "en2")
sig_matrix <- meas_matrix
sig_matrix["a2", ] <- c(FALSE, TRUE)

expected_significant_annotations <- new("significant_annotations",
                                        significant = sig_matrix,
                                        measured = meas_matrix,
                                        sig_calls = "pvalues < 0.05")

expected_sig_combined <- test_combined
expected_sig_combined@statistics@significant <- expected_significant_annotations

test_that("get correct sig annotations from combined_enrichment", {
  expect_equal(expected_sig_combined, get_significant_annotations(test_combined, pvalues < 0.05))
})
