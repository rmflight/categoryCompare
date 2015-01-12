# tests getting statistics from different types of objects

library("categoryComparev2")

context("getting statistics from objects")

c1 <- new("statistical_results",
          statistics = list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
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