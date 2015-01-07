# runs the tests for getting significant results from combined objects

library("categoryComparev2")

context("basic logical indexing of list elements")

c1 <- list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
           counts = c(a1 = 5, a2 = 10, a3 = 1),
           odds = c(a1 = 20, a2 = 100, a3 = 0))

pvalue1 <- c(TRUE, FALSE, TRUE)

test_that("list indexing works properly",
          expect_equal(pvalue1, multi_query_list(c1, pvalues < 0.05))
          )