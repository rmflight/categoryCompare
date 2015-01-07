# runs the tests for getting significant results from combined objects

library("categoryComparev2")

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

test_that("list indexing works properly",
          expect_equal(pvalue1, multi_query_list(c1, pvalues < 0.05))
          expect_equal(pvalue2, multi_query_list(c1, pvalues < 0.001))
          expect_equal(counts1, multi_query_list(c1, counts > 1))
          expect_equal(counts2, multi_query_list(c1, counts < 10))
          expect_equal(odds1, multi_query_list(c1, odds > 0))
          expect_equal(pvalue_counts1, multi_query_list(c1, pvalues < 0.5, counts >= 1))
          expect_error(multi_query_list(c1, crap)))