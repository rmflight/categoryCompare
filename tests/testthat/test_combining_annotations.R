library("categoryComparev2")

context("combining annotations")

c1 <- list(a = c("1", "2"),
           b = c("5", "6"))

c2 <- list(a = c("2", "3"),
           b = c("6", "7"),
           c = c("9", "10"))

c3 <- list(b = c("8", "11"))

c1_c2 <- list(a = c("1", "2", "3"),
              b = c("5", "6", "7"),
              c = c("9", "10"))

c1_c3 <- list(a = c("1", "2"),
              b = c("5", "6", "8", "11"))

c2_c3 <- list(a = c("2", "3"),
              b = c("6", "7", "8", "11"),
              c = c("9", "10"))

c1_c2_c3 <- list(a = c("1", "2", "3"),
                 b = c("5", "6", "7", "8", "11"),
                 c = c("9", "10"))

test_that("combining annotation features works correctly", {
  expect_equal(c1_c2, combine_annotation_features(list(c1 = c1, c2 = c2)))
  expect_equal(c1_c3, combine_annotation_features(list(c1 = c1, c3 = c3)))
  expect_equal(c2_c3, combine_annotation_features(list(c2 = c2, c3 = c3)))
  expect_equal(c1_c2_c3, combine_annotation_features(list(c1 = c1, c2 = c2, c3 = c3)))
})

context("combining text lists")

c1 <- c(a = "a", b = "b", c = "c")
c2 <- c(a = "a", b = "b", c = "c", d = "d")
c3 <- c(b = "b", c = "c", d = "d")
c4 <- c(a = "b", b = "b", c = "c")

out_result <- c(a = "a", b = "b", c = "c", d = "d")
use_names <- c("a", "b", "c", "d")

test_that("combining text lists works correctly", {
  expect_equal(out_result, combine_text(list(c1 = c1, c2 = c2), use_names, "test"))
  expect_equal(out_result, combine_text(list(c1 = c1, c3 = c3), use_names, "test"))
  expect_equal(out_result, combine_text(list(c1 = c2, c3 = c3), use_names, "test"))
  expect_equal(out_result, combine_text(list(c1 = c1, c2 = c2, c3 = c3), use_names, "test"))
  expect_error(combine_text(list(c1 = c1, c4 = c4), use_names, "test"))
})