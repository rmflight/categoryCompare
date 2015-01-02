# testing the basic enriched class

library(categoryComparev2)
library(GO.db)
load("test_data.RData")

bp_annotation_list <- estrogen_10_enrichment$annotation
bp_annotation <- new("annotation",
                     annotation_features = bp_annotation_list,
                     description = Term(names(bp_annotation_list)),
                     counts = sapply(bp_annotation_list, length),
                     type = "GO.BP")

e10 <- new("enriched_result",
           features = estrogen_10_enrichment$features,
           universe = estrogen_10_enrichment$universe,
           annotation = bp_annotation,
           statistics = estrogen_10_enrichment$enrich)

e48 <- new("enriched_result",
           features = estrogen_48_enrichment$features,
           universe = estrogen_48_enrichment$universe,
           annotation = bp_annotation,
           statistics = estrogen_48_enrichment$enrich)

bp_combined <- combine_enrichments(e10, e48)