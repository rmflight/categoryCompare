# testing the basic enriched class

load("test_data.RData")

e10 <- new("enriched_result",
           features = estrogen_10_enrichment$features,
           universe = estrogen_10_enrichment$universe,
           annotation = estrogen_10_enrichment$annotation,
           statistics = estrogen_10_enrichment$enrich)

e48 <- new("enriched_result",
           features = estrogen_48_enrichment$features,
           universe = estrogen_48_enrichment$universe,
           annotation = estrogen_48_enrichment$annotation,
           statistics = estrogen_48_enrichment$enrich)