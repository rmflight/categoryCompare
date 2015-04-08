# testing the basic enriched class

library(categoryComparev2)
library(GO.db)
library(parallel)
options(mc.cores = 4)
load("test_data.RData")

bp_annotation_list <- estrogen_10_enrichment$annotation
bp_annotation <- new("annotation",
                     annotation_features = bp_annotation_list,
                     description = Term(names(bp_annotation_list)),
                     counts = sapply(bp_annotation_list, length),
                     type = "GO.BP")

e10 <- enriched_result(estrogen_10_enrichment$features, 
                       estrogen_10_enrichment$universe,
                       bp_annotation, 
                       new("statistical_results",
                            statistic_data = estrogen_10_enrichment$enrich[c("pvalues", "odds_ratio", "counts")],
                            annotation_id = estrogen_10_enrichment$enrich$stat_names))

e48 <- enriched_result(estrogen_48_enrichment$features,
                       estrogen_48_enrichment$universe,
                       bp_annotation,
                       new("statistical_results",
                            statistic_data = estrogen_48_enrichment$enrich[c("pvalues", "odds_ratio", "counts")],
                            annotation_id = estrogen_48_enrichment$enrich$stat_names))

enriched <- list(e10 = e10, e48 = e48)

bp_combined <- combine_enrichments(e10 = e10, e48 = e48)

# if there are no significant things set, then generate the full graph
bp_graph <- generate_annotation_graph(bp_combined) # default graph gen

bp_sig <- get_significant_annotations(bp_combined, pvalues <= 0.05, counts >= 2) # add significant annotations

bp_graph_overlap <- generate_annotation_graph(bp_combined, "overlap") # adding a different type of similarity metric

bp_graph_sig <- generate_annotation_graph(bp_sig) # generate graph from sig annotations

bp_graph_sig2 <- filter_annotation_graph(bp_graph, bp_sig) # filter a previously generated annotation graph

bp_combined_2 <- combine_enrichments("combined", e10 = e10, e48 = e48)

bp_sig <- get_significant_annotations(bp_combined, pvalues <= 0.05, counts >= 2)
