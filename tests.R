# testing the basic enriched class

library(categoryComparev2)
library(GO.db)
library(parallel)
options(mc.cores = 4)
load("test_data2.RData")

bp_annotation_obj <- new("annotation",
                     annotation_features = bp_annotation,
                     description = Term(names(bp_annotation)),
                     counts = sapply(bp_annotation, length),
                     type = "GO.BP")

c1 <- enriched_result(c1_enrichment$features,
                       c1_enrichment$universe,
                       bp_annotation_obj, 
                       new("statistical_results",
                            statistic_data = c1_enrichment$enrich[c("pvalues", "odds_ratio", "counts")],
                            annotation_id = c1_enrichment$enrich$stat_names))

c2 <- enriched_result(c2_enrichment$features,
                       c2_enrichment$universe,
                       bp_annotation_obj,
                       new("statistical_results",
                            statistic_data = c2_enrichment$enrich[c("pvalues", "odds_ratio", "counts")],
                            annotation_id = c2_enrichment$enrich$stat_names))

enriched <- list(c1 = c1, c2 = c2)

bp_combined <- combine_enrichments(c1 = c1, c2 = c2)

# if there are no significant things set, then generate the full graph
bp_graph <- generate_annotation_graph(bp_combined) # default graph gen

bp_sig <- get_significant_annotations(bp_combined, pvalues <= 0.05, counts >= 2) # add significant annotations

bp_graph_sig <- generate_annotation_graph(bp_sig) # generate graph from sig annotations

bp_table <- generate_table(bp_combined)
bp_table_sig <- generate_table(bp_sig)

bp_graph_sig2 <- filter_annotation_graph(bp_graph, bp_sig) # filter a previously generated annotation graph

bp_combined_2 <- combine_enrichments("combined", e10 = e10, e48 = e48)

bp_sig <- get_significant_annotations(bp_combined, pvalues <= 0.05, counts >= 2)
