library("sleuth")

dir <- "/home/manager/course_data/rna_seq/data/"
s2c <- read.table(paste0(dir, "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)

sample_id <- c("MT1", "MT2", "SBP1", "SBP2", "SBP3")
kal_dirs <- sapply(paste0(dir, sample_id), function(id) file.path(id))
s2c <- dplyr::select(s2c, sample = sample, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

t2g <- read.table(paste0(dir, "PccAS_v3.desc"), header=TRUE, sep="\t")

so <- sleuth_prep(s2c, ~ condition, extra_bootstrap_summary = TRUE, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_wt(so, 'conditionSBP', 'full')

results_table <- sleuth_results(so, 'conditionSBP', test_type = 'wt')

write.table(results_table, file=paste0(dir, "kallisto.results"), quote=FALSE, sep="\t", row.names=FALSE)

sleuth_live(so)