
## Reviewer 1 Feedback, Comment 5: Getting Exact Percentages of Immune and Epithelium in Mixed_3D and PD_3D

load("eng.subset.integrated_NodeAligned.Robj")

# Getting condition-level % 
md <- eng.subset.integrated_HK@meta.data
res <- md %>%
  dplyr::filter(Condition %in% c("Mixed_3D","PD_3D")) %>%
  dplyr::count(Condition, CellClass.NodeAligned, name = "n") %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(total = sum(n), pct = 100 * n / total) %>%
  dplyr::filter(CellClass.NodeAligned %in% c("Immune","Mesenchyme")) %>%
  dplyr::select(Condition, CellClass = CellClass.NodeAligned, n, pct) %>%
  tidyr::pivot_wider(
    names_from = CellClass,
    values_from = c(n, pct),
    values_fill = 0) %>%
  dplyr::ungroup()
res

# Per-replicate % and summary across replicates
by_rep <- md %>%
  dplyr::filter(Condition %in% c("Mixed_3D","PD_3D")) %>%
  dplyr::count(Condition, Orig_ID, CellClass.NodeAligned, name = "n") %>%
  dplyr::group_by(Condition, Orig_ID) %>%
  dplyr::mutate(total = sum(n), pct = 100 * n / total) %>%
  dplyr::filter(CellClass.NodeAligned %in% c("Immune","Mesenchyme")) %>%
  dplyr::select(Condition, Orig_ID, CellClass = CellClass.NodeAligned, pct)
summary_by_cond <- by_rep %>%
  dplyr::group_by(Condition, CellClass) %>%
  dplyr::summarise(mean_pct = mean(pct), sd_pct = sd(pct), .groups = "drop")
by_rep
summary_by_cond
