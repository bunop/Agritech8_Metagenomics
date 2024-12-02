
custom_heatmap <- function(ampvis2_object, ...) {
  # Create the heatmap object
  heatmap_obj <- ampvis2::amp_heatmap(
    ampvis2_object,
    ...
  )

  # Extract current y-axis labels
  current_labels <- levels(heatmap_obj$data$Display)

  # Pattern to match 'Unclassified Bacteria' labels
  pattern_bacteria <- "k__Bacteria_[a-f0-9]{32}"

  # Pattern to match 'Unclassified' labels
  pattern_unclassified <- "k__Unclassified_[a-f0-9]{32}"

  # Apply case_when function
  new_labels <- dplyr::case_when(
    grepl(pattern_bacteria, current_labels) ~ "Unclassified Bacteria",
    grepl(pattern_unclassified, current_labels) ~ "Unclassified",
    TRUE ~ current_labels
  )

  # Set names
  new_labels <- stats::setNames(new_labels, current_labels)

  # replace the y-axis labels
  heatmap_plot <- heatmap_obj +
    ggplot2::scale_y_discrete(labels = new_labels)

  return(heatmap_plot)
}
