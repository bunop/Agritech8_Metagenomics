
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

plot_alpha_diversity <- function(data, x, y, fill, facet, title, xlab, ylab, custom_order = NULL) {
  # Check if custom_order is provided
  if (!is.null(custom_order)) {
    data[[x]] <- forcats::fct_relevel(data[[x]], custom_order)
  }

  # Create the faceting formula
  facet_formula <- as.formula(paste("~", facet))

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill))) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    ggplot2::facet_wrap(facet_formula, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(plot)
}

# like the function above, but plot histograms with sd bars
plot_alpha_diversity_with_stddev <- function(data, x, y, fill, facet, title, xlab, ylab) {
  # Ensure the data has columns for mean, SD, and Metric
  data <- data %>%
    mutate(
      ymin = mean - sd,
      ymax = mean + sd
    )

  # Create the faceting formula
  facet_formula <- as.formula(paste("~", facet))

  # Generate the plot
  plot <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill))) +
    geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, position = position_dodge(0.9)) +
    facet_wrap(facet_formula, scales = "free_y") +
    theme_minimal() +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(plot)
}

plot_distances <- function(distance_df, column) {
  # define column names
  column.1 <- sym(paste0(column, ".1"))
  column.2 <- sym(paste0(column, ".2"))
  facet_formula <- as.formula(paste("~", column.1))

  ggplot2::ggplot(distance_df, ggplot2::aes(x = !!column.2, y = Distance, fill = !!column.2)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(facet_formula, scales = "free_y") +
    ggplot2::labs(title = "Between-Group Distances",
         x = column,
         y = "Distance") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = column,
      fill = column
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )}
