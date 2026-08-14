#!/usr/bin/env Rscript

# Render the compact, self-contained BRCA1 figure embedded by the repository
# README.  The detailed all-node-support figure remains a separate audit
# artifact.  This renderer is deliberately example-specific: it validates the
# fixed tree/metadata/clade contracts, never reroots or ladderizes the tree,
# and uses no network access.

fail <- function(message) {
  stop(message, call. = FALSE)
}

usage_text <- function() {
  paste(
    "Usage: render_brca1_readme.R --tree <rooted-tree.nwk>",
    "--metadata <sequence_metadata.tsv> --clade-support <clade_support.tsv>",
    "--out-prefix <new-prefix> [--width <inches>] [--height <inches>]",
    "[--preview-png <new-preview.png>]"
  )
}

parse_cli <- function(arguments) {
  if (length(arguments) == 1L && arguments[[1L]] %in% c("--help", "-h")) {
    cat(usage_text(), "\n")
    quit(status = 0L, save = "no")
  }
  values <- list(width = "11.5", height = "7", preview_png = "")
  aliases <- c(
    "--tree" = "tree",
    "--metadata" = "metadata",
    "--clade-support" = "clade_support",
    "--out-prefix" = "out_prefix",
    "--width" = "width",
    "--height" = "height",
    "--preview-png" = "preview_png"
  )
  index <- 1L
  while (index <= length(arguments)) {
    flag <- arguments[[index]]
    key <- unname(aliases[flag])
    if (length(key) == 0L || is.na(key)) {
      fail(paste("Unknown argument:", flag, "\n", usage_text()))
    }
    if (index == length(arguments)) {
      fail(paste("Missing value for", flag))
    }
    values[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  required <- c("tree", "metadata", "clade_support", "out_prefix")
  missing <- required[vapply(
    required,
    function(key) is.null(values[[key]]) || values[[key]] == "",
    logical(1)
  )]
  if (length(missing) > 0L) {
    fail(paste("Missing required arguments:", paste(missing, collapse = ", ")))
  }
  values
}

require_packages <- function(packages) {
  missing <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]
  if (length(missing) > 0L) {
    fail(paste(
      "Missing local R packages:", paste(missing, collapse = ", "),
      "Install versions compatible with the active R/Bioconductor release before rerunning."
    ))
  }
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = character()
  )
}

read_selected_metadata <- function(path) {
  metadata <- read_tsv(path)
  required <- c(
    "tip_id", "analysis_role", "inclusion_status", "species", "accession"
  )
  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0L) {
    fail(paste("Metadata is missing columns:", paste(missing, collapse = ", ")))
  }
  selected <- metadata[metadata$inclusion_status == "selected", , drop = FALSE]
  if (
    nrow(selected) == 0L ||
      any(selected$tip_id == "") ||
      any(selected$species == "") ||
      any(selected$accession == "") ||
      anyDuplicated(selected$tip_id)
  ) {
    fail("Selected metadata require unique non-empty tip, species, and accession values.")
  }
  allowed_roles <- c("study", "expanded", "outgroup")
  invalid_roles <- setdiff(unique(selected$analysis_role), allowed_roles)
  if (length(invalid_roles) > 0L) {
    fail(paste("Unsupported analysis_role values:", paste(invalid_roles, collapse = ", ")))
  }
  if (
    sum(selected$analysis_role == "study") != 1L ||
      sum(selected$analysis_role == "outgroup") < 2L
  ) {
    fail("The README figure requires one study tip and at least two outgroup tips.")
  }
  selected$tip_fontface <- ifelse(
    selected$analysis_role == "study", "bold.italic", "italic"
  )
  selected$label_color <- ifelse(
    selected$analysis_role == "study",
    "#E69F00",
    ifelse(selected$analysis_role == "outgroup", "#666666", "#262626")
  )
  # ggtree metadata attachment uses the first column as the join key.
  selected[, c("tip_id", setdiff(colnames(selected), "tip_id")), drop = FALSE]
}

descendant_tips <- function(tree, node) {
  tip_count <- length(tree$tip.label)
  if (node <= tip_count) {
    return(tree$tip.label[[node]])
  }
  children <- tree$edge[tree$edge[, 1L] == node, 2L]
  if (length(children) == 0L) {
    fail(paste("Internal node has no descendants:", node))
  }
  unlist(lapply(children, function(child) descendant_tips(tree, child)), use.names = FALSE)
}

validate_root_split <- function(tree, outgroup_tips) {
  if (!ape::is.rooted(tree)) {
    fail("The README figure requires a structurally rooted input tree.")
  }
  roots <- setdiff(unique(tree$edge[, 1L]), unique(tree$edge[, 2L]))
  if (length(roots) != 1L) {
    fail("The input tree must contain exactly one structural root.")
  }
  root_children <- tree$edge[tree$edge[, 1L] == roots[[1L]], 2L]
  child_tip_sets <- lapply(root_children, function(node) sort(descendant_tips(tree, node)))
  if (!any(vapply(child_tip_sets, identical, logical(1), sort(outgroup_tips)))) {
    fail("The structural root split does not isolate the selected outgroups.")
  }
}

read_broad_clades <- function(path, tree) {
  clades <- read_tsv(path)
  required <- c(
    "clade_id", "target_tip_count", "target_tips", "recovery_status", "node_id"
  )
  missing <- setdiff(required, colnames(clades))
  if (length(missing) > 0L) {
    fail(paste("Clade-support table is missing columns:", paste(missing, collapse = ", ")))
  }
  broad_names <- c("Mammalia", "Sauropsida", "Amphibia")
  broad <- clades[match(broad_names, clades$clade_id), , drop = FALSE]
  if (any(is.na(broad$clade_id)) || any(broad$recovery_status != "recovered")) {
    fail("Mammalia, Sauropsida, and Amphibia must all be recovered before highlighting.")
  }
  tip_count <- length(tree$tip.label)
  maximum_node <- tip_count + tree$Nnode
  broad$node <- suppressWarnings(as.integer(broad$node_id))
  if (any(is.na(broad$node)) || any(broad$node <= tip_count | broad$node > maximum_node)) {
    fail("Every broad-clade node_id must identify an internal node in the supplied tree.")
  }
  broad$tip_list <- strsplit(broad$target_tips, ",", fixed = TRUE)
  for (index in seq_len(nrow(broad))) {
    expected <- sort(broad$tip_list[[index]])
    if (
      length(expected) != as.integer(broad$target_tip_count[[index]]) ||
        any(expected == "") ||
        anyDuplicated(expected)
    ) {
      fail(paste("Invalid target-tip declaration for", broad$clade_id[[index]]))
    }
    observed <- sort(descendant_tips(tree, broad$node[[index]]))
    if (!identical(observed, expected)) {
      fail(paste(
        "Clade node/tip declaration does not match the supplied tree:",
        broad$clade_id[[index]]
      ))
    }
  }
  broad
}

parse_support <- function(label) {
  if (is.na(label) || label == "" || label == "Root") {
    return(c(sh_alrt = NA_real_, ufboot = NA_real_))
  }
  pieces <- strsplit(label, "/", fixed = TRUE)[[1L]]
  if (length(pieces) != 2L) {
    fail(paste("Unexpected internal-node label:", label))
  }
  values <- suppressWarnings(as.numeric(pieces))
  if (any(is.na(values)) || any(values < 0 | values > 100)) {
    fail(paste("Invalid SH-aLRT/UFBoot support label:", label))
  }
  c(sh_alrt = values[[1L]], ufboot = values[[2L]])
}

support_table <- function(tree) {
  if (is.null(tree$node.label) || length(tree$node.label) != tree$Nnode) {
    fail("The README figure requires one node-label slot per internal node.")
  }
  parsed <- t(vapply(tree$node.label, parse_support, numeric(2)))
  support <- data.frame(
    node = length(tree$tip.label) + seq_len(tree$Nnode),
    original_label = tree$node.label,
    sh_alrt = parsed[, "sh_alrt"],
    ufboot = parsed[, "ufboot"],
    stringsAsFactors = FALSE
  )
  support$joint_strong <- with(
    support,
    !is.na(sh_alrt) & !is.na(ufboot) & sh_alrt >= 80 & ufboot >= 95
  )
  support$display_class <- ifelse(
    is.na(support$sh_alrt) | is.na(support$ufboot),
    "not-displayed",
    ifelse(support$joint_strong, "joint-strong", "below-joint-threshold")
  )
  support$weak_label <- ifelse(
    support$display_class == "below-joint-threshold",
    support$original_label,
    NA_character_
  )
  support
}

file_sha256 <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  unname(as.character(openssl::sha256(connection)))
}

write_settings <- function(path, settings) {
  table <- data.frame(
    setting = names(settings),
    value = unname(vapply(settings, as.character, character(1))),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    table,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

main <- function() {
  options <- parse_cli(commandArgs(trailingOnly = TRUE))
  width <- suppressWarnings(as.numeric(options$width))
  height <- suppressWarnings(as.numeric(options$height))
  if (!is.finite(width) || !is.finite(height) || width <= 0 || height <= 0) {
    fail("width and height must be positive finite numbers in inches.")
  }
  require_packages(c("ape", "ggplot2", "ggtree", "openssl", "svglite"))
  if (!isTRUE(capabilities("cairo"))) {
    fail("This renderer requires an R build with Cairo PDF support.")
  }
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggtree))

  tree <- ape::read.tree(options$tree)
  if (
    is.null(tree) || inherits(tree, "multiPhylo") ||
      length(tree$tip.label) == 0L || any(tree$tip.label == "") ||
      anyDuplicated(tree$tip.label)
  ) {
    fail("Newick input must contain exactly one tree with unique non-empty tips.")
  }
  if (
    is.null(tree$edge.length) || length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length) | tree$edge.length < 0) ||
      !any(tree$edge.length > 0)
  ) {
    fail("The README phylogram requires complete finite non-negative branch lengths.")
  }
  metadata <- read_selected_metadata(options$metadata)
  if (!identical(sort(tree$tip.label), sort(metadata$tip_id))) {
    fail("Tree tips must equal selected metadata tip IDs exactly.")
  }
  validate_root_split(
    tree,
    metadata$tip_id[metadata$analysis_role == "outgroup"]
  )
  broad_clades <- read_broad_clades(options$clade_support, tree)
  support <- support_table(tree)

  role_palette <- c(study = "#E69F00", expanded = "#009E73", outgroup = "#999999")
  role_shapes <- c(study = 23, expanded = 21, outgroup = 22)
  role_labels <- c(
    study = "Study sequence", expanded = "Expanded reference", outgroup = "Outgroup"
  )
  band_palette <- c(
    Mammalia = "#E8F1F8", Sauropsida = "#F7EAD7", Amphibia = "#EEEAF4"
  )
  band_text <- c(
    Mammalia = "Mammalia",
    Sauropsida = "Sauropsida",
    Amphibia = "Amphibia\n(outgroup)"
  )

  plot <- ggtree::ggtree(
    tree,
    layout = "rectangular",
    ladderize = FALSE,
    branch.length = "branch.length"
  ) %<+% metadata
  positions <- plot$data[, c("node", "x", "y", "isTip", "label"), drop = FALSE]
  if (any(!is.finite(positions$x)) || any(!is.finite(positions$y))) {
    fail("ggtree produced non-finite plotting coordinates.")
  }
  maximum_tree_x <- max(positions$x)
  species_offset <- 0.04
  accession_offset <- 0.69
  species_label_x <- maximum_tree_x + species_offset
  accession_label_x <- maximum_tree_x + accession_offset
  band_label_x <- maximum_tree_x + 1.36
  right_limit <- maximum_tree_x + 1.52

  tip_positions <- merge(
    positions[positions$isTip, c("label", "x", "y"), drop = FALSE],
    metadata,
    by.x = "label",
    by.y = "tip_id",
    all.x = TRUE,
    sort = FALSE
  )
  if (
    nrow(tip_positions) != nrow(metadata) ||
      any(is.na(tip_positions$species)) ||
      any(is.na(tip_positions$accession))
  ) {
    fail("Not every tree tip received aligned display metadata.")
  }

  band_labels <- data.frame(
    clade_id = broad_clades$clade_id,
    x = band_label_x,
    y = vapply(
      broad_clades$tip_list,
      function(tips) {
        tip_y <- positions$y[positions$isTip & positions$label %in% tips]
        if (length(tip_y) != length(tips)) {
          fail("Not every broad-clade tip received a plotting coordinate.")
        }
        mean(range(tip_y))
      },
      numeric(1)
    ),
    label = unname(band_text[broad_clades$clade_id]),
    stringsAsFactors = FALSE
  )
  for (index in seq_len(nrow(broad_clades))) {
    tips <- broad_clades$tip_list[[index]]
    clade_tip_max <- max(positions$x[positions$isTip & positions$label %in% tips])
    plot <- plot + ggtree::geom_hilight(
      node = broad_clades$node[[index]],
      type = "rect",
      to.bottom = TRUE,
      fill = unname(band_palette[broad_clades$clade_id[[index]]]),
      color = NA,
      alpha = 0.42,
      extend = right_limit - clade_tip_max
    )
  }

  support_positions <- merge(
    support,
    positions[, c("node", "x", "y"), drop = FALSE],
    by = "node",
    all.x = TRUE,
    sort = FALSE
  )
  if (any(!is.finite(support_positions$x)) || any(!is.finite(support_positions$y))) {
    fail("Not every internal support slot received a plotting coordinate.")
  }
  strong <- support_positions[
    support_positions$display_class == "joint-strong", , drop = FALSE
  ]
  weak <- support_positions[
    support_positions$display_class == "below-joint-threshold", , drop = FALSE
  ]

  plot <- plot +
    ggplot2::geom_point(
      data = strong,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.0,
      stroke = 0.35,
      color = "#303030",
      fill = "#303030"
    ) +
    ggplot2::geom_point(
      data = weak,
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.0,
      stroke = 0.48,
      color = "#6E6E6E",
      fill = "white"
    ) +
    ggplot2::geom_text(
      data = weak,
      ggplot2::aes(x = x, y = y, label = weak_label),
      inherit.aes = FALSE,
      nudge_x = 0.025,
      nudge_y = 0.28,
      hjust = 0,
      vjust = 0.5,
      family = "Arial",
      size = 2.45,
      color = "#555555"
    ) +
    ggtree::geom_tippoint(
      ggplot2::aes(shape = analysis_role, fill = analysis_role),
      color = "#333333",
      size = 2.55,
      stroke = 0.42
    ) +
    ggplot2::geom_segment(
      data = tip_positions,
      ggplot2::aes(x = x, xend = species_label_x - 0.015, y = y, yend = y),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.24,
      color = "#B8B8B8"
    ) +
    ggplot2::geom_text(
      data = tip_positions,
      ggplot2::aes(
        x = species_label_x,
        y = y,
        label = species,
        color = analysis_role,
        fontface = tip_fontface
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      family = "Arial",
      size = 3.2
    ) +
    ggplot2::geom_text(
      data = tip_positions,
      ggplot2::aes(x = accession_label_x, y = y, label = accession),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      family = "Arial",
      size = 2.45,
      color = "#6B6B6B"
    ) +
    ggplot2::geom_text(
      data = band_labels,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      family = "Arial",
      fontface = "bold",
      lineheight = 0.92,
      size = 2.85,
      color = "#52606D"
    ) +
    ggplot2::scale_fill_manual(
      name = "Sequence role",
      values = role_palette,
      breaks = names(role_palette),
      labels = unname(role_labels[names(role_palette)]),
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = "Sequence role",
      values = role_shapes,
      breaks = names(role_shapes),
      labels = unname(role_labels[names(role_shapes)]),
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(study = "#E69F00", expanded = "#262626", outgroup = "#666666"),
      guide = "none"
    ) +
    ggplot2::annotate(
      "segment",
      x = 0, xend = 0.5, y = 0.18, yend = 0.18,
      color = "#555555",
      linewidth = 0.35
    ) +
    ggplot2::annotate(
      "segment",
      x = 0, xend = 0, y = 0.10, yend = 0.26,
      color = "#555555",
      linewidth = 0.35
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.5, xend = 0.5, y = 0.10, yend = 0.26,
      color = "#555555",
      linewidth = 0.35
    ) +
    ggplot2::annotate(
      "text",
      x = 0.25, y = -0.12, label = "0.5 substitutions/site",
      family = "Arial",
      size = 2.65,
      color = "#555555"
    ) +
    ggplot2::expand_limits(x = right_limit, y = -0.35) +
    ggtree::theme_tree() +
    ggplot2::labs(
      title = "BRCA1 protein gene tree",
      subtitle = paste0(
        length(tree$tip.label),
        " vertebrate proteins | IQ-TREE maximum-likelihood phylogram | ",
        "provisional amphibian-outgroup display root"
      ),
      caption = paste(
        "Internal nodes: filled = SH-aLRT >= 80 and UFBoot >= 95;",
        "open = below the joint threshold (exact SH-aLRT/UFBoot shown).",
        "Tip order was not ladderized."
      )
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        override.aes = list(shape = unname(role_shapes), size = 3)
      ),
      shape = "none"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial", color = "#262626", size = 8.5),
      plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 8.2, color = "#555555", hjust = 0),
      plot.caption = ggplot2::element_text(
        size = 7, color = "#555555", hjust = 0, lineheight = 1.12,
        margin = ggplot2::margin(t = 5)
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.title = ggplot2::element_text(size = 7.3, face = "bold"),
      legend.text = ggplot2::element_text(size = 7.1),
      legend.key.height = grid::unit(3.6, "mm"),
      legend.key.width = grid::unit(5.0, "mm"),
      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 8, l = 10),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  output_prefix <- path.expand(options$out_prefix)
  output_directory <- dirname(output_prefix)
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }
  final_paths <- c(
    svg = paste0(output_prefix, ".svg"),
    pdf = paste0(output_prefix, ".pdf"),
    settings = paste0(output_prefix, ".settings.tsv")
  )
  existing <- final_paths[file.exists(final_paths)]
  if (length(existing) > 0L) {
    fail(paste("Refusing to overwrite existing outputs:", paste(existing, collapse = ", ")))
  }
  preview_path <- if (options$preview_png == "") "" else path.expand(options$preview_png)
  if (preview_path != "") {
    if (file.exists(preview_path)) {
      fail(paste("Refusing to overwrite existing preview:", preview_path))
    }
    preview_directory <- dirname(preview_path)
    if (!dir.exists(preview_directory)) {
      dir.create(preview_directory, recursive = TRUE, showWarnings = FALSE)
    }
  }
  temporary_prefix <- tempfile(
    pattern = paste0(".", basename(output_prefix), "."),
    tmpdir = output_directory
  )
  temporary_paths <- c(
    svg = paste0(temporary_prefix, ".svg"),
    pdf = paste0(temporary_prefix, ".pdf"),
    settings = paste0(temporary_prefix, ".settings.tsv")
  )
  on.exit(unlink(temporary_paths, force = TRUE), add = TRUE)

  svglite::svglite(
    filename = temporary_paths[["svg"]],
    width = width,
    height = height,
    bg = "white"
  )
  tryCatch(print(plot), finally = grDevices::dev.off())
  grDevices::cairo_pdf(
    filename = temporary_paths[["pdf"]],
    width = width,
    height = height,
    family = "Arial",
    bg = "white"
  )
  tryCatch(print(plot), finally = grDevices::dev.off())
  if (preview_path != "") {
    ggplot2::ggsave(
      preview_path,
      plot = plot,
      device = "png",
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
  settings <- c(
    renderer = "ggtree+ggplot2 BRCA1 README hero",
    display_style = "compact-rectangular-phylogram",
    root_state = "provisional-outgroup-rooted-display-derivative",
    branch_length = "substitutions-per-site",
    support_format = "sh-alrt/ufboot",
    support_display = "filled-joint-threshold_open-exact-below-threshold",
    sh_alrt_joint_threshold = 80,
    ufboot_joint_threshold = 95,
    joint_strong_node_count = sum(support$display_class == "joint-strong"),
    below_joint_threshold_node_count = sum(
      support$display_class == "below-joint-threshold"
    ),
    undisplayed_internal_label_slot_count = sum(
      support$display_class == "not-displayed"
    ),
    ladderized = "false",
    tip_count = length(tree$tip.label),
    study_tip_count = sum(metadata$analysis_role == "study"),
    expanded_tip_count = sum(metadata$analysis_role == "expanded"),
    outgroup_tip_count = sum(metadata$analysis_role == "outgroup"),
    label_mode = "aligned-italic-species_plus_gray-accession",
    width_inches = width,
    height_inches = height,
    font_family = "Arial",
    species_label_size_mm = 3.2,
    accession_label_size_mm = 2.45,
    support_label_size_mm = 2.45,
    broad_clade_bands = paste(broad_clades$clade_id, collapse = ","),
    Mammalia_band_color = band_palette[["Mammalia"]],
    Sauropsida_band_color = band_palette[["Sauropsida"]],
    Amphibia_band_color = band_palette[["Amphibia"]],
    study_color = role_palette[["study"]],
    expanded_color = role_palette[["expanded"]],
    outgroup_color = role_palette[["outgroup"]],
    tree_sha256 = file_sha256(options$tree),
    metadata_sha256 = file_sha256(options$metadata),
    clade_support_sha256 = file_sha256(options$clade_support),
    R = as.character(getRversion()),
    ape = as.character(utils::packageVersion("ape")),
    ggplot2 = as.character(utils::packageVersion("ggplot2")),
    ggtree = as.character(utils::packageVersion("ggtree")),
    openssl = as.character(utils::packageVersion("openssl")),
    svglite = as.character(utils::packageVersion("svglite"))
  )
  write_settings(temporary_paths[["settings"]], settings)

  finalized <- character()
  for (name in names(final_paths)) {
    if (!file.rename(temporary_paths[[name]], final_paths[[name]])) {
      unlink(finalized, force = TRUE)
      fail(paste("Failed to finalize output:", final_paths[[name]]))
    }
    finalized <- c(finalized, final_paths[[name]])
  }
  message(paste("Wrote", paste(final_paths, collapse = ", ")))
}

tryCatch(
  main(),
  error = function(error) {
    message(paste0("ERROR [BRCA1_README_FIGURE] ", conditionMessage(error)))
    quit(status = 2L, save = "no")
  }
)
