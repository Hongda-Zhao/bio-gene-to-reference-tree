#!/usr/bin/env Rscript

# Render a reviewed gene tree with ggtree + ggplot2 without network access.
# The script performs exact tip-ID joins, never reroots or ladderizes the tree,
# refuses overwrite, and emits vector SVG/PDF plus a settings audit table.

fail <- function(message) {
  stop(message, call. = FALSE)
}

usage_text <- function() {
  paste(
    "Usage: render_tree_ggtree.R --tree <tree.nwk> --metadata <sequence_metadata.tsv>",
    "--out-prefix <new-prefix> --root-state <unrooted|outgroup-rooted>",
    "[--itol-roles <itol_roles.txt>] [--layout <rectangular|circular>]",
    "[--branch-length <auto|phylogram|cladogram>]",
    "[--support-format <none|fasttree-sh-like|sh-alrt/ufboot|sh-alrt/bootstrap>]",
    "[--show-tip-labels <true|false>] [--width <inches>] [--height <inches>]"
  )
}

parse_cli <- function(arguments) {
  if (length(arguments) == 1L && arguments[[1L]] %in% c("--help", "-h")) {
    cat(usage_text(), "\n")
    quit(status = 0L, save = "no")
  }
  values <- list(
    layout = "rectangular",
    branch_length = "auto",
    support_format = "none",
    show_tip_labels = "true",
    width = "10",
    height = "8",
    itol_roles = ""
  )
  aliases <- c(
    "--tree" = "tree",
    "--metadata" = "metadata",
    "--itol-roles" = "itol_roles",
    "--out-prefix" = "out_prefix",
    "--root-state" = "root_state",
    "--layout" = "layout",
    "--branch-length" = "branch_length",
    "--support-format" = "support_format",
    "--show-tip-labels" = "show_tip_labels",
    "--width" = "width",
    "--height" = "height"
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
  required <- c("tree", "metadata", "out_prefix", "root_state")
  missing <- required[vapply(required, function(key) is.null(values[[key]]) || values[[key]] == "", logical(1))]
  if (length(missing) > 0L) {
    fail(paste("Missing required arguments:", paste(missing, collapse = ", ")))
  }
  values
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    fail(paste(
      "Missing local R packages:", paste(missing, collapse = ", "),
      "Install versions compatible with the active R/Bioconductor release before rerunning."
    ))
  }
}

parse_boolean <- function(value, name) {
  normalized <- tolower(value)
  if (!normalized %in% c("true", "false")) {
    fail(paste(name, "must be true or false."))
  }
  normalized == "true"
}

read_selected_metadata <- function(path) {
  metadata <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = character()
  )
  required <- c("tip_id", "analysis_role", "inclusion_status", "species", "accession")
  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0L) {
    fail(paste("Metadata is missing columns:", paste(missing, collapse = ", ")))
  }
  selected <- metadata[metadata$inclusion_status == "selected", , drop = FALSE]
  if (nrow(selected) == 0L) {
    fail("Metadata contains no selected tips.")
  }
  if (
    any(selected$tip_id == "") ||
      any(selected$tip_id != trimws(selected$tip_id)) ||
      anyDuplicated(selected$tip_id)
  ) {
    fail("Selected metadata tip_id values must be non-empty, whitespace-exact, and unique.")
  }
  allowed_roles <- c("study", "expanded", "outgroup")
  invalid_roles <- setdiff(unique(selected$analysis_role), allowed_roles)
  if (length(invalid_roles) > 0L) {
    fail(paste("Unsupported analysis_role values:", paste(invalid_roles, collapse = ", ")))
  }
  selected$display_label <- ifelse(
    selected$species == "",
    selected$tip_id,
    paste0(selected$species, " | ", selected$tip_id)
  )
  # ggtree's %<+% operator treats the first metadata column as the join key.
  # Make that key tip_id explicitly so the contract does not depend on the
  # user's TSV column order.
  selected[, c("tip_id", setdiff(colnames(selected), "tip_id")), drop = FALSE]
}

parse_itol_roles <- function(path, expected_tips, expected_roles) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) < 3L || lines[[1L]] != "DATASET_COLORSTRIP" || !"SEPARATOR TAB" %in% lines) {
    fail("iTOL roles must use the official DATASET_COLORSTRIP format with TAB separation.")
  }
  data_index <- match("DATA", lines)
  if (is.na(data_index) || data_index == length(lines)) {
    fail("iTOL roles file has no DATA rows.")
  }
  data_lines <- lines[(data_index + 1L):length(lines)]
  data_lines <- data_lines[nzchar(data_lines)]
  pieces <- strsplit(data_lines, "\t", fixed = TRUE)
  if (any(lengths(pieces) != 3L)) {
    fail("Every iTOL role row must contain exactly tip, color, and role label.")
  }
  roles <- data.frame(
    tip_id = vapply(pieces, `[[`, character(1), 1L),
    color = vapply(pieces, `[[`, character(1), 2L),
    role_label = vapply(pieces, `[[`, character(1), 3L),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(roles$tip_id) || !identical(sort(roles$tip_id), sort(expected_tips))) {
    fail("iTOL DATA tips must equal the selected metadata and Newick tip set exactly.")
  }
  if (any(!grepl("^#[0-9A-Fa-f]{6}$", roles$color))) {
    fail("Every iTOL role color must be a six-digit hexadecimal value.")
  }
  label_to_role <- c("Study" = "study", "Expanded" = "expanded", "Outgroup" = "outgroup")
  mapped <- unname(label_to_role[roles$role_label])
  if (any(is.na(mapped))) {
    fail("iTOL role labels must be Study, Expanded, or Outgroup.")
  }
  expected <- unname(expected_roles[match(roles$tip_id, names(expected_roles))])
  if (!identical(mapped, expected)) {
    fail("iTOL role labels disagree with sequence_metadata.tsv analysis_role values.")
  }
  colors_by_role <- split(roles$color, mapped)
  if (any(vapply(colors_by_role, function(colors) length(unique(toupper(colors))) != 1L, logical(1)))) {
    fail("Each analysis role must use one consistent iTOL color.")
  }
  palette <- vapply(colors_by_role, function(colors) unique(colors)[[1L]], character(1))
  default_palette <- c(study = "#E69F00", expanded = "#009E73", outgroup = "#999999")
  default_palette[names(palette)] <- palette
  present_roles <- unique(mapped)
  if (anyDuplicated(toupper(unname(default_palette[present_roles])))) {
    fail("Present analysis roles must use visually distinct colors.")
  }
  default_palette
}

file_sha256 <- function(path) {
  requireNamespace("openssl", quietly = TRUE)
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  unname(as.character(openssl::sha256(connection)))
}

descendant_tip_labels <- function(tree, node) {
  tip_count <- length(tree$tip.label)
  if (node <= tip_count) {
    return(tree$tip.label[[node]])
  }
  children <- tree$edge[tree$edge[, 1L] == node, 2L]
  if (length(children) == 0L) {
    fail(paste("Internal node has no descendants:", node))
  }
  unlist(lapply(children, function(child) descendant_tip_labels(tree, child)), use.names = FALSE)
}

validate_outgroup_root_split <- function(tree, outgroup_tips) {
  candidate_roots <- setdiff(unique(tree$edge[, 1L]), unique(tree$edge[, 2L]))
  if (length(candidate_roots) != 1L) {
    fail("A rooted tree must contain exactly one structural root node.")
  }
  root_children <- tree$edge[tree$edge[, 1L] == candidate_roots[[1L]], 2L]
  child_tip_sets <- lapply(root_children, function(child) sort(descendant_tip_labels(tree, child)))
  if (!any(vapply(child_tip_sets, identical, logical(1), sort(outgroup_tips)))) {
    fail("The structural root split does not isolate the selected outgroup tip set.")
  }
}

structural_root_marker_index <- function(tree) {
  if (!ape::is.rooted(tree) || is.null(tree$node.label)) {
    return(NA_integer_)
  }
  candidate_roots <- setdiff(unique(tree$edge[, 1L]), unique(tree$edge[, 2L]))
  if (length(candidate_roots) != 1L) {
    return(NA_integer_)
  }
  label_index <- candidate_roots[[1L]] - length(tree$tip.label)
  if (
    label_index < 1L ||
      label_index > length(tree$node.label) ||
      is.na(tree$node.label[[label_index]]) ||
      !identical(tree$node.label[[label_index]], "Root")
  ) {
    return(NA_integer_)
  }
  as.integer(label_index)
}

validate_support <- function(tree, support_format) {
  allowed <- c("none", "fasttree-sh-like", "sh-alrt/ufboot", "sh-alrt/bootstrap")
  if (!support_format %in% allowed) {
    fail(paste("support-format must be one of", paste(allowed, collapse = ", ")))
  }
  if (support_format == "none") {
    return(FALSE)
  }
  labels <- tree$node.label
  root_marker_index <- structural_root_marker_index(tree)
  if (!is.na(root_marker_index)) {
    # postprocess_brca1_tree.R writes this exact marker on the actual
    # structural root. It records rooting provenance and is not support.
    labels[[root_marker_index]] <- NA_character_
  }
  labels <- labels[!is.na(labels) & labels != ""]
  if (length(labels) == 0L) {
    fail("Support display was requested but the tree has no internal-node labels.")
  }
  if (support_format == "fasttree-sh-like") {
    numeric_labels <- suppressWarnings(as.numeric(labels))
    if (any(is.na(numeric_labels)) || any(numeric_labels < 0 | numeric_labels > 1)) {
      fail("FastTree SH-like support labels must be numeric values from 0 to 1.")
    }
  } else {
    pieces <- strsplit(labels, "/", fixed = TRUE)
    if (any(lengths(pieces) != 2L)) {
      fail("Dual IQ-TREE support labels must have the explicit value/value form.")
    }
    numeric_labels <- suppressWarnings(as.numeric(unlist(pieces, use.names = FALSE)))
    if (any(is.na(numeric_labels)) || any(numeric_labels < 0 | numeric_labels > 100)) {
      fail("IQ-TREE support values must be numeric values from 0 to 100.")
    }
  }
  TRUE
}

choose_branch_mode <- function(tree, requested) {
  if (!requested %in% c("auto", "phylogram", "cladogram")) {
    fail("branch-length must be auto, phylogram, or cladogram.")
  }
  lengths <- tree$edge.length
  if (!is.null(lengths) && any(!is.finite(lengths) | lengths < 0)) {
    fail("Branch lengths must be finite and non-negative when present.")
  }
  valid <- !is.null(lengths) && length(lengths) == nrow(tree$edge) &&
    all(is.finite(lengths)) && all(lengths >= 0) && any(lengths > 0)
  if (requested == "phylogram" && !valid) {
    fail("A phylogram requires complete, finite, non-negative, non-zero branch lengths.")
  }
  if (requested == "auto") {
    return(if (valid) "phylogram" else "cladogram")
  }
  requested
}

write_settings <- function(path, settings) {
  table <- data.frame(
    setting = names(settings),
    value = unname(vapply(settings, as.character, character(1))),
    stringsAsFactors = FALSE
  )
  utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

main <- function() {
  options <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (!options$layout %in% c("rectangular", "circular")) {
    fail("layout must be rectangular or circular.")
  }
  if (!options$root_state %in% c("unrooted", "outgroup-rooted")) {
    fail("root-state must be unrooted or outgroup-rooted; this renderer never reroots a tree.")
  }
  show_tip_labels <- parse_boolean(options$show_tip_labels, "show-tip-labels")
  width <- suppressWarnings(as.numeric(options$width))
  height <- suppressWarnings(as.numeric(options$height))
  if (!is.finite(width) || !is.finite(height) || width <= 0 || height <= 0) {
    fail("width and height must be positive finite numbers in inches.")
  }

  require_packages(c("ape", "ggplot2", "ggtree", "openssl", "svglite"))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggtree))

  tree <- ape::read.tree(options$tree)
  if (inherits(tree, "multiPhylo")) {
    fail("Newick input must contain exactly one tree.")
  }
  if (
    is.null(tree) ||
      length(tree$tip.label) == 0L ||
      any(tree$tip.label == "") ||
      anyDuplicated(tree$tip.label)
  ) {
    fail("Newick must contain non-empty, unique tip labels.")
  }
  if (options$root_state == "outgroup-rooted" && !ape::is.rooted(tree)) {
    fail("root-state outgroup-rooted requires a structurally rooted Newick tree.")
  }
  metadata <- read_selected_metadata(options$metadata)
  if (!identical(sort(tree$tip.label), sort(metadata$tip_id))) {
    missing_metadata <- setdiff(tree$tip.label, metadata$tip_id)
    missing_tree <- setdiff(metadata$tip_id, tree$tip.label)
    fail(paste0(
      "Newick and selected metadata tip sets differ. Missing metadata: ",
      paste(missing_metadata, collapse = ","), "; missing tree tips: ",
      paste(missing_tree, collapse = ",")
    ))
  }
  roles_by_tip <- stats::setNames(metadata$analysis_role, metadata$tip_id)
  if (options$root_state == "outgroup-rooted" && !any(metadata$analysis_role == "outgroup")) {
    fail("root-state outgroup-rooted requires at least one selected outgroup tip.")
  }
  if (options$root_state == "outgroup-rooted") {
    validate_outgroup_root_split(
      tree,
      metadata$tip_id[metadata$analysis_role == "outgroup"]
    )
  }
  palette <- c(study = "#E69F00", expanded = "#009E73", outgroup = "#999999")
  palette_source <- "default"
  if (options$itol_roles != "") {
    palette <- parse_itol_roles(options$itol_roles, tree$tip.label, roles_by_tip)
    palette_source <- "itol_roles.txt"
  }
  branch_mode <- choose_branch_mode(tree, options$branch_length)
  show_support <- validate_support(tree, options$support_format)
  root_marker_index <- structural_root_marker_index(tree)
  if (show_support && !is.na(root_marker_index)) {
    # Do not draw the rooting-provenance marker as if it were branch support.
    tree$node.label[[root_marker_index]] <- NA_character_
  }
  branch_argument <- if (branch_mode == "phylogram") "branch.length" else "none"
  if (length(tree$tip.label) > 150L && show_tip_labels) {
    fail("More than 150 tips require --show-tip-labels false or a deliberately specialized renderer.")
  }
  if (length(tree$tip.label) > 50L && show_tip_labels && identical(options$height, "8")) {
    warning("More than 50 labelled tips are being rendered at the default height; increase --height after visual review.")
  }

  plot <- ggtree::ggtree(
    tree,
    layout = options$layout,
    ladderize = FALSE,
    branch.length = branch_argument
  ) %<+% metadata
  plot <- plot +
    ggtree::geom_tippoint(ggplot2::aes(color = analysis_role), size = 2) +
    ggplot2::scale_color_manual(
      values = palette,
      breaks = c("study", "expanded", "outgroup"),
      drop = FALSE,
      na.translate = FALSE
    ) +
    ggplot2::labs(
      color = "Sequence role",
      x = if (branch_mode == "phylogram") "Substitutions per site" else NULL,
      caption = paste0(
        "Gene tree; root state declared as ", options$root_state,
        if (options$root_state == "unrooted") "; display orientation/root position is arbitrary" else "",
        ". Branches: ",
        if (branch_mode == "phylogram") "substitutions per site" else "topology only",
        ". Tip order was not ladderized. Support: ", options$support_format, "."
      )
    ) +
    ggplot2::theme(legend.position = "right", plot.caption = ggplot2::element_text(hjust = 0))
  if (show_tip_labels) {
    plot <- plot + ggtree::geom_tiplab(ggplot2::aes(label = display_label), size = 2.5)
    if (options$layout == "rectangular") {
      finite_x <- plot$data$x[is.finite(plot$data$x)]
      if (length(finite_x) > 0L) {
        x_span <- diff(range(finite_x))
        if (is.finite(x_span) && x_span > 0) {
          # Reserve deterministic room for species + accession labels on the
          # longest tip. Canvas width alone does not expand the data coordinate
          # range and can still leave long labels clipped at the panel edge.
          plot <- plot + ggplot2::expand_limits(x = max(finite_x) + 0.30 * x_span)
        }
      }
    }
  }
  if (show_support) {
    plot <- plot + ggtree::geom_text2(
      ggplot2::aes(subset = !isTip & !is.na(label) & label != "", label = label),
      size = 2.2,
      hjust = -0.15
    )
  }
  if (branch_mode == "phylogram" && options$layout == "rectangular") {
    plot <- plot + ggtree::theme_tree2()
  }
  if (branch_mode == "phylogram" && options$layout != "rectangular") {
    plot <- plot + ggtree::geom_treescale()
  }

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
  temporary_prefix <- tempfile(pattern = paste0(".", basename(output_prefix), "."), tmpdir = output_directory)
  temporary_paths <- c(
    svg = paste0(temporary_prefix, ".svg"),
    pdf = paste0(temporary_prefix, ".pdf"),
    settings = paste0(temporary_prefix, ".settings.tsv")
  )
  on.exit(unlink(temporary_paths, force = TRUE), add = TRUE)

  ggplot2::ggsave(
    temporary_paths[["svg"]], plot = plot, device = svglite::svglite,
    width = width, height = height, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    temporary_paths[["pdf"]], plot = plot, device = "pdf",
    width = width, height = height, units = "in", bg = "white"
  )
  settings <- c(
    renderer = "ggtree+ggplot2",
    root_state = options$root_state,
    layout = options$layout,
    requested_branch_length = options$branch_length,
    rendered_branch_mode = branch_mode,
    support_format = options$support_format,
    ladderized = "false",
    show_tip_labels = tolower(options$show_tip_labels),
    tip_count = length(tree$tip.label),
    internal_node_count = tree$Nnode,
    width_inches = width,
    height_inches = height,
    palette_source = palette_source,
    study_color = palette[["study"]],
    expanded_color = palette[["expanded"]],
    outgroup_color = palette[["outgroup"]],
    tree_sha256 = file_sha256(options$tree),
    metadata_sha256 = file_sha256(options$metadata),
    itol_roles_sha256 = if (options$itol_roles == "") "not-supplied" else file_sha256(options$itol_roles),
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
    message(paste0("ERROR [GGTREE_RENDER] ", conditionMessage(error)))
    quit(status = 2L, save = "no")
  }
)
