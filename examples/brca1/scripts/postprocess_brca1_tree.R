#!/usr/bin/env Rscript

# Validate, preserve, and root an IQ-TREE gene tree without network access.
#
# The input topology is treated as unrooted.  It is written back without any
# topology-changing operation, while a separate derivative is rooted on the
# two outgroups approved in sequence_metadata.tsv.  Optional exploratory trees
# are compared with the primary topology using the unrooted Robinson-Foulds
# (Penny-Hendy 1985) distance implemented by ape::dist.topo.

fail <- function(message) {
  stop(message, call. = FALSE)
}

usage_text <- function() {
  paste(
    "Usage: postprocess_brca1_tree.R --tree <IQ-TREE.treefile>",
    "--metadata <sequence_metadata.tsv> --out-prefix <new-prefix>",
    "[--profile-tree <LABEL=unrooted-tree.nwk>] ..."
  )
}

parse_cli <- function(arguments) {
  if (length(arguments) == 1L && arguments[[1L]] %in% c("--help", "-h")) {
    cat(usage_text(), "\n")
    quit(status = 0L, save = "no")
  }
  values <- list(profile_trees = character())
  index <- 1L
  while (index <= length(arguments)) {
    flag <- arguments[[index]]
    if (!flag %in% c("--tree", "--metadata", "--out-prefix", "--profile-tree")) {
      fail(paste("Unknown argument:", flag, "\n", usage_text()))
    }
    if (index == length(arguments)) {
      fail(paste("Missing value for", flag))
    }
    value <- arguments[[index + 1L]]
    if (flag == "--profile-tree") {
      values$profile_trees <- c(values$profile_trees, value)
    } else {
      key <- switch(
        flag,
        "--tree" = "tree",
        "--metadata" = "metadata",
        "--out-prefix" = "out_prefix"
      )
      if (!is.null(values[[key]])) {
        fail(paste("Argument supplied more than once:", flag))
      }
      values[[key]] <- value
    }
    index <- index + 2L
  }
  required <- c("tree", "metadata", "out_prefix")
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

require_ape <- function() {
  if (!requireNamespace("ape", quietly = TRUE)) {
    fail("Missing local R package: ape. Install it before rerunning; this script never installs packages.")
  }
}

read_selected_metadata <- function(path) {
  if (!file.exists(path) || file.info(path)$isdir || file.info(path)$size == 0) {
    fail(paste("Metadata file is missing, is a directory, or is empty:", path))
  }
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
  required <- c("tip_id", "analysis_role", "inclusion_status")
  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0L) {
    fail(paste("Metadata is missing columns:", paste(missing, collapse = ", ")))
  }
  if (nrow(metadata) == 0L) {
    fail("Metadata contains no records.")
  }
  if (
    any(is.na(metadata$tip_id)) ||
      any(metadata$tip_id == "") ||
      any(metadata$tip_id != trimws(metadata$tip_id)) ||
      anyDuplicated(metadata$tip_id)
  ) {
    fail("All metadata tip_id values must be non-empty, whitespace-exact, and unique.")
  }
  if (any(!metadata$inclusion_status %in% c("selected", "rejected"))) {
    fail("Metadata inclusion_status values must be selected or rejected.")
  }
  selected <- metadata[metadata$inclusion_status == "selected", , drop = FALSE]
  if (nrow(selected) < 4L) {
    fail("Metadata must contain at least four selected tips.")
  }
  if (any(!selected$analysis_role %in% c("study", "expanded", "outgroup"))) {
    fail("Selected metadata analysis_role values must be study, expanded, or outgroup.")
  }
  outgroups <- selected$tip_id[selected$analysis_role == "outgroup"]
  if (length(outgroups) != 2L) {
    fail(paste("Exactly two selected outgroup-role tips are required; found", length(outgroups)))
  }
  list(selected = selected, outgroups = outgroups)
}

read_one_tree <- function(path, label) {
  if (!file.exists(path) || file.info(path)$isdir || file.info(path)$size == 0) {
    fail(paste(label, "tree file is missing, is a directory, or is empty:", path))
  }
  tree <- suppressWarnings(ape::read.tree(file = path))
  if (is.null(tree) || !inherits(tree, "phylo") || inherits(tree, "multiPhylo")) {
    fail(paste(label, "must contain exactly one readable Newick tree."))
  }
  tree
}

validate_unrooted_binary_tree <- function(tree, label) {
  tip_count <- length(tree$tip.label)
  if (tip_count < 4L) {
    fail(paste(label, "must contain at least four tips."))
  }
  if (
    any(is.na(tree$tip.label)) ||
      any(tree$tip.label == "") ||
      any(tree$tip.label != trimws(tree$tip.label)) ||
      anyDuplicated(tree$tip.label)
  ) {
    fail(paste(label, "tip labels must be non-empty, whitespace-exact, and unique."))
  }
  if (is.null(tree$edge) || !is.matrix(tree$edge) || ncol(tree$edge) != 2L) {
    fail(paste(label, "has no valid two-column edge matrix."))
  }
  if (any(!is.finite(tree$edge)) || any(tree$edge != as.integer(tree$edge))) {
    fail(paste(label, "contains non-finite or non-integer node identifiers."))
  }
  if (!is.null(tree$root.edge)) {
    fail(paste(label, "has a root.edge and is not an untouched unrooted IQ-TREE topology."))
  }
  if (isTRUE(ape::is.rooted(tree))) {
    fail(paste(label, "is structurally rooted; an unrooted binary IQ-TREE topology is required."))
  }
  if (is.null(tree$Nnode) || tree$Nnode != tip_count - 2L) {
    fail(paste(label, "is not a fully resolved unrooted binary topology."))
  }
  if (nrow(tree$edge) != 2L * tip_count - 3L) {
    fail(paste(label, "has an invalid edge count for an unrooted binary topology."))
  }
  if (
    is.null(tree$edge.length) ||
      length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0)
  ) {
    fail(paste(label, "must have one finite, non-negative branch length per edge."))
  }
  if (any(tree$edge[, 1L] == tree$edge[, 2L]) || anyDuplicated(data.frame(tree$edge))) {
    fail(paste(label, "contains a self-edge or duplicated edge."))
  }

  maximum_node <- tip_count + tree$Nnode
  if (any(tree$edge < 1L) || any(tree$edge > maximum_node)) {
    fail(paste(label, "contains a node identifier outside the expected ape range."))
  }
  if (any(tree$edge[, 1L] <= tip_count)) {
    fail(paste(label, "uses a tip as an edge parent."))
  }
  child_counts <- tabulate(tree$edge[, 2L], nbins = maximum_node)
  parent_counts <- tabulate(tree$edge[, 1L], nbins = maximum_node)
  roots <- which(parent_counts > 0L & child_counts == 0L)
  if (length(roots) != 1L) {
    fail(paste(label, "must have exactly one structural Newick root."))
  }
  root <- roots[[1L]]
  if (parent_counts[[root]] != 3L) {
    fail(paste(label, "structural root must have degree three for an unrooted binary tree."))
  }
  internal <- seq.int(tip_count + 1L, maximum_node)
  nonroot_internal <- setdiff(internal, root)
  if (
    any(child_counts[seq_len(tip_count)] != 1L) ||
      any(child_counts[nonroot_internal] != 1L) ||
      any(parent_counts[nonroot_internal] != 2L)
  ) {
    fail(paste(label, "contains disconnected, repeated, or non-binary nodes."))
  }
  invisible(tree)
}

assert_tip_set_equal <- function(observed, expected, label) {
  missing <- sort(setdiff(expected, observed))
  extra <- sort(setdiff(observed, expected))
  if (length(missing) > 0L || length(extra) > 0L) {
    fail(paste0(
      label,
      " tip set does not exactly equal selected metadata; missing=[",
      paste(missing, collapse = ","),
      "]; extra=[",
      paste(extra, collapse = ","),
      "]"
    ))
  }
}

component_tip_labels <- function(tree, start, blocked_from, blocked_to) {
  edges <- tree$edge
  neighbours <- split(
    c(edges[, 2L], edges[, 1L]),
    c(edges[, 1L], edges[, 2L])
  )
  visited <- integer()
  stack <- start
  while (length(stack) > 0L) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (node %in% visited) {
      next
    }
    visited <- c(visited, node)
    adjacent <- neighbours[[as.character(node)]]
    if (is.null(adjacent)) {
      adjacent <- integer()
    }
    for (next_node in adjacent) {
      if (
        (node == blocked_from && next_node == blocked_to) ||
          (node == blocked_to && next_node == blocked_from)
      ) {
        next
      }
      if (!next_node %in% visited) {
        stack <- c(stack, next_node)
      }
    }
  }
  sort(tree$tip.label[visited[visited <= length(tree$tip.label)]])
}

find_isolating_edge <- function(tree, target_tips) {
  target <- sort(target_tips)
  all_tips <- sort(tree$tip.label)
  complement <- sort(setdiff(all_tips, target))
  for (index in seq_len(nrow(tree$edge))) {
    from <- tree$edge[index, 1L]
    to <- tree$edge[index, 2L]
    side <- component_tip_labels(tree, to, from, to)
    if (identical(side, target) || identical(side, complement)) {
      return(index)
    }
  }
  NA_integer_
}

descendant_tip_labels <- function(tree, node) {
  tip_count <- length(tree$tip.label)
  if (node <= tip_count) {
    return(tree$tip.label[[node]])
  }
  children <- tree$edge[tree$edge[, 1L] == node, 2L]
  if (length(children) == 0L) {
    fail(paste("Rooted derivative contains an internal node without descendants:", node))
  }
  unlist(lapply(children, function(child) descendant_tip_labels(tree, child)), use.names = FALSE)
}

encode_tip_set <- function(tips) {
  ordered <- sort(enc2utf8(tips))
  encoded <- vapply(
    ordered,
    function(tip) paste0(nchar(tip, type = "bytes"), ":", tip),
    character(1)
  )
  paste0(length(ordered), "|", paste(encoded, collapse = "|"))
}

canonical_split_key <- function(side, all_tips) {
  side <- sort(side)
  complement <- sort(setdiff(all_tips, side))
  if (
    length(side) == 0L ||
      length(complement) == 0L ||
      length(side) + length(complement) != length(all_tips)
  ) {
    fail("Cannot encode an empty, complete, or inconsistent tree split.")
  }
  side_key <- encode_tip_set(side)
  complement_key <- encode_tip_set(complement)
  if (
    length(side) < length(complement) ||
      (length(side) == length(complement) && side_key < complement_key)
  ) {
    side_key
  } else {
    complement_key
  }
}

structural_root <- function(tree, label) {
  roots <- setdiff(unique(tree$edge[, 1L]), unique(tree$edge[, 2L]))
  if (length(roots) != 1L) {
    fail(paste(label, "must contain exactly one structural Newick root."))
  }
  roots[[1L]]
}

normalized_node_labels <- function(tree, label) {
  if (is.null(tree$node.label)) {
    return(rep(NA_character_, tree$Nnode))
  }
  if (length(tree$node.label) != tree$Nnode) {
    fail(paste(label, "must contain one node-label slot per internal node."))
  }
  labels <- as.character(tree$node.label)
  labels[labels == ""] <- NA_character_
  labels
}

capture_unrooted_split_labels <- function(tree) {
  tip_count <- length(tree$tip.label)
  internal <- seq.int(tip_count + 1L, tip_count + tree$Nnode)
  root <- structural_root(tree, "Primary unrooted tree")
  labels <- normalized_node_labels(tree, "Primary unrooted tree")
  root_label <- labels[[root - tip_count]]
  if (!is.na(root_label)) {
    fail(
      "The arbitrary structural root of the unrooted tree carries a label that cannot be assigned to a bipartition."
    )
  }

  nodes <- setdiff(internal, root)
  keys <- character(length(nodes))
  split_labels <- rep(NA_character_, length(nodes))
  for (index in seq_along(nodes)) {
    node <- nodes[[index]]
    parent_edges <- which(tree$edge[, 2L] == node)
    if (length(parent_edges) != 1L) {
      fail("Every non-root internal node must have exactly one parent edge.")
    }
    edge_index <- parent_edges[[1L]]
    side <- component_tip_labels(
      tree,
      tree$edge[edge_index, 2L],
      tree$edge[edge_index, 1L],
      tree$edge[edge_index, 2L]
    )
    keys[[index]] <- canonical_split_key(side, tree$tip.label)
    split_labels[[index]] <- labels[[node - tip_count]]
  }
  if (length(keys) != tip_count - 3L || anyDuplicated(keys)) {
    fail("Primary unrooted internal edges did not yield one unique key per bipartition.")
  }
  ordering <- order(keys)
  list(keys = keys[ordering], labels = split_labels[ordering])
}

root_child_nodes <- function(tree, outgroups) {
  root <- structural_root(tree, "Outgroup-rooted derivative")
  children <- tree$edge[tree$edge[, 1L] == root, 2L]
  if (length(children) != 2L) {
    fail("The outgroup-rooted derivative root must have exactly two children.")
  }
  child_tip_sets <- lapply(
    children,
    function(child) sort(descendant_tip_labels(tree, child))
  )
  outgroup_matches <- which(vapply(
    child_tip_sets,
    identical,
    logical(1),
    sort(outgroups)
  ))
  complement_matches <- which(vapply(
    child_tip_sets,
    identical,
    logical(1),
    sort(setdiff(tree$tip.label, outgroups))
  ))
  if (length(outgroup_matches) != 1L || length(complement_matches) != 1L) {
    fail("The rooted child nodes do not uniquely identify outgroup and complement sides.")
  }
  list(
    root = root,
    outgroup = children[[outgroup_matches[[1L]]]],
    complement = children[[complement_matches[[1L]]]]
  )
}

capture_rooted_split_labels <- function(tree, outgroups) {
  tip_count <- length(tree$tip.label)
  root_nodes <- root_child_nodes(tree, outgroups)
  labels <- normalized_node_labels(tree, "Outgroup-rooted derivative")
  root_label <- labels[[root_nodes$root - tip_count]]
  if (!identical(root_label, "Root")) {
    fail("The structural root label must be exactly 'Root' after support remapping.")
  }
  if (!is.na(labels[[root_nodes$complement - tip_count]])) {
    fail("The complement-side root child must be the single unlabeled duplicate of the root split.")
  }

  internal <- seq.int(tip_count + 1L, tip_count + tree$Nnode)
  nodes <- setdiff(internal, c(root_nodes$root, root_nodes$complement))
  keys <- vapply(
    nodes,
    function(node) {
      canonical_split_key(descendant_tip_labels(tree, node), tree$tip.label)
    },
    character(1)
  )
  if (length(keys) != tip_count - 3L || anyDuplicated(keys)) {
    fail("Rooted internal nodes did not yield one unique retained label slot per unrooted bipartition.")
  }
  split_labels <- labels[nodes - tip_count]
  ordering <- order(keys)
  list(keys = keys[ordering], labels = split_labels[ordering])
}

assert_split_label_maps_equal <- function(expected, observed) {
  if (!identical(expected$keys, observed$keys)) {
    fail("Rerooting changed the set of internal unrooted bipartitions.")
  }
  expected_missing <- is.na(expected$labels)
  observed_missing <- is.na(observed$labels)
  if (
    !identical(expected_missing, observed_missing) ||
      !identical(expected$labels[!expected_missing], observed$labels[!observed_missing])
  ) {
    fail("Rerooting did not preserve the exact support label assigned to every unrooted bipartition.")
  }
  invisible(TRUE)
}

restore_rooted_split_labels <- function(tree, split_map, outgroups) {
  tip_count <- length(tree$tip.label)
  root_nodes <- root_child_nodes(tree, outgroups)
  internal <- seq.int(tip_count + 1L, tip_count + tree$Nnode)
  # ape::write.tree() serializes NA_character_ node labels as the literal
  # token "NA".  Use an empty slot for genuinely unlabeled edges, then
  # normalize empty slots back to NA only while comparing the support map.
  labels <- rep("", tree$Nnode)
  labels[[root_nodes$root - tip_count]] <- "Root"

  assigned_keys <- character()
  for (node in setdiff(internal, c(root_nodes$root, root_nodes$complement))) {
    key <- canonical_split_key(descendant_tip_labels(tree, node), tree$tip.label)
    match_index <- match(key, split_map$keys)
    if (is.na(match_index)) {
      fail("The rooted derivative contains a bipartition absent from the unrooted support map.")
    }
    if (key %in% assigned_keys) {
      fail("A support label would be assigned more than once after rerooting.")
    }
    mapped_label <- split_map$labels[[match_index]]
    if (!is.na(mapped_label)) {
      labels[[node - tip_count]] <- mapped_label
    }
    assigned_keys <- c(assigned_keys, key)
  }
  if (!identical(sort(assigned_keys), sort(split_map$keys))) {
    fail("Not every unrooted support label was assigned exactly once after rerooting.")
  }
  tree$node.label <- labels
  assert_split_label_maps_equal(
    split_map,
    capture_rooted_split_labels(tree, outgroups)
  )
  tree
}

validate_rooted_derivative <- function(tree, outgroups) {
  tip_count <- length(tree$tip.label)
  if (!isTRUE(ape::is.rooted(tree))) {
    fail("The outgroup-rooted derivative is not structurally rooted.")
  }
  if (is.null(tree$Nnode) || tree$Nnode != tip_count - 1L) {
    fail("The outgroup-rooted derivative is not fully binary.")
  }
  if (
    is.null(tree$edge.length) ||
      length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0)
  ) {
    fail("The outgroup-rooted derivative has missing, non-finite, or negative branch lengths.")
  }
  root <- structural_root(tree, "Outgroup-rooted derivative")
  root_children <- tree$edge[tree$edge[, 1L] == root, 2L]
  if (length(root_children) != 2L) {
    fail("The outgroup-rooted derivative root must have exactly two children.")
  }
  child_tip_sets <- lapply(
    root_children,
    function(child) sort(descendant_tip_labels(tree, child))
  )
  expected_sets <- list(
    sort(outgroups),
    sort(setdiff(tree$tip.label, outgroups))
  )
  observed_keys <- sort(vapply(child_tip_sets, paste, collapse = "\r", character(1)))
  expected_keys <- sort(vapply(expected_sets, paste, collapse = "\r", character(1)))
  if (!identical(observed_keys, expected_keys)) {
    fail("The structural root split does not exactly isolate the approved outgroup set.")
  }
  invisible(tree)
}

parse_profile_specs <- function(specifications) {
  if (length(specifications) == 0L) {
    return(data.frame(label = character(), path = character(), stringsAsFactors = FALSE))
  }
  labels <- character(length(specifications))
  paths <- character(length(specifications))
  for (index in seq_along(specifications)) {
    pieces <- strsplit(specifications[[index]], "=", fixed = TRUE)[[1L]]
    if (length(pieces) < 2L) {
      fail("--profile-tree values must use LABEL=PATH.")
    }
    labels[[index]] <- pieces[[1L]]
    paths[[index]] <- paste(pieces[-1L], collapse = "=")
    if (
      labels[[index]] == "" ||
        !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", labels[[index]]) ||
        paths[[index]] == ""
    ) {
      fail("Profile labels must be non-empty safe identifiers and paths must be non-empty.")
    }
  }
  if (anyDuplicated(labels) || "primary" %in% labels) {
    fail("Profile labels must be unique; primary is reserved.")
  }
  data.frame(label = labels, path = paths, stringsAsFactors = FALSE)
}

make_sensitivity_table <- function(primary, primary_path, profile_specs, outgroups) {
  tip_count <- length(primary$tip.label)
  rf_maximum <- 2L * (tip_count - 3L)
  rows <- list(data.frame(
    comparison_id = "primary",
    reference_tree = primary_path,
    profile_tree = primary_path,
    n_tips = tip_count,
    rf_distance = 0L,
    rf_maximum = rf_maximum,
    rf_normalized = "0.000000",
    topology_identical = "true",
    outgroup_isolating_split = "true",
    root_screen_status = "pass",
    interpretation = "distance 0: identical unrooted topology (self-check)",
    stringsAsFactors = FALSE
  ))
  if (nrow(profile_specs) == 0L) {
    return(rows[[1L]])
  }
  for (index in seq_len(nrow(profile_specs))) {
    label <- profile_specs$label[[index]]
    path <- profile_specs$path[[index]]
    profile <- read_one_tree(path, paste0("Profile '", label, "'"))
    validate_unrooted_binary_tree(profile, paste0("Profile '", label, "'"))
    assert_tip_set_equal(profile$tip.label, primary$tip.label, paste0("Profile '", label, "'"))
    profile_isolating_edge <- find_isolating_edge(profile, outgroups)
    if (is.na(profile_isolating_edge)) {
      fail(paste0(
        "Profile '", label,
        "' does not contain the approved outgroups as an isolating unrooted split; ",
        "no rooted derivative was written. Outgroups=[",
        paste(sort(outgroups), collapse = ","),
        "]"
      ))
    }
    rf <- as.numeric(ape::dist.topo(primary, profile, method = "PH85"))
    if (length(rf) != 1L || !is.finite(rf) || rf < 0 || rf > rf_maximum) {
      fail(paste("Profile", label, "produced an invalid unrooted Robinson-Foulds distance."))
    }
    identical_topology <- rf == 0
    rows[[length(rows) + 1L]] <- data.frame(
      comparison_id = label,
      reference_tree = primary_path,
      profile_tree = path,
      n_tips = tip_count,
      rf_distance = as.integer(rf),
      rf_maximum = rf_maximum,
      rf_normalized = sprintf("%.6f", rf / rf_maximum),
      topology_identical = ifelse(identical_topology, "true", "false"),
      outgroup_isolating_split = "true",
      root_screen_status = "pass",
      interpretation = ifelse(
        identical_topology,
        "distance 0: identical unrooted topology",
        "distance > 0: unrooted topology differs"
      ),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

assert_outputs_are_new <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L) {
    fail(paste("Refusing to overwrite existing output(s):", paste(existing, collapse = ", ")))
  }
}

main <- function() {
  arguments <- parse_cli(commandArgs(trailingOnly = TRUE))
  require_ape()
  metadata <- read_selected_metadata(arguments$metadata)
  primary <- read_one_tree(arguments$tree, "Primary")
  validate_unrooted_binary_tree(primary, "Primary")
  assert_tip_set_equal(primary$tip.label, metadata$selected$tip_id, "Primary")

  isolating_edge <- find_isolating_edge(primary, metadata$outgroups)
  if (is.na(isolating_edge)) {
    fail(paste0(
      "The two approved outgroups are not monophyletic in the unrooted topology; ",
      "no rooted derivative was written. Outgroups=[",
      paste(sort(metadata$outgroups), collapse = ","),
      "]"
    ))
  }

  split_labels <- capture_unrooted_split_labels(primary)
  rooted <- ape::root(primary, outgroup = metadata$outgroups, resolve.root = TRUE)
  assert_tip_set_equal(rooted$tip.label, primary$tip.label, "Outgroup-rooted derivative")
  validate_rooted_derivative(rooted, metadata$outgroups)
  rooted <- restore_rooted_split_labels(rooted, split_labels, metadata$outgroups)

  profile_specs <- parse_profile_specs(arguments$profile_trees)
  sensitivity <- make_sensitivity_table(
    primary,
    arguments$tree,
    profile_specs,
    metadata$outgroups
  )

  output_paths <- c(
    unrooted = paste0(arguments$out_prefix, ".unrooted.nwk"),
    rooted = paste0(arguments$out_prefix, ".outgroup-rooted.nwk"),
    sensitivity = paste0(arguments$out_prefix, ".topology_sensitivity.tsv")
  )
  assert_outputs_are_new(output_paths)
  output_parent <- dirname(arguments$out_prefix)
  if (!dir.exists(output_parent) && !dir.create(output_parent, recursive = TRUE)) {
    fail(paste("Could not create output directory:", output_parent))
  }

  # No topology-changing function has been applied to primary.  write.tree
  # only normalizes its Newick serialization; rooted is a separate derivative.
  ape::write.tree(primary, file = output_paths[["unrooted"]])
  ape::write.tree(rooted, file = output_paths[["rooted"]])
  utils::write.table(
    sensitivity,
    file = output_paths[["sensitivity"]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  cat("Validated unrooted tree:", output_paths[["unrooted"]], "\n")
  cat("Verified outgroup-rooted derivative:", output_paths[["rooted"]], "\n")
  cat("Topology sensitivity:", output_paths[["sensitivity"]], "\n")
}

tryCatch(
  main(),
  error = function(condition) {
    cat("ERROR:", conditionMessage(condition), "\n", file = stderr())
    quit(status = 2L, save = "no")
  }
)
