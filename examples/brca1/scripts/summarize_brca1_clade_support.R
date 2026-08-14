#!/usr/bin/env Rscript

# Audit named BRCA1 comparison clades in the structurally outgroup-rooted tree.
#
# This is deliberately example-specific.  It accepts exactly the 18 accession
# tips in the public BRCA1 example, validates IQ-TREE's SH-aLRT/UFBoot label
# contract, and writes one row for each predeclared biological comparison.  A
# clade that is absent from the inferred gene tree is reported as such; absence
# is not an execution error and is never converted into a support value.

fail <- function(message) {
  stop(message, call. = FALSE)
}

usage_text <- function() {
  paste(
    "Usage: summarize_brca1_clade_support.R",
    "--tree <gene-tree.outgroup-rooted.nwk> --output <clade_support.tsv>"
  )
}

parse_cli <- function(arguments) {
  if (length(arguments) == 1L && arguments[[1L]] %in% c("--help", "-h")) {
    cat(usage_text(), "\n")
    quit(status = 0L, save = "no")
  }
  values <- list()
  index <- 1L
  while (index <= length(arguments)) {
    flag <- arguments[[index]]
    if (!flag %in% c("--tree", "--output")) {
      fail(paste("Unknown argument:", flag, "\n", usage_text()))
    }
    if (index == length(arguments)) {
      fail(paste("Missing value for", flag))
    }
    key <- switch(flag, "--tree" = "tree", "--output" = "output")
    if (!is.null(values[[key]])) {
      fail(paste("Argument supplied more than once:", flag))
    }
    values[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  missing <- c("tree", "output")[vapply(
    c("tree", "output"),
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

tip_sets <- list(
  Amphibia = c("NP_001107963.1", "XP_029429046.1"),
  Mammalia = c(
    "NP_009225.1", "NP_001038958.1", "NP_033894.3",
    "XP_017204702.2", "NP_001013434.1", "NP_848668.1",
    "XP_014595447.2", "XP_003414318.3", "XP_058140293.1",
    "NP_001029141.1"
  ),
  Sauropsida = c(
    "NP_989500.1", "XP_072775070.1", "XP_019406054.1",
    "XP_023967135.2", "XP_008111382.1", "XP_026576759.1"
  ),
  Primates = c("NP_009225.1", "NP_001038958.1"),
  Glires = c("NP_033894.3", "XP_017204702.2"),
  Aves = c("NP_989500.1", "XP_072775070.1"),
  Lepidosauria = c("XP_008111382.1", "XP_026576759.1"),
  Archelosauria = c(
    "XP_023967135.2", "XP_019406054.1", "NP_989500.1",
    "XP_072775070.1"
  ),
  Archosauria = c("XP_019406054.1", "NP_989500.1", "XP_072775070.1"),
  `Crocodylia+Testudines` = c("XP_019406054.1", "XP_023967135.2")
)

hypothesis_roles <- c(
  Amphibia = "comparison_clade",
  Mammalia = "comparison_clade",
  Sauropsida = "comparison_clade",
  Primates = "comparison_clade",
  Glires = "comparison_clade",
  Aves = "comparison_clade",
  Lepidosauria = "comparison_clade",
  Archelosauria = "comparison_clade",
  Archosauria = "comparison_clade",
  `Crocodylia+Testudines` = "explicit_alternative"
)

expected_tips <- sort(unique(unlist(tip_sets[c("Amphibia", "Mammalia", "Sauropsida")])))

read_one_tree <- function(path) {
  if (!file.exists(path) || file.info(path)$isdir || file.info(path)$size == 0) {
    fail(paste("Tree file is missing, is a directory, or is empty:", path))
  }
  tree <- suppressWarnings(ape::read.tree(file = path))
  if (is.null(tree) || !inherits(tree, "phylo") || inherits(tree, "multiPhylo")) {
    fail("--tree must contain exactly one readable Newick tree.")
  }
  tree
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
  unlist(
    lapply(children, function(child) descendant_tip_labels(tree, child)),
    use.names = FALSE
  )
}

assert_exact_tip_set <- function(observed) {
  missing <- sort(setdiff(expected_tips, observed))
  extra <- sort(setdiff(observed, expected_tips))
  if (length(missing) > 0L || length(extra) > 0L) {
    fail(paste0(
      "Tree tips do not exactly equal the fixed 18-accession BRCA1 set; missing=[",
      paste(missing, collapse = ","), "]; extra=[",
      paste(extra, collapse = ","), "]"
    ))
  }
}

validate_rooted_tree <- function(tree) {
  tip_count <- length(tree$tip.label)
  if (
    any(is.na(tree$tip.label)) || any(tree$tip.label == "") ||
      any(tree$tip.label != trimws(tree$tip.label)) || anyDuplicated(tree$tip.label)
  ) {
    fail("Tip labels must be non-empty, whitespace-exact, and unique.")
  }
  assert_exact_tip_set(tree$tip.label)
  if (!isTRUE(ape::is.rooted(tree))) {
    fail("Tree must be structurally rooted, not merely displayed with an arbitrary Newick root.")
  }
  if (!is.null(tree$root.edge)) {
    fail("Tree must not carry a separate root.edge branch.")
  }
  if (
    is.null(tree$Nnode) || tree$Nnode != tip_count - 1L ||
      is.null(tree$edge) || !is.matrix(tree$edge) || ncol(tree$edge) != 2L ||
      nrow(tree$edge) != 2L * tip_count - 2L
  ) {
    fail("Tree must be a fully resolved rooted binary topology.")
  }
  if (
    is.null(tree$edge.length) || length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)
  ) {
    fail("Every edge must have one finite, non-negative branch length.")
  }
  if (any(!is.finite(tree$edge)) || any(tree$edge != as.integer(tree$edge))) {
    fail("Tree contains non-finite or non-integer node identifiers.")
  }
  if (any(tree$edge[, 1L] == tree$edge[, 2L]) || anyDuplicated(data.frame(tree$edge))) {
    fail("Tree contains a self-edge or duplicated edge.")
  }
  maximum_node <- tip_count + tree$Nnode
  if (any(tree$edge < 1L) || any(tree$edge > maximum_node)) {
    fail("Tree contains a node identifier outside the expected ape range.")
  }
  if (any(tree$edge[, 1L] <= tip_count)) {
    fail("Tree uses a tip as an edge parent.")
  }
  child_counts <- tabulate(tree$edge[, 2L], nbins = maximum_node)
  parent_counts <- tabulate(tree$edge[, 1L], nbins = maximum_node)
  roots <- which(parent_counts > 0L & child_counts == 0L)
  if (length(roots) != 1L || parent_counts[[roots[[1L]]]] != 2L) {
    fail("Tree must have one structural root with exactly two children.")
  }
  internal <- seq.int(tip_count + 1L, maximum_node)
  nonroot_internal <- setdiff(internal, roots[[1L]])
  if (
    any(child_counts[seq_len(tip_count)] != 1L) ||
      any(child_counts[nonroot_internal] != 1L) ||
      any(parent_counts[nonroot_internal] != 2L)
  ) {
    fail("Tree contains disconnected, repeated, or non-binary nodes.")
  }

  root_children <- tree$edge[tree$edge[, 1L] == roots[[1L]], 2L]
  root_sides <- lapply(
    root_children,
    function(child) sort(descendant_tip_labels(tree, child))
  )
  expected_sides <- list(
    sort(tip_sets$Amphibia),
    sort(setdiff(expected_tips, tip_sets$Amphibia))
  )
  observed_keys <- sort(vapply(root_sides, paste, collapse = "\r", character(1)))
  expected_keys <- sort(vapply(expected_sides, paste, collapse = "\r", character(1)))
  if (!identical(observed_keys, expected_keys)) {
    fail("Structural root split must exactly isolate the two approved amphibian outgroups.")
  }
  roots[[1L]]
}

parse_support_pair <- function(label) {
  pattern <- "^([0-9]+(?:\\.[0-9]+)?)/([0-9]+(?:\\.[0-9]+)?)$"
  match <- regexec(pattern, label, perl = TRUE)
  pieces <- regmatches(label, match)[[1L]]
  if (length(pieces) != 3L) {
    return(NULL)
  }
  values <- suppressWarnings(as.numeric(pieces[2:3]))
  if (any(!is.finite(values)) || any(values < 0) || any(values > 100)) {
    return(NULL)
  }
  list(sh_alrt = pieces[[2L]], ufboot = pieces[[3L]])
}

validate_support_labels <- function(tree, root) {
  if (is.null(tree$node.label) || length(tree$node.label) != tree$Nnode) {
    fail("Tree must carry one IQ-TREE label slot for every internal node.")
  }
  internal_nodes <- seq.int(length(tree$tip.label) + 1L, length(tree$tip.label) + tree$Nnode)
  labels <- tree$node.label
  labels[is.na(labels)] <- ""
  root_index <- match(root, internal_nodes)
  if (labels[[root_index]] != "Root") {
    fail("Structural root label must be exactly 'Root'.")
  }
  nonroot <- internal_nodes != root
  blank_nodes <- internal_nodes[nonroot & labels == ""]
  if (length(blank_nodes) > 1L) {
    fail("At most one unlabeled reroot artifact is allowed.")
  }
  if (length(blank_nodes) == 1L) {
    root_children <- tree$edge[tree$edge[, 1L] == root, 2L]
    if (!blank_nodes[[1L]] %in% root_children) {
      fail("The unlabeled reroot artifact must be an immediate child of the structural root.")
    }
    artifact_tips <- sort(descendant_tip_labels(tree, blank_nodes[[1L]]))
    expected_artifact_tips <- sort(setdiff(expected_tips, tip_sets$Amphibia))
    if (!identical(artifact_tips, expected_artifact_tips)) {
      fail("The unlabeled reroot artifact must subtend the 16-tip amniote side of the root split.")
    }
  }
  informative_indices <- which(nonroot & labels != "")
  malformed <- labels[informative_indices][vapply(
    labels[informative_indices],
    function(label) is.null(parse_support_pair(label)),
    logical(1)
  )]
  if (length(malformed) > 0L) {
    fail(paste0(
      "Every informative internal label must be an exact numeric SH-aLRT/UFBoot pair; malformed=[",
      paste(unique(malformed), collapse = ","), "]"
    ))
  }
  list(labels = labels, internal_nodes = internal_nodes, root = root, blank_nodes = blank_nodes)
}

find_exact_clade_node <- function(tree, target_tips) {
  internal_nodes <- seq.int(
    length(tree$tip.label) + 1L,
    length(tree$tip.label) + tree$Nnode
  )
  target <- sort(target_tips)
  matches <- internal_nodes[vapply(
    internal_nodes,
    function(node) identical(sort(descendant_tip_labels(tree, node)), target),
    logical(1)
  )]
  if (length(matches) > 1L) {
    fail(paste("Internal error: more than one node recovered target tip set", paste(target, collapse = ",")))
  }
  if (length(matches) == 0L) NA_integer_ else matches[[1L]]
}

make_summary <- function(tree, label_audit) {
  rows <- lapply(names(tip_sets), function(clade_id) {
    target <- tip_sets[[clade_id]]
    node <- find_exact_clade_node(tree, target)
    if (is.na(node)) {
      return(data.frame(
        clade_id = clade_id,
        hypothesis_role = unname(hypothesis_roles[[clade_id]]),
        target_tip_count = length(target),
        target_tips = paste(target, collapse = ","),
        recovery_status = "not-recovered",
        node_id = "",
        node_label = "",
        support_status = "not_recovered",
        sh_alrt = "",
        ufboot = "",
        stringsAsFactors = FALSE
      ))
    }
    label_index <- match(node, label_audit$internal_nodes)
    label <- label_audit$labels[[label_index]]
    if (node == label_audit$root || node %in% label_audit$blank_nodes) {
      support_status <- "not_applicable_root_artifact"
      sh_alrt <- ""
      ufboot <- ""
    } else {
      support <- parse_support_pair(label)
      if (is.null(support)) {
        fail(paste("Recovered comparison clade lacks an exact SH-aLRT/UFBoot pair:", clade_id))
      }
      support_status <- "available"
      sh_alrt <- support$sh_alrt
      ufboot <- support$ufboot
    }
    data.frame(
      clade_id = clade_id,
      hypothesis_role = unname(hypothesis_roles[[clade_id]]),
      target_tip_count = length(target),
      target_tips = paste(target, collapse = ","),
      recovery_status = "recovered",
      node_id = as.character(node),
      node_label = label,
      support_status = support_status,
      sh_alrt = sh_alrt,
      ufboot = ufboot,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

main <- function() {
  arguments <- parse_cli(commandArgs(trailingOnly = TRUE))
  require_ape()
  if (file.exists(arguments$output)) {
    fail(paste("Refusing to overwrite existing output:", arguments$output))
  }
  tree <- read_one_tree(arguments$tree)
  root <- validate_rooted_tree(tree)
  label_audit <- validate_support_labels(tree, root)
  summary <- make_summary(tree, label_audit)

  output_parent <- dirname(arguments$output)
  if (!dir.exists(output_parent) && !dir.create(output_parent, recursive = TRUE)) {
    fail(paste("Could not create output directory:", output_parent))
  }
  utils::write.table(
    summary,
    file = arguments$output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  cat("Clade-support audit:", arguments$output, "\n")
}

tryCatch(
  main(),
  error = function(condition) {
    cat("ERROR:", conditionMessage(condition), "\n", file = stderr())
    quit(status = 2L, save = "no")
  }
)
