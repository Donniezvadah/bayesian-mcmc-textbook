# Check internal links in the rendered site.
#
# Usage (from the repository root, after `quarto render`):
#   Rscript scripts/check-links.R [site_dir]
#
# For every HTML file in the site directory, extracts href and src attributes,
# keeps relative links (external URLs, mailto:, javascript: and data: URIs are
# skipped), and checks that (i) the target file exists and (ii) any #fragment
# matches an id in the target page. Exits with status 1 if any link is broken.
# Uses base R only, so it runs in CI before any optional packages are installed.

args <- commandArgs(trailingOnly = TRUE)
site <- if (length(args) >= 1) args[1] else "_book"
if (!dir.exists(site)) stop(sprintf("Site directory '%s' not found; render the book first.", site))

html_files <- list.files(site, pattern = "[.]html$", recursive = TRUE, full.names = TRUE)
read_html <- function(f) paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
ids_cache <- new.env()
page_ids <- function(f) {
  key <- normalizePath(f)
  if (is.null(ids_cache[[key]])) {
    txt <- read_html(f)
    m <- regmatches(txt, gregexpr("\\sid=\"([^\"]+)\"", txt))[[1]]
    ids_cache[[key]] <- sub("^\\sid=\"([^\"]+)\"$", "\\1", m)
  }
  ids_cache[[key]]
}

broken <- character(0); n_checked <- 0L
for (f in html_files) {
  txt <- read_html(f)
  attrs <- regmatches(txt, gregexpr("(href|src)=\"[^\"]*\"", txt))[[1]]
  links <- unique(sub("^(href|src)=\"([^\"]*)\"$", "\\2", attrs))
  links <- links[nzchar(links)]
  links <- links[!grepl("^([a-zA-Z][a-zA-Z0-9+.-]*:|//)", links)]   # skip external and special schemes
  for (l in links) {
    l <- gsub("&amp;", "&", l, fixed = TRUE)
    parts <- strsplit(sub("\\?[^#]*", "", l), "#", fixed = TRUE)[[1]]
    path <- utils::URLdecode(if (length(parts)) parts[1] else "")
    frag <- if (length(parts) > 1) utils::URLdecode(parts[2]) else ""
    target <- if (nzchar(path)) normalizePath(file.path(dirname(f), path), mustWork = FALSE) else normalizePath(f)
    if (dir.exists(target)) target <- file.path(target, "index.html")
    n_checked <- n_checked + 1L
    if (!file.exists(target)) {
      broken <- c(broken, sprintf("%s -> %s (missing file)", sub(paste0("^", site, "/"), "", f), l))
    } else if (nzchar(frag) && grepl("[.]html$", target) && !(frag %in% page_ids(target))) {
      broken <- c(broken, sprintf("%s -> %s (missing anchor)", sub(paste0("^", site, "/"), "", f), l))
    }
  }
}

cat(sprintf("Checked %d internal links in %d HTML files.\n", n_checked, length(html_files)))
if (length(broken)) {
  cat(sprintf("%d broken link(s):\n", length(broken)))
  cat(paste0("  ", unique(broken)), sep = "\n")
  quit(status = 1)
}
cat("No broken internal links.\n")
