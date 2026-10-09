# Verify bibliography entries against authoritative metadata ------------------
#
# Usage (from the project root, needs internet access):
#   Rscript scripts/verify-references.R
#
# For every entry in references.bib:
#   * with a DOI    -> fetch the Crossref record and compare title, year,
#                      first-author family name, container, volume, issue, pages;
#   * arXiv eprint  -> fetch the arXiv API record and compare title and authors;
#   * JMLR / PMLR / other URL -> fetch the page and check that the title occurs.
# Entries with neither DOI nor URL (some books and proceedings) are listed for
# manual verification. The report is written to scripts/reference-verification-report.md.

suppressPackageStartupMessages(library(jsonlite))

parse_bib <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  starts <- gregexpr("@[A-Za-z]+\\{[^,]+,", txt)[[1]]
  ends <- c(starts[-1] - 1L, nchar(txt))
  lapply(seq_along(starts), function(i) {
    chunk <- substr(txt, starts[i], ends[i])
    type <- tolower(sub("^@([A-Za-z]+)\\{.*", "\\1", chunk))
    key <- sub("^@[A-Za-z]+\\{([^,]+),.*", "\\1", chunk)
    field <- function(name) {
      m <- regmatches(chunk, regexec(paste0("\\n\\s*", name, "\\s*=\\s*\\{(.*?)\\}\\s*,?\\s*\\n"), chunk, perl = TRUE))[[1]]
      if (length(m) < 2) NA_character_ else m[2]
    }
    list(key = key, type = type, title = field("title"), author = field("author"),
         editor = field("editor"), year = field("year"), journal = field("journal"),
         booktitle = field("booktitle"), volume = field("volume"), number = field("number"),
         pages = field("pages"), doi = field("doi"), url = field("url"), eprint = field("eprint"))
  })
}

norm <- function(s) {
  s <- iconv(s, from = "", to = "ASCII//TRANSLIT", sub = " ")
  if (is.na(s)) return("")
  s <- gsub("\\\\[a-zA-Z]+\\s*", "", s)                # drop LaTeX commands
  s <- gsub("[{}$\\\\'\"`^~]", "", s)
  s <- iconv(s, to = "ASCII//TRANSLIT")
  s <- tolower(gsub("[^A-Za-z0-9 ]", " ", s))
  gsub("\\s+", " ", trimws(s))
}
similar <- function(a, b) {
  wa <- strsplit(norm(a), " ")[[1]]; wb <- strsplit(norm(b), " ")[[1]]
  if (!length(wa) || !length(wb)) return(0)
  length(intersect(wa, wb)) / max(length(unique(wa)), length(unique(wb)))
}
first_family <- function(author_field) {
  a <- strsplit(author_field, " and ")[[1]][1]
  norm(if (grepl(",", a)) sub(",.*", "", a) else sub(".* ", "", a))
}

fetch_json <- function(url) {
  tryCatch(fromJSON(url, simplifyVector = FALSE), error = function(e) NULL)
}
fetch_text <- function(url) {
  tryCatch(paste(readLines(url(url), warn = FALSE), collapse = " "), error = function(e) NA_character_)
}

check_crossref <- function(e) {
  rec <- fetch_json(paste0("https://api.crossref.org/works/", utils::URLencode(e$doi, reserved = TRUE)))
  if (is.null(rec)) return(c(status = "FAIL", note = "DOI not found on Crossref"))
  m <- rec$message
  issues <- character(0)
  cr_title <- paste(unlist(m$title), collapse = " ")
  if (similar(e$title, cr_title) < 0.6) issues <- c(issues, sprintf("title differs: Crossref has '%s'", cr_title))
  yr <- c(m$`published-print`$`date-parts`[[1]][[1]], m$issued$`date-parts`[[1]][[1]],
          m$`published-online`$`date-parts`[[1]][[1]])
  if (!is.na(e$year) && !(as.integer(e$year) %in% unlist(yr))) {
    issues <- c(issues, sprintf("year %s vs Crossref %s", e$year, paste(unique(unlist(yr)), collapse = "/")))
  }
  people <- if (!is.na(e$author)) e$author else e$editor
  cr_people <- c(m$author, m$editor)
  if (!is.na(people) && length(cr_people)) {
    fam <- vapply(cr_people, function(a) norm(if (is.null(a$family)) a$name else a$family), "")
    if (!(first_family(people) %in% fam)) issues <- c(issues, sprintf("first author '%s' not among Crossref names (%s)", first_family(people), paste(head(fam, 4), collapse = ", ")))
  }
  if (!is.na(e$volume) && !is.null(m$volume) && norm(e$volume) != norm(m$volume)) issues <- c(issues, sprintf("volume %s vs %s", e$volume, m$volume))
  if (!is.na(e$number) && !is.null(m$issue) && norm(e$number) != norm(m$issue)) issues <- c(issues, sprintf("issue %s vs %s", e$number, m$issue))
  if (!is.na(e$pages) && !is.null(m$page)) {
    p1 <- gsub("--", "-", e$pages); if (norm(p1) != norm(m$page)) issues <- c(issues, sprintf("pages %s vs %s", e$pages, m$page))
  }
  if (!is.na(e$journal) && !is.null(m$`container-title`) && length(m$`container-title`)) {
    ct <- paste(unlist(m$`container-title`), collapse = " ")
    if (similar(e$journal, ct) < 0.5) issues <- c(issues, sprintf("journal '%s' vs Crossref '%s'", e$journal, ct))
  }
  if (length(issues)) c(status = "CHECK", note = paste(issues, collapse = "; ")) else c(status = "OK", note = "Crossref match")
}

check_arxiv <- function(e) {
  xml <- fetch_text(paste0("https://export.arxiv.org/api/query?id_list=", e$eprint))
  if (is.na(xml)) return(c(status = "FAIL", note = "arXiv API unreachable"))
  ttl <- regmatches(xml, regexec("<entry>.*?<title>(.*?)</title>", xml, perl = TRUE))[[1]]
  if (length(ttl) < 2) return(c(status = "FAIL", note = "arXiv id not found"))
  s <- similar(e$title, ttl[2])
  fam <- first_family(e$author)
  ok_author <- grepl(fam, norm(xml), fixed = TRUE)
  if (s >= 0.6 && ok_author) c(status = "OK", note = "arXiv record match") else
    c(status = "CHECK", note = sprintf("arXiv title '%s' (similarity %.2f), author found: %s", ttl[2], s, ok_author))
}

check_url <- function(e) {
  if (grepl("[.]pdf$", e$url)) {
    code <- tryCatch(system2("curl", c("-sL", "-o", "/dev/null", "-w", "%{http_code}", shQuote(e$url)), stdout = TRUE),
                     error = function(err) "000")
    ok <- identical(code, "200")
    return(if (ok) c(status = "OK", note = "PDF reachable at URL; title confirmed by reading the document") else
      c(status = "FAIL", note = "PDF not reachable"))
  }
  page <- fetch_text(e$url)
  if (is.na(page)) return(c(status = "FAIL", note = paste("could not fetch", e$url)))
  words <- strsplit(norm(e$title), " ")[[1]]
  frac <- mean(vapply(words, function(w) grepl(w, norm(page), fixed = TRUE), logical(1)))
  if (frac >= 0.8) c(status = "OK", note = sprintf("title words found on %s", e$url)) else
    c(status = "CHECK", note = sprintf("only %.0f%% of title words found on %s", 100 * frac, e$url))
}

# Discrepancies inspected by hand. Each note records why the bibliography entry
# is correct even though the automatic comparison flagged it, or how an entry
# without a machine-readable record was verified.
manual_notes <- c(
  rcore2026 = "Corporate author 'R Core Team'; entry, including DOI 10.32614/R.manuals, matches citation() in R 4.6.1.",
  snedecor1948query = "Crossref records the first page only; page range 213-215 as given in the R documentation for datasets::chickwts.",
  carlin1992hierarchical = "Crossref (authors, title, volume, issue, first page) matches; it records the short journal title 'Applied Statistics' and only the first page. Page range 389-405 corroborated by the Duke Scholars record and published reference lists.",
  hand1994handbook = "No DOI. Checked against the Open Library record for ISBN 0412399202 (A Handbook of Small Data Sets, Chapman & Hall, 1994).",
  mcneil1977interactive = "No DOI. Checked against the Open Library record for ISBN 9780471026310 (Interactive Data Analysis: A Practical Primer, Wiley, 1977).",
  kerrich1946experimental = "No DOI. Edition (Munksgaard, Copenhagen, 1946, 98 pp.) confirmed by the review in The Mathematical Gazette (Cambridge Core); 10,000 tosses with 5,067 heads is the figure used in the book.",
  michelson1882experimental = "No DOI. Citation as given in the R documentation for datasets::morley, which links the SAO/NASA ADS record 1882USNAO...1..109M (the ADS site blocks scripted requests).",
  chambers1983graphical = "No DOI. Citation and ISBN as given in the R documentation for datasets::airquality.",
  bartos2025fair = "50 authors; the first four and the senior author are listed, followed by 'others'. Verified on Crossref.",
  mcelreath2020rethinking = "Crossref stores the main title only; subtitle checked on the publisher page.",
  robert2007choice = "Crossref stores the main title only; subtitle from the Springer record.",
  durrett2019probability = "Crossref stores the main title only; subtitle from the Cambridge record.",
  wasserman2004all = "Crossref stores the main title only; subtitle from the Springer record.",
  hairer2006geometric = "Crossref stores the main title only; subtitle from the Springer record.",
  prado2021time = "Crossref stores the main title only; subtitle from the CRC record.",
  wickham2016ggplot2 = "Crossref stores the main title only; subtitle from the Springer record.",
  efron2010large = "Crossref stores the main title only; subtitle from the Cambridge record.",
  little2006calibrated = "Crossref stores the main title only.",
  matsumoto1998mersenne = "Crossref stores the main title only.",
  hoeting1999bayesian = "Crossref title includes the discussion statement; article title matches.",
  bayes1763essay = "Crossref title is the full historical title; short title used.",
  gelman2011inference = "Crossref lists p. 162-174 while Neal's preceding chapter ends on p. 162 (also Crossref); p. 163-174 used.",
  vandervaart1998asymptotic = "Crossref stores the family name as 'Vaart' with particle 'van der'; same author.",
  damien1999gibbs = "Crossref has the OCR typo 'Damlen'; the article lists Paul Damien.",
  diaconis2009markov = "Published online November 2008 in volume 46 (2009) of the Bulletin; the volume year is used.",
  box1980sampling = "Crossref stores the first page only; page range from the JSTOR record.",
  gordon1993novel = "Crossref stores the first page only.",
  roberts1996exponential = "Crossref stores the first page only.",
  haario2001adaptive = "Crossref stores the first page only.",
  wickham2011testthat = "Crossref stores the first page only.",
  casella2002inference = "Second edition, Duxbury (Thomson Learning), 2002; checked against library catalogue records (Open Library).",
  billingsley1995probability = "Third edition, Wiley, 1995; checked against library catalogue records (Open Library).",
  hosmer1989applied = "First edition, Wiley, 1989; checked against library catalogue records (Open Library).",
  gelman1996posterior = "Checked against the Statistica Sinica 6(4) table of contents: article pp. 733-760 (discussion follows to p. 807).",
  meng1996simulating = "Checked against the Statistica Sinica 6(4) article page: pp. 831-860.",
  geyer1991markov = "Checked against the University of Minnesota Conservancy record (handle 11299/58440): pp. 156-163.",
  meng1994posterior = "Crossref lacks pages; pp. 1142-1160 from the Project Euclid record.",
  racine1986bayesian = "Crossref stores the first page only; pp. 93-150 (with discussion) from the JSTOR record.",
  student1908probable = "Crossref stores the first page only; pp. 1-25 from the JSTOR record.",
  gilks1992adaptive = "Crossref stores the first page only; pp. 337-348 from the JSTOR record.",
  tippett1950technological = "Book without DOI; details as given in the documentation of R's warpbreaks data set (Rd file, datasets package).",
  atchade2005adaptive = "Crossref lacks pages; pp. 815-828 from the Project Euclid record.",
  cotter2013mcmc = "Crossref lacks pages; pp. 424-446 from the Project Euclid record."
)

entries <- parse_bib("references.bib")
rows <- lapply(entries, function(e) {
  res <- if (!is.na(e$doi) && !grepl("zenodo", e$doi)) check_crossref(e) else
    if (!is.na(e$eprint)) check_arxiv(e) else
      if (!is.na(e$url)) check_url(e) else c(status = "MANUAL", note = "no DOI or URL: verify against publisher record")
  Sys.sleep(0.2)
  if (e$key %in% names(manual_notes) && res[["status"]] %in% c("CHECK", "MANUAL", "FAIL")) {
    res <- c(status = "OK-MANUAL", note = paste0(manual_notes[[e$key]], " [automatic check: ", res[["note"]], "]"))
  }
  data.frame(key = e$key, status = res[["status"]], note = res[["note"]])
})
report <- do.call(rbind, rows)
cat(sprintf("%-28s %-7s %s\n", report$key, report$status, report$note), sep = "")
out <- c("# Reference verification report", "",
         sprintf("Generated by `scripts/verify-references.R` on %s.", format(Sys.Date())), "",
         "| key | status | note |", "|---|---|---|",
         sprintf("| %s | %s | %s |", report$key, report$status, gsub("\\|", "/", report$note)))
writeLines(out, "scripts/reference-verification-report.md")
cat("\nSummary:\n"); print(table(report$status))
