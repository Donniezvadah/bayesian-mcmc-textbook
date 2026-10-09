# Shared set-up, sourced (with echo off) at the start of every chapter.
#
# 1. Loads every course function from R/ so that later chapters can use the
#    samplers and diagnostics developed earlier.
# 2. Sets knitr and graphics defaults so figures look consistent.
# Working directory is the project root (execute-dir: project in _quarto.yml).

course_source_all <- function(root = ".") {
  files <- sort(list.files(file.path(root, "R"), pattern = "[.]R$",
                           recursive = TRUE, full.names = TRUE))
  for (f in files) sys.source(f, envir = globalenv())
  invisible(files)
}
course_source_all()

options(digits = 4, scipen = 6, width = 78, warnPartialMatchDollar = TRUE)

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  fig.align = "center",
  fig.width = 7,
  fig.height = 4,
  out.width = "92%",
  dpi = 150,
  dev = "png",
  dev.args = list(bg = "white"),
  course_par = TRUE
)

# Apply the course graphical parameters before every chunk that draws a plot.
knitr::knit_hooks$set(course_par = function(before, options, envir) {
  if (before) course_par()
})
