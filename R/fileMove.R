#' This function is used to move files to other locations.
#'
#' @param from The current path to the file.
#' @param to The path to the new location of the file.
#' @return nothing to return
#' @examples
#' fileMove("~/source", "~/newLocation")

fileMove <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}
