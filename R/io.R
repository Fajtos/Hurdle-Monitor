.OUT_STATE <- new.env(parent = emptyenv())
.OUT_STATE$root <- "out"

set_out_root <- function(path) {
  .OUT_STATE$root <- path
  invisible(path)
}

get_out_root <- function() {
  .OUT_STATE$root
}

write_table <- function(df, name, caption = NULL, digits = 3) {
  root <- get_out_root()
  dir.create(file.path(root, "tables"), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(root, "tables", paste0(name, ".csv")), row.names = FALSE)
  tex <- knitr::kable(df, format = "latex", booktabs = TRUE,
                      caption = caption, digits = digits)
  writeLines(as.character(tex), file.path(root, "tables", paste0(name, ".tex")))
  invisible(df)
}

save_figure <- function(plot, name, width = 8, height = 4.5) {
  root <- get_out_root()
  dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)
  path <- file.path(root, "figures", paste0(name, ".pdf"))
  
  if (file.exists(path)) try(file.remove(path), silent = TRUE)
  
  
  try_device <- function(dev, ...) {
    tryCatch({
      ggplot2::ggsave(path, plot, width = width, height = height,
                      device = dev, ...)
      TRUE
    }, error = function(e) {
      message("  device failed for ", name, ": ", conditionMessage(e))
      FALSE
    })
  }
  
  if (try_device(grDevices::cairo_pdf)) return(invisible(path))
  if (try_device(grDevices::pdf))       return(invisible(path))
  
  stop("Could not save figure to ", path)
}

save_raw <- function(obj, name) {
  root <- get_out_root()
  dir.create(file.path(root, "raw"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, file.path(root, "raw", paste0(name, ".rds")))
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

safe_median_delay <- function(times, flagged, shift_time, min_post = 5) {
  if (length(times) == 0) return(NA_integer_)
  post_shift <- !is.na(times) & flagged & (times >= shift_time)
  if (sum(post_shift) < min_post) return(NA_integer_)
  as.integer(median(times[post_shift] - shift_time, na.rm = TRUE))
}