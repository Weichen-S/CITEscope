#' Launch the CITEscope Shiny app
#' @param workers Number of parallel workers (default: 4)
#' @param upload_max_mb Max single-file upload size in MB (default: 1024 = 1GB)
#' @param future_max_gb Max future globals size in GB (default: 15)
#' @return shiny.appobj
#' @export
run_app <- function(workers = 4, upload_max_mb = 1024, future_max_gb = 15) {
  if (!requireNamespace("shiny", quietly = TRUE))  stop("Package 'shiny' is required.",  call. = FALSE)
  if (!requireNamespace("future", quietly = TRUE)) stop("Package 'future' is required.", call. = FALSE)

  # session-wide options
  options(
    shiny.maxRequestSize   = upload_max_mb * 1024^2,  # bytes
    future.globals.maxSize = future_max_gb * 1024^3   # bytes
  )

  # defensive cleanup
  try(future::plan(future::sequential), silent = TRUE)
  try(future:::ClusterRegistry("stop"), silent = TRUE)
  gc()

  # Set plan per OS
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::multicore, workers = workers)
  }

  shiny::shinyApp(ui = ui(), server = server)
}
