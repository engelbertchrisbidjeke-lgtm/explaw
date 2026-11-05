prepare_data <- function(data, time_col, event_col, group_col = NULL) {
  if (!is.data.frame(data)) stop("'data' must be a data.frame.")
  if (!time_col %in% names(data)) stop("Column '", time_col, "' not found.")
  if (!event_col %in% names(data)) stop("Column '", event_col, "' not found.")
  if (!is.null(group_col) && !group_col %in% names(data)) {
    stop("Column '", group_col, "' not found.")
  }
  if (time_col == event_col) {
    stop("'time_col' and 'event_col' must be different.")
  }
  if (!is.null(group_col) &&
      (group_col == time_col || group_col == event_col)) {
    stop("'group_col' must differ from 'time_col' and 'event_col'.")
  }
  if (!is.numeric(data[[time_col]])) {
    stop("'", time_col, "' must be numeric.")
  }

  event_vals <- unique(data[[event_col]])
  event_vals <- event_vals[!is.na(event_vals)]
  if (!all(event_vals %in% c(0, 1))) {
    stop("'", event_col, "' must contain only 0 and 1 (NA allowed).")
  }

  d <- data
  d$time <- as.numeric(data[[time_col]])
  d$event <- as.integer(data[[event_col]])
  if (time_col != "time") d[[time_col]] <- NULL
  if (event_col != "event") d[[event_col]] <- NULL

  if (!is.null(group_col)) {
    d$group <- data[[group_col]]
    if (group_col != "group") d[[group_col]] <- NULL
  }

  d <- d[!is.na(d$time) & !is.na(d$event) & d$time >= 0, , drop = FALSE]

  others <- setdiff(names(d), c("time", "event", "group"))
  order_cols <- c(intersect(c("time", "event", "group"), names(d)), others)
  d <- d[, order_cols, drop = FALSE]
  rownames(d) <- NULL
  d
}

summarize_data <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  data.frame(
    n_total     = nrow(d),
    n_events    = sum(d$event == 1L),
    n_censored  = sum(d$event == 0L),
    pct_censored = round(100 * mean(d$event == 0L), 1),
    time_min    = min(d$time),
    time_max    = max(d$time)
  )
}
