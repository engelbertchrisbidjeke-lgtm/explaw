#' Prepare raw data for survival analysis
#'
#' Standardizes a raw data frame so that all other functions in the package
#' can use a common format. Renames the time, event and (optional) group
#' columns to `time`, `event` and `group`, checks their types, removes
#' invalid rows (missing, infinite or negative times, missing events), and
#' puts the key columns first.
#'
#' @param data A `data.frame` containing at least a time column and an
#'   event column.
#' @param time_col Character. Name of the column holding the observed
#'   durations. Must be numeric and non-negative. Rows with a missing,
#'   infinite or negative duration are removed.
#' @param event_col Character. Name of the column holding the event
#'   indicator. Must contain only `0` (censored) and `1` (event observed);
#'   numeric, logical, factor or character columns whose values are `0`
#'   and `1` are accepted (a factor is read through its labels, not its
#'   internal codes). `NA` values are tolerated and removed.
#' @param group_col Character or `NULL` (default). Optional name of a
#'   grouping column (e.g. production lot, test condition). When supplied,
#'   it is renamed to `group`.
#'
#' @return A `data.frame` with at least the columns `time` and `event`
#'   (and `group` if `group_col` is provided), followed by all other
#'   original columns. Rows with missing or infinite `time`, negative
#'   `time`, or missing `event` are removed.
#'
#' @details
#' If `data` already contains a column named `time`, `event` (or `group`
#' when `group_col` is supplied) that is not the one you designated, it is
#' replaced by the standardized column and a warning is issued, so that no
#' column is overwritten silently.
#'
#' The function stops with an error if the inputs are invalid (missing
#' column, non-numeric time, `event` outside `{0, 1}`) or if no valid row
#' remains after cleaning.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(
#'   data      = motors,
#'   time_col  = "time",
#'   event_col = "cens",
#'   group_col = "temp"
#' )
#' head(d)
#'
#' @seealso [summarize_data()], [estimate_lambda()]
#' @export
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

  # Event indicator: factors are read through their labels, never their codes
  ev_raw <- data[[event_col]]
  if (is.factor(ev_raw)) ev_raw <- as.character(ev_raw)
  ev <- suppressWarnings(as.numeric(ev_raw))
  if (any(is.na(ev) & !is.na(ev_raw)) || !all(ev[!is.na(ev)] %in% c(0, 1))) {
    stop("'", event_col, "' must contain only 0 and 1 (NA allowed).")
  }

  # Other columns; existing 'time' / 'event' / 'group' columns that were not
  # designated by the user are replaced, with a warning
  others    <- setdiff(names(data), c(time_col, event_col, group_col))
  new_names <- c("time", "event", if (!is.null(group_col)) "group")
  clash     <- intersect(new_names, others)
  if (length(clash) > 0) {
    warning("Column(s) ", paste0("'", clash, "'", collapse = ", "),
            " already present in 'data' and not used as time/event/group: ",
            "replaced by the standardized column(s).")
    others <- setdiff(others, clash)
  }

  # Build the result from scratch (no in-place renaming)
  d <- data.frame(
    time  = as.numeric(data[[time_col]]),
    event = as.integer(ev)
  )
  if (!is.null(group_col)) d$group <- data[[group_col]]
  for (nm in others) d[[nm]] <- data[[nm]]

  keep <- is.finite(d$time) & d$time >= 0 & !is.na(d$event)
  if (!any(keep)) {
    stop("No valid row left after removing missing, infinite or negative ",
         "times and missing events.")
  }
  d <- d[keep, , drop = FALSE]
  rownames(d) <- NULL
  d
}

#' Summarize a prepared dataset
#'
#' Produces a one-row summary of the prepared data frame: total number of
#' units, number of observed events, number of censored units, censoring
#' percentage, and minimum / maximum observed times. Useful as a quick
#' sanity check before fitting any model.
#'
#' @param d A prepared `data.frame`, as returned by [prepare_data()].
#'
#' @return A one-row `data.frame` with columns `n_total`, `n_events`,
#'   `n_censored`, `pct_censored`, `time_min` and `time_max`.
#'
#' @examples
#' data(motors, package = "MASS")
#' d <- prepare_data(motors, "time", "cens", "temp")
#' summarize_data(d)
#'
#' @seealso [prepare_data()]
#' @export
summarize_data <- function(d) {
  if (!is.data.frame(d) || !all(c("time", "event") %in% names(d))) {
    stop("'d' must contain 'time' and 'event' columns. ",
         "Use prepare_data() first.")
  }
  data.frame(
    n_total      = nrow(d),
    n_events     = sum(d$event == 1L),
    n_censored   = sum(d$event == 0L),
    pct_censored = round(100 * mean(d$event == 0L), 1),
    time_min     = min(d$time),
    time_max     = max(d$time)
  )
}
