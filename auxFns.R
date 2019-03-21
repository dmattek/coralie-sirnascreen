# Custom functions to accompany scripts to porcess Coralie's siRNA screen

## Time-series processing ----

# Returns dt with trajectories that last from the first till last frame
# and include at most in.max.break interruptions
LOCtrajExtr = function(in.dt,
                      in.max.break = 1,
                      in.met.series = 'Metadata_Series',
                      in.met.t = 'Metadata_T',
                      in.met.tracklabel = 'TrackObjects_Label',
                      in.met.tracklabeluni = 'TrackObjects_LabelUni',
                      in.aggr.cols = NULL) {
  
  require(data.table)
  
  loc.dt = copy(in.dt)
  
  # TrackObjects assigns the same label for different cells from IdentifyPrimaryObjects
  # The following aggregation makes sure there's a unique TrackObjects_Label
  # for every Site at every time point: it takes the mean intensity of duplicated labels.
  # Make sure it makes sense!!!
  # Roughly 10% of objects affected
  
  if (!is.null(in.aggr.cols)) {
    loc.dt = loc.dt[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                    by = c(in.met.series, in.met.t, in.met.tracklabel), .SDcols = in.aggr.cols]
  }
  
  # cells from different sites have the same TrackObjects_Label
  # make it unique accross the experiment and pad numbers with zeros
  loc.dt[, c(in.met.tracklabeluni) := paste(sprintf("%03d", get(in.met.series)), sprintf("%04d", get(in.met.tracklabel)), sep = "_")]
  
  ####
  ## Allow for single-breaks in the middle of the trajectory
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  
  # build vector with unique timepoints for the entire experiment
  loc.t.range = unique(loc.dt[, c(in.met.t), with = FALSE])
  
  # Select cells with number of timepoints equal to
  # the range nrow(loc.t.range)  AND
  # with first and last frame at the beginning and end of the movie, respectively.
  # If no aggregation columns (in.aggr.cols) are provided,
  # tracks with forks will be omitted here because the number of timepoints exceeds nrow(loc.t.range)
  loc.dt.tmp1 = loc.dt[, .(
    Ntpt = .N,
    T.start = first(get(in.met.t)),
    T.end = last(get(in.met.t))
  ),
  by = c(in.met.tracklabeluni)][Ntpt <=  nrow(loc.t.range) &
                                 T.start == min(loc.t.range) &
                                 T.end == max(loc.t.range)]
  # loc.dt.tmp1 = loc.dt
  
  # With cells selected above,
  # select tcourses with at most 1 break.
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  # Calculation based on differences between consecutive Metadata_T.
  # Metadata_T increases by 1, thus a single-frame break would result in a difference of 2.
  
  # select cells recorded at the beginning and end of the experiment
  loc.dt.tmp2 = loc.dt[get(in.met.tracklabeluni) %in% loc.dt.tmp1[[in.met.tracklabeluni]]]
  
  # set order
  setkeyv(loc.dt.tmp2, c(in.met.series, in.met.tracklabel, in.met.t))
  
  # calculate difference in consecutive time
  loc.dt.tmp2[, Metadata_T.diff := c(-1, diff(get(in.met.t))), by = c(in.met.tracklabeluni)]
  
  # identify cells with at least one break longer than 1 frame
  # column Nframes stores the number of instances with break longer than 1 frame
  loc.dt.tmp2 = loc.dt.tmp2[Metadata_T.diff > 1 + in.max.break, .(Nframes = .N), by = c(in.met.tracklabeluni)]
  
  # Selected trajectories with frames at 1st and last time points AND with at most 1-frame break
  loc.out = loc.dt[get(in.met.tracklabeluni) %in% setdiff(loc.dt.tmp1[[in.met.tracklabeluni]],
                                                       loc.dt.tmp2[[in.met.tracklabeluni]])]
  return(loc.out)
}


# Returns original dt with RealTime column added
# Real time is based on acquisition frequency
# Input parameters:
# in.dt - data.table
# in.met.t - string with the name of a column in in.dt with frame number, e.g. Metadata_T
# in.acq.freq - acquisition frequency in minutes (integer)
LOCaddRealTime = function(in.dt, in.met.t, in.acq.freq, in.rtcol = 'RealTime') {
  
  loc.dt = copy(in.dt)
  
  loc.dt.t.trans = data.table(meta.tmp = unique(loc.dt[[in.met.t]]), 
                              xxx = seq(min(loc.dt[[in.met.t]])*in.acq.freq, max(loc.dt[[in.met.t]])*in.acq.freq, in.acq.freq))

  setnames(loc.dt.t.trans, c('meta.tmp', 'xxx'), c(in.met.t, in.rtcol))

  loc.dt = merge(loc.dt, loc.dt.t.trans, by = in.met.t)
  
  return(loc.dt)
}

## File reading ----
#' Check whether a string consists only from digits
#'
#' @param x Input string
#'
#' @return True if the input string consists only from digits, False otherwise.
#' @export
#'
#' @examples
#' checkDigits('1111')
#' checkDigits('1111cccc')
LOCcheckDigits <- function(x) {
  grepl('^[-]?[0-9]+[.]?[0-9]*$' , x)
}

#' Check whether a string matches TRUE/FALSE, T/F, or T.../F...
#'
#' @param x Input string
#'
#' @return True if the input string matches the pattern, False otherwise.
#' @export
#'
#' @examples
#' checkLogical('TRUE')
#' checkLogical('xxxTxxx')
LOCcheckLogical <- function(x) {
  grepl('^TRUE$|^FALSE$|^T$|^F$' , x)
}



#' Converts string elements of a named list to apporpriate types
#'
#' Strings that consist of digits are converted to type \code{numeric}, strings with TRUE/FALSE, T/F, or T.../F... to \code{logical}.
#'
#' @param in.l Named list fo strings.
#'
#' @return Named list with elements converted to appropriate types.
#' @export
#' @import xlsx
#'
#' @examples
#'  l.tst = list()
#'  l.tst$aaa = '1000'
#'  l.tst$bbb = '1000xx'
#'  l.tst$ccc = 'True'
#'  l.tst$ddd = 'xxxTrue'
#'  l.res = convertStringListToTypes(l.tst)
#'  str(l.res)

LOCconvertStringList2Types <- function(in.l) {
  # convert strings with digits to numeric
  # uses logical indexing: http://stackoverflow.com/questions/42207235/replace-list-elements-by-name-with-another-list
  loc.l = LOCcheckDigits(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.numeric)
  
  # convert strings with TRUE/FALSE to logical
  loc.l = LOCcheckLogical(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.logical)
  
  return(in.l)
}


#' Return a list with parameter names and their values read from xlsx file
#'
#' The Excel xlsx file has to contain at least two columns:
#' 1st column with parameter names
#' 2nd column with parameter values
#'
#' In case of rJava error when loading, run: sudo R CMD javareconf
#'
#' @param in.fname Name of the xlsx file.
#' @param in.cols Vector with column names to read. Has to be of length 2.
#' @param in.sheet.idx Integer with the sheet number in the xlsx to process.
#'
#' @return Named list with parameters and their values.
#' @export

LOCreadPar = function(in.fname, in.cols = 1:2, in.sheet.idx = 1) {
  
  if(length(in.cols) != 2)
    stop('Parameter in.cols has to be of length 2.')
  
  df.par = as.data.frame(readxl::read_xlsx(path = in.fname, 
                                           sheet = in.sheet.idx, 
                                           col_names = F
  ))
  
  # convert data frame with parameters to a list
  l.par = split(df.par[, 2], df.par[, 1])
  
  # convert strings with digits to numeric and strings with TRUE/FALSE to logical
  l.par = LOCconvertStringList2Types(l.par)
  
  return(l.par)
}



## Misc ----
# keep desired number of significant digits in a data.table
LOCsignif_dt <- function(dt, digits) {
  loc.dt = copy(dt)
  
  loc.cols <- vapply(loc.dt, is.double, FUN.VALUE = logical(1))
  loc.cols = names(loc.cols[loc.cols])
  
  loc.dt[, (loc.cols) := signif(.SD, digits), .SDcols = loc.cols]
  
  return(loc.dt)
}

# functions to floor or ceil a number with decimal level
LOCfloor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
LOCceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

#' Put a factor first
#'
#' @param in.dt Input data table
#' @param in.col Name of the column to work on
#' @param in.item String with a factor to put as first
#'
#' @return
#' @export
#'
#' @examples
LOCputItemFirst = function(in.dt, in.col, in.item) {
  loc.v = unique(in.dt[[in.col]])
  loc.v.1st = as.vector(loc.v[loc.v %like% in.item])
  loc.v.other = as.vector(setdiff(loc.v, loc.v.1st))
  
  in.dt[, (in.col) := factor(get(in.col), levels = c(loc.v.1st, loc.v.other))]
}


