# Process data from Coralie's siRNA screen
# 1. Identifies long trajectories from single-cell time-series data from LAP track (merged, 1-line-header csv.gz)
# 2. Assign track-ids to cells in receptor snaphot (via fuzzy merge)

require(data.table, quietly = T)
require(R.utils, quietly = T)
require(fuzzyjoin, quietly = T)
require(optparse, quietly = T)


## Custom timeseries processing ----

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


## Custom file reading functions ----

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


#' Return a list with parameter names and their values read from xlsx or csv file
#'
#' The xlsx or csv file has to contain at least two columns:
#' 1st column with parameter names
#' 2nd column with parameter values
#'
#' @param in.fname Name of the xlsx file.
#' @param in.sheet.idx Integer with the sheet number in the xlsx to process.
#'
#' @return Named list with parameters and their values.
#' @export

LOCreadPar = function(in.fname, in.sheet.idx = 1) {
  
  if (length(grep('.xlsx', in.fname)) > 0) {
    require(readxl)
    
    df.par = as.data.frame(readxl::read_xlsx(path = in.fname, 
                                             sheet = in.sheet.idx, 
                                             col_names = F))
  } else if( length(grep('.csv', in.fname)) > 0 ) {
    df.par = read.csv(file = in.fname,  header = F, sep = ',', stringsAsFactors = F, strip.white = T)
  } else
    stop('LOCreadPar: file with parameters is neither a csv, nor an xlsx', call. = F)
  
  
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



## Parse command-line arguments ----
# Code taken from:
# https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/

option_list = list(
  optparse::make_option(c("-r", "--rootdir"), type="character", default=NULL, 
                        help="path to root directory", metavar="character"),
  optparse::make_option(c("-f", "--plotformat"), type="character", default=NULL, 
                        help="path to file with parameters of the analysis, csv or xlsx e.g. cp.out/plotFormat.xlsx", metavar="character"),
  optparse::make_option(c("-d", "--opt$debug"), action="store_true", default = FALSE, dest = 'debug', 
                        help="Print extra output")
); 

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

if (is.null(opt$rootdir)){
  optparse::print_help(opt_parser)
  stop("Please supply path to root directory with data to analyse (input rootdir).\n", call.=FALSE)
}

if (is.null(opt$plotformat)){
  optparse::print_help(opt_parser)
  stop("Please supply path to parameter file, csv or xlsx (input plotformat).\n", call.=FALSE)
}

if (opt$debug)
  cat('\nEntering debug mode. More output...\n\n')

## Read parameter file ----
# The file can be either in csv or xlsx format.
cat(sprintf('Reading parameters of the analysis from:\n%s\n\n', path.expand(opt$plotformat)))


#opt$rootdir = '/Volumes/MacintoshHD/tmp/Coralie/20190311_systII_siPOOLs_plate1_and_2_singlePulse/20190311_102308_731/'
#opt$plotformat = 'plotFormat.csv'

#opt$rootdir = '/Volumes/MacintoshHD/tmp/Coralie/20190620_113936_274'
#opt$plotformat = '/Volumes/MacintoshHD/tmp/Coralie/20190620_113936_274/plotFormat.csv'


l.par = LOCreadPar(file.path(opt$plotformat))

l.par$dir.root = opt$rootdir
cat(sprintf('Working in:\n%s\n\n', l.par$dir.root))

## Set other variables ----

l.col = list()

l.col$imerk.cyto.ring.meanint = l.par$imerk.cyto.ring.meanint
l.col$imerk.nuc.meanint       = l.par$imerk.nuc.meanint
l.col$imerk.cell.meanint      = l.par$imerk.cell.meanint

l.col$imnuc.nuc.meanint       = l.par$imnuc.nuc.meanint

l.col$imrec.rec.meanint       = l.par$imrec.rec.meanint

l.col$met.objnum = l.par$met.objnum
l.col$met.fov = l.par$met.fov
l.col$met.frame = l.par$met.frame
l.col$met.trackid = l.par$met.trackid
l.col$pos.x = l.par$pos.x
l.col$pos.y = l.par$pos.y

# col names to be assigned in this script
l.col$joindist = 'join_dist'
l.col$met.trackiduni = 'track_id_uni'

# Load data ----

dt.ts = data.table::fread(file.path(l.par$dir.root, l.par$f.ts))
if (opt$debug)
  cat(sprintf('Reading time series data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.ts)))

dt.rec = data.table::fread(file.path(l.par$dir.root, l.par$f.rec))
if (opt$debug)
  cat(sprintf('Reading receptor data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.rec)))


# Process & clean timeseries data ----
cat(sprintf('\nCleaning data\n'))

# select only columns with usefull measurements (Intensity, AreaShape, Position)

l.col$meas = names(dt.ts)[grep('.*_Intensity_.*|.*_AreaShape_.*|.*_Neighbors_.*|.*_Location_.*', names(dt.ts))]

# 1. Remove all cells at various time points (not tracks!) without cytosol identified

dt.ts = dt.ts[get(l.col$imerk.cyto.ring.meanint) > 0]

# 2. Extract uninterrupted timecourses
## max 1-frame break (can be many breaks)
## Trajectories that last entire experiment

dt.ts.sel = LOCtrajExtr(in.dt = dt.ts, 
                    in.max.break = l.par$ts.maxbreak, 
                    in.met.series = l.col$met.fov, 
                    in.met.t = l.col$met.frame, 
                    in.met.tracklabel = l.col$met.trackid, 
                    in.met.tracklabeluni = l.col$met.trackiduni,
                    in.aggr.cols = l.col$meas)

if (opt$debug)
  cat(sprintf("%d full tracks selected\n", length(unique(dt.ts.sel[[l.col$met.trackiduni]]))))


# Add track ids to receptor data ----
cat(sprintf('\nMerging receptor data with track IDs\n'))

# column names based on which the fuzzy join will operate
v.cols.join = c(l.col$pos.x, 
                l.col$pos.y)

# dt with last frame from time series data; only columns necessary for merge with receptor data
dt.ts.sel.lastframe = dt.ts.sel[get(l.col$met.frame) == max(dt.ts.sel[[l.col$met.frame]]), 
                                c(l.col$met.fov, 
                                  l.col$met.trackid, 
                                  v.cols.join), with = F]
# rescale data; necessary if intensities also used together with position for fuzzy merge
# dt.ts.sel.lastframe[, (v.cols.join) := lapply(.SD, function(x) as.vector(scale(x))), by = c(l.col$met.fov), .SDcols = v.cols.join]

# dt with receptor data and necessary columns
dt.rec.tmp = dt.rec[, c(l.col$met.fov, 
                        l.col$met.objnum,
                        v.cols.join), with = F]

# rescale data; necessary if intensities also used together with position for fuzzy merge
#dt.rec.tmp[, (v.cols.join) := lapply(.SD, function(x) as.vector(scale(x))), by = c(l.col$met.fov), .SDcols = v.cols.join]

# apply fuzzy join individually to every FOV
l.mer = lapply(sort(unique(dt.ts.sel.lastframe[[l.col$met.fov]])), function(x) {
  loc.dt = fuzzyjoin::distance_join(dt.ts.sel.lastframe[get(l.col$met.fov) == x], 
                                    dt.rec.tmp[get(l.col$met.fov) == x, c( l.col$met.objnum,
                                                                           v.cols.join), with = F],
                                    by = v.cols.join,
                                    max_dist = l.par$max.dist,
                                    distance_col = l.col$joindist, 
                                    mode = 'left')
  
  if (opt$debug)
    cat(sprintf("FOV %d, nrows %d\n", x, nrow(loc.dt)))
  
  return(loc.dt)
})

l.mer = Filter(nrow, l.mer)

dt.rec.mer = do.call(rbind, l.mer)

# report the number of unmacthed cells
# failure to match comes from too large of a distance between closest matches (see l.par$max.dist)
# such tracks will be removed from the dataset
dt.rec.mer.fail = dt.rec.mer[is.na(dt.rec.mer[[l.col$joindist]])]
n.failedmatches = nrow(dt.rec.mer.fail)

if (n.failedmatches > 0 ) {
  if (opt$debug)
    cat(sprintf("%d track(s) failed to match with receptor snapshot\n", n.failedmatches))
  
  # remove failed matches from receptor dataset
  dt.rec.mer = dt.rec.mer[!is.na(dt.rec.mer[[l.col$joindist]])]
  
  # remove failed matches from time course data
  setkeyv(dt.rec.mer.fail, c(l.col$met.fov, l.col$met.trackid))
  setkeyv(dt.ts.sel, c(l.col$met.fov, l.col$met.trackid))
  dt.ts.sel = dt.ts.sel[!dt.rec.mer.fail]
  
}

# if two objects are close within the threshold, fuzzy_join will join two objects from y
# here, we take the closest match

dt.rec.mer.dupl = dt.rec.mer[(duplicated(dt.rec.mer[, c(l.col$met.fov, l.col$met.trackid), with = F]))]

if (nrow(dt.rec.mer.dupl) > 0 ) {
  if (opt$debug)
    cat(sprintf('%d duplicate match(es)', nrow(dt.rec.mer.dupl)))

  # From: https://stackoverflow.com/a/55279497/1898713
  # Order by join_dist for every unique FOV & trackID, take first element
  dt.rec.mer = setorderv(dt.rec.mer, c(l.col$met.fov, l.col$met.trackid, l.col$joindist))[, .SD[1L], by = c(l.col$met.fov, l.col$met.trackid)]
    
}


# Save ----  
dir.tmp = file.path(l.par$dir.root, l.par$dir.data)
cat(sprintf('\nSaving data in:\n%s\n\n', dir.tmp))

# Create directory for data output
if(!dir.exists(dir.tmp))
   dir.create(dir.tmp)

# receptor
dt.rec.final = merge(dt.rec[, c(l.col$met.fov, 
                                l.col$met.objnum, 
                                l.col$imrec.rec.meanint ), with = F], 
                     dt.rec.mer[, c(l.col$met.fov, 
                                    l.col$met.objnum, 
                                    l.col$met.trackid), with = F])

fname.tmp = file.path(dir.tmp, l.par$f.out.rec)
if (opt$debug)
  cat(sprintf("Writing receptor data to\n%s\n\n", fname.tmp))

fwrite(x = LOCsignif_dt(dt.rec.final, 6), file = fname.tmp)
gzip(fname.tmp, overwrite = T)

# time series
dt.ts.final = dt.ts.sel[, names(dt.ts.sel)[!(names(dt.ts.sel) %in% l.col$met.trackiduni)], with = F]

fname.tmp = file.path(dir.tmp, l.par$f.out.ts)
if (opt$debug)
  cat(sprintf("Writing time series data to\n%s\n\n", fname.tmp))

fwrite(x = LOCsignif_dt(dt.ts.final, 6), file = fname.tmp)
gzip(fname.tmp, overwrite = T)

cat(sprintf('Analysis finished successfully!\n'))
