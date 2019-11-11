# Process data from Coralie's siRNA screen
# 1. Remove cells without cytosol identified
# 2. Remove dark cells in single-cell time series based on the average of first 5 timepoints of nuclear ERK-KTR.
# 3. Identifies long trajectories from single-cell time-series data from LAP track (merged, 1-line-header csv.gz)
# 4. Remove dark and bright cells (below 5% and above 95% filter (5-95%)) based on receptor level.
# 5. Assigns track-ids to cells in receptor snaphots (via fuzzy merge)
# 6. Plot single-cell time series
# 7. Plot population averages
# 8. Plot densities of ERK-KTR and receptor with bounds for data removal

require(data.table, quietly = T)
require(R.utils, quietly = T)
require(fuzzyjoin, quietly = T)
require(optparse, quietly = T)
require(ggplot2, quietly = T)
require(readxl, quietly = T)

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
  loc.dt[, c(in.met.tracklabeluni) := paste0(sprintf("%03d", get(in.met.series)), sprintf("_%04d", get(in.met.tracklabel)))]
  
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

#' Check for the existence of a list element, add it to another list
#'
#' Given two lists and the name of the element, check for the existance of the element in the first list.
#' If it exists, add it and its value to the second list.
#' If it doesn't exist, add it to the second list and assign the default value.
#'
#' @param inLfrom List to check from.
#' @param inLto List to add en element to.
#' @param inCheckVal String with the name of a list element to check for, and to add the value stored by that element to the new list.
#' @param inDefVal The default value. Assign it to inCheckVal if it doesn't exist the inLfrom list.
#'
#' @return List with an element added
#' @export
#'
#' @examples
#' 
#' list1 = list(a = "text", b = 1)
#' list2 = list()
#' 
#' # The element exists in list1 and will be added to list2
#' LOCcheckAndAddListElem(list1, list2, "a", "default text")
#' 
#' # Element doesn't exist in list1. This element will be added to list2 with a default value
#' LOCcheckAndAddListElem(list1, list2, "c", 0)
#' 
#' # Element doesn't exist in list1. This element will be added to list1 with a default value
#' LOCcheckAndAddListElem(list1, inLto = NULL, "c", 0)
#' 

LOCcheckAndAddListElem = function(inLfrom, inLto = NULL, inCheckVal, inDefVal, inDeb = F) {
  
  locLreturn = NULL
  
  if (equals(inLfrom, inLto) | is.null(inLto)) {
    # if checking in one list only:
    # - if the element exists in inLfrom, do nothing
    # - if the element doesn't exist in inLfrom, add the element to inLfrom and assign the default value
    if (!exists(inCheckVal, where = inLfrom)) {
      inLfrom[[inCheckVal]] = inDefVal
      
      if (inDeb)
        cat(sprintf("Element %s doesn't exist in the inLfrom list. Adding it to inLfrom with a default value: %s\n",
                    inCheckVal, inDefVal))
    }
    locLreturn = inLfrom
  } else {
    # if both lists are different objects, 
    # - if the element exists in inLfrom, transfer the value from inLfrom to inLto, 
    # - if the element doesn't exist in inLfrom, add the element to inLto and assign the default value
    if (exists(inCheckVal, where = inLfrom)) {
      inLto[[inCheckVal]] = inLfrom[[inCheckVal]]    
      
      if (inDeb)
        cat(sprintf("Element %s exists in the inLfrom list. Transfering from inLfrom to inLto\n",
                    inCheckVal))
    } else {
      inLto[[inCheckVal]] = inDefVal
      
      if (inDeb)
        cat(sprintf("Element %s doesn't exist in the inLfrom list. Adding it to inLto with a default value: %s\n",
                    inCheckVal, inDefVal))
    }
    locLreturn = inLto
  }
  
  return(locLreturn)
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
                        help="full path to root directory", metavar="character"),
  optparse::make_option(c("-f", "--plotformat"), type="character", default=NULL, 
                        help="full path to file with parameters of the analysis, csv or xlsx e.g. cp.out/plotFormat.xlsx", metavar="character"),
  optparse::make_option(c("-d", "--opt$debug"), action="store_true", default = FALSE, dest = 'debug', 
                        help="Print extra output")
); 

opt_parser = optparse::OptionParser(option_list=option_list);
opt = optparse::parse_args(opt_parser);

# Temporary for executing from with RStudio
# opt$rootdir = "~/Projects/Olivier/Coralie/20191014_NIH3T3_syst_optoFGFR1_siPOOL_100ms100per"
# opt$plotformat = "~/Projects/Olivier/Coralie/20191014_NIH3T3_syst_optoFGFR1_siPOOL_100ms100per/plotFormat-5q50q.csv"
# opt$debug = TRUE

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

l.par = LOCreadPar(file.path(opt$plotformat))

l.par$dir.root = opt$rootdir
cat(sprintf('Working in:\n%s\n\n', l.par$dir.root))

## Set other variables ----

l.col = list()

# check for parameters present in l.par; add default values, if missing from the plotFormat file
l.col = LOCcheckAndAddListElem(l.par, l.col, "imerk.cyto.ring.meanint", "objCyto_ring_Intensity_MeanIntensity_imErk")
l.col = LOCcheckAndAddListElem(l.par, l.col, "imerk.nuc.meanint", "objNuc_Intensity_MeanIntensity_imErk")
l.col = LOCcheckAndAddListElem(l.par, l.col, "imerk.cell.meanint", "obj_ring_Intensity_MeanIntensity_imErk")
l.col = LOCcheckAndAddListElem(l.par, l.col, "imnuc.nuc.meanint", "objNuc_Intensity_MeanIntensity_imNuc")
l.col = LOCcheckAndAddListElem(l.par, l.col, "imrec.rec.meanint", "obj_Rec_Intensity_MeanIntensity_imRec")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.objnum", "objNuc_ObjectNumber")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.fov", "Image_Metadata_Site")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.treat", "Stimulation_treatment")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.group", "Grouping")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.rt", "RealTime")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.frame", "Image_Metadata_T")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.trackid", "track_id")
l.col = LOCcheckAndAddListElem(l.par, l.col, "pos.x", "objNuc_Location_Center_X")
l.col = LOCcheckAndAddListElem(l.par, l.col, "pos.y", "objNuc_Location_Center_Y")

# col names to be assigned in this script
l.col$joindist = 'join_dist'
l.col$met.trackiduni = 'track_id_uni'
l.col$ratioERK = 'ratioERK'

# checking other parameters; assign defaults if absent from the plotFormt file
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "f.ts",  "objNuc_1line_clean_tracks.csv")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "f.rec", "objNuc.csv")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "f.pm",  "platemap.xlsx")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "ts.maxbreak", 1)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "ts.maxframebase", 5)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "ts.stimt", 9)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "ts.freq", 1)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.stim.h", 0.25)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.yaxis.min", 0.0)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.yaxis.max", 1.5)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.col.summary", "#D65252")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.col.stim", "#4879EF")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.xlab", 'Time (min)')
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.ylab", 'C/N ERK-KTR')
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "rec.int.min", 0.05)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "rec.int.max",  0.5)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "max.dist", 10)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "f.out.ts", "tCoursesSelected_cleaned.csv")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "f.out.rec", "receptorMatchedToTracks_cleaned.csv")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "cytoNucRatio", TRUE)


# initialsie plot list
l.p = list()

# Load data ----

# single-cell timeseries
dt.ts = data.table::fread(file.path(l.par$dir.root, l.par$f.ts))
if (opt$debug)
  cat(sprintf('Reading time series data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.ts)))

# receptor snapshots 
dt.rec = data.table::fread(file.path(l.par$dir.root, l.par$f.rec))
if (opt$debug)
  cat(sprintf('Reading receptor data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.rec)))

# platemap
dt.pm = as.data.table(readxl::read_xlsx(file.path(l.par$dir.root, l.par$f.pm), sheet = 1, skip = 2, col_names = T))
if (opt$debug)
  cat(sprintf('Reading receptor data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.pm)))


# Clean timeseries data ----
cat(sprintf('\nCleaning data\n'))

# 1. Remove all cells at various time points (not tracks!) without cytosol identified
dt.ts = dt.ts[get(l.col$imerk.cyto.ring.meanint) > 0]


# 2. Extract uninterrupted timecourses
# max 1-frame break (can be many breaks)
# Trajectories that last entire experiment

# select only columns with usefull measurements (Intensity, AreaShape, Position)
# used for interpolation?
l.col$meas = names(dt.ts)[grep('.*_Intensity_.*|.*_AreaShape_.*|.*_Neighbors_.*|.*_Location_.*', names(dt.ts))]

dt.ts = LOCtrajExtr(in.dt = dt.ts, 
                        in.max.break = l.par$ts.maxbreak, 
                        in.met.series = l.col$met.fov, 
                        in.met.t = l.col$met.frame, 
                        in.met.tracklabel = l.col$met.trackid, 
                        in.met.tracklabeluni = l.col$met.trackiduni,
                        in.aggr.cols = l.col$meas)

if (opt$debug)
  cat(sprintf("%d full tracks selected\n", length(unique(dt.ts[[l.col$met.trackiduni]]))))

# 3. Remove time-series based on ERK-KTR level
# calculate the mean of first 4 time points
dt.ts.base = dt.ts[get(l.col$met.frame) < l.par$ts.maxframebase, .(xxx = mean(get(l.col$imerk.cell.meanint))), by = c(l.col$met.trackiduni)]

l.p$dens.ts = ggplot(data = dt.ts.base, aes_string(x = paste0('log10(xxx)'))) +
  geom_density() +
  geom_vline(xintercept = log10(l.par$imerk.cell.meanint.thresh), color = 'red', linetype = 2) +
  xlab(paste0('log10 ', l.col$imerk.cell.meanint)) +
  ggtitle('Distribution of ERK-KTR mean fl.int. from the whole\narea bounded by the outer edge of the ring') +
  theme_bw()

#l.p$dens.ts

# keep tracks with baseline above the threshold
setkeyv(dt.ts, l.col$met.trackiduni)
setkeyv(dt.ts.base, l.col$met.trackiduni)
dt.ts = dt.ts[dt.ts.base[xxx > l.par$imerk.cell.meanint.thresh]]
dt.ts[, xxx := NULL]

# clean
rm(dt.ts.base)


# Clean receptor data ----

# Calculate quantiles by taking the region defined in the input parameter file
v.quantiles = quantile(x = dt.rec[[l.col$imrec.rec.meanint]], probs = c(l.par$rec.int.min, l.par$rec.int.max))

# make a plot for QC; save later to a file
l.p$dens.rec = ggplot(data = dt.rec, aes_string(x = paste0('log10(', l.col$imrec.rec.meanint, ')'))) +
  geom_density() +
  geom_vline(xintercept = log10(v.quantiles), colour = 'red', linetype = 2) +
  xlab(paste0('log10 ', l.col$imrec.rec.meanint)) +
  ggtitle('Distribution of Receptor mean fl.int. from\nthe whole-cell area') +
  theme_bw()

# select a percentile region of receptor intensities
dt.rec = dt.rec[get(l.col$imrec.rec.meanint) > v.quantiles[1] & get(l.col$imrec.rec.meanint) < v.quantiles[2]]

# clean
rm(v.quantiles)

# Add track ids to receptor data ----
cat(sprintf('\nMerging receptor data with track IDs\n'))

# column names based on which the fuzzy join will operate
v.cols.join = c(l.col$pos.x, 
                l.col$pos.y)

# dt with last frame from time series data; only columns necessary for merge with receptor data
dt.ts.lastframe = dt.ts[get(l.col$met.frame) == max(dt.ts[[l.col$met.frame]]), 
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
# dt.rec.tmp[, (v.cols.join) := lapply(.SD, function(x) as.vector(scale(x))), by = c(l.col$met.fov), .SDcols = v.cols.join]


# Fuzzy join on individual FOVs ----
# Add this stage, receptor data is reduced but time series data from the last frame is full, 
# hence left-join.

# list of all FOVs available in the data
vFOVs = sort(unique(dt.ts.lastframe[[l.col$met.fov]]))

l.mer = lapply(vFOVs, function(x) {
  # data for a single frame
  loc.dt.ts  = dt.ts.lastframe[get(l.col$met.fov) == x]
  loc.dt.rec = dt.rec.tmp[get(l.col$met.fov) == x, 
                          c( l.col$met.objnum,
                             v.cols.join), with = F]
  
  # checking number of rows
  loc.nrow.dt.ts  = nrow(loc.dt.ts)
  loc.nrow.dt.rec = nrow(loc.dt.rec)
  
  if (opt$debug) {
    cat(sprintf("FOV %d, nrows dt.ts.lastframe = %d, dt.rec.tmp = %d\n", 
                x, loc.nrow.dt.ts, loc.nrow.dt.rec))
  }
    
  # Perform fuzzyjoin only if both tables have measurements
  if (loc.nrow.dt.ts & loc.nrow.dt.rec) {
    loc.dt = fuzzyjoin::distance_join(x = loc.dt.ts, 
                                      y = loc.dt.rec,
                                      method = "euclidean",
                                      by = v.cols.join,
                                      max_dist = l.par$max.dist,
                                      distance_col = l.col$joindist, 
                                      mode = 'left')
    if (opt$debug)
      cat(sprintf("   nrows after fuzzyjoin = %d\n", nrow(loc.dt)))
    
    # If no matches found from y, loc.dt WILL NOT HAVE join_dist column.
    # Columns from y will be filled with NAs.
    # Check for that and add join_dist column filled with NAs
    
    if (!(l.col$joindist %in% names(loc.dt))) {
      loc.dt[, (l.col$joindist) := NA]
      
      if (opt$debug)
        cat("   fuzzyjoin did not find any matching records from y\n")
    }
    
  } else {
    loc.dt = NULL
    if (opt$debug)
      cat("   fuzzyjoin not performed\n")
  }
  
  
  return(loc.dt)
})


# select list elements with nrows > 0
l.mer = Filter(nrow, l.mer)

# convert list to data.table
dt.rec.mer = do.call(rbind, l.mer)

# clean
rm(l.mer, dt.rec.tmp, dt.ts.lastframe, v.cols.join, vFOVs)

# report the number of unmatched cells
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
  setkeyv(dt.ts, c(l.col$met.fov, l.col$met.trackid))
  dt.ts = dt.ts[!dt.rec.mer.fail]
  
  rm(dt.rec.mer.fail)
}

rm(n.failedmatches)


# If two objects are close within the threshold, fuzzy_join will join two objects from y
# here, we take the closest match

dt.rec.mer.dupl = dt.rec.mer[(duplicated(dt.rec.mer[, c(l.col$met.fov, l.col$met.trackid), with = F]))]

if (nrow(dt.rec.mer.dupl) ) {
  if (opt$debug)
    cat(sprintf('%d duplicate match(es)', nrow(dt.rec.mer.dupl)))
  
  # From: https://stackoverflow.com/a/55279497/1898713
  # Order by join_dist for every unique FOV & trackID, take first element
  dt.rec.mer = setorderv(dt.rec.mer, c(l.col$met.fov, l.col$met.trackid, l.col$joindist))[, .SD[1L], by = c(l.col$met.fov, l.col$met.trackid)]
  
  rm(dt.rec.mer.dupl)
}

# create final receptor dt
dt.rec.final = merge(dt.rec[, c(l.col$met.fov, 
                                l.col$met.objnum, 
                                l.col$imrec.rec.meanint ), with = F], 
                     dt.rec.mer[, c(l.col$met.fov, 
                                    l.col$met.objnum, 
                                    l.col$met.trackid), with = F])

rm(dt.rec, dt.rec.mer)


# Save data ----  
dir.tmp = file.path(l.par$dir.root, l.par$dir.data)
cat(sprintf('\nSaving data in:\n%s\n\n', dir.tmp))

# Create directory for data output
if(!dir.exists(dir.tmp))
  dir.create(dir.tmp)


fname.tmp = file.path(dir.tmp, l.par$f.out.rec)
if (opt$debug)
  cat(sprintf("Writing receptor data to\n%s\n\n", fname.tmp))

fwrite(x = LOCsignif_dt(dt.rec.final, 6), file = fname.tmp)
gzip(fname.tmp, overwrite = T)

# time series
dt.ts.final = dt.ts[, names(dt.ts)[!(names(dt.ts) %in% l.col$met.trackiduni)], with = F]

fname.tmp = file.path(dir.tmp, l.par$f.out.ts)
if (opt$debug)
  cat(sprintf("Writing time series data to\n%s\n\n", fname.tmp))

fwrite(x = LOCsignif_dt(dt.ts.final, 6), file = fname.tmp)
gzip(fname.tmp, overwrite = T)
rm(dt.ts.final)

rm(dir.tmp, fname.tmp)


## Plot single-cell trime series ----

# Add metadata to single-cell timeseries

if (l.par$cytoNucRatio) {
  dt.ts[, (l.col$ratioERK) := get(l.col$imerk.cyto.ring.meanint) / get(l.col$imerk.nuc.meanint)]
} else {
  dt.ts[, (l.col$ratioERK) := get(l.col$imerk.nuc.meanint)]
}

dt.ts.plot = dt.ts[, c(l.col$met.fov, l.col$met.frame, l.col$met.trackiduni, l.col$ratioERK), with = F]
dt.ts.plot = merge(dt.ts.plot, 
                   dt.pm[, c('Position', 'Well', 'Stimulation_treatment', 'Grouping'), with = F], 
                   by.x = l.col$met.fov, 
                   by.y = 'Position')

# Add realtime based on acquisition frequency
dt.ts.plot[, (l.col$met.rt) := get(l.col$met.frame) * l.par$ts.freq]
dt.ts.plot[, (l.col$met.frame) := NULL]

# create & save plots per functional group
v.grouping = unique(dt.ts.plot[[l.col$met.group]])
s.dir.plots = file.path(l.par$dir.root, l.par$dir.plots)

cat(sprintf('\nSaving plots in:\n%s\n\n', s.dir.plots))

if(!dir.exists(s.dir.plots))
  dir.create(s.dir.plots)


# temp dt with coordinates of a segment to indicate stimulation
dt.stim = data.table(x = l.par$ts.stimt,
                     xend = l.par$ts.stimt,
                     y = l.par$plot.yaxis.min,
                     yend = l.par$plot.yaxis.min + l.par$plot.stim.h,
                     group = 1)

l.dummy = lapply(v.grouping, function(x) {
  locP = ggplot(data = dt.ts.plot[get(l.col$met.group) == x], 
                aes_string(x = l.col$met.rt, 
                           y = l.col$ratioERK, 
                           group = l.col$met.trackiduni)) +
    geom_line(alpha = 0.1) +
    geom_segment(data = dt.stim, 
                 aes(x = x, xend = xend, 
                     y = y, yend = yend, 
                     group = group), 
                 color = l.par$plot.col.stim) +
    facet_wrap(as.formula(paste0('~', l.col$met.treat, '+', l.col$met.fov))) +
    stat_summary(fun.y=mean, 
                 aes(group=1),
                 geom="line", 
                 lwd=1, 
                 color = l.par$plot.col.summary) +
    xlab(l.par$plot.xlab) +
    ylab(l.par$plot.ylab) +
    coord_cartesian(ylim = c(l.par$plot.yaxis.min, 
                             l.par$plot.yaxis.max)) +
    theme_bw()
  
  locFname = file.path(s.dir.plots, sprintf('tcourses_group%02d.pdf', x))
  ggsave(filename = locFname, plot = locP, width = 8, height = 6)
  
  return(NULL)
})

# Plot population averages ----

dt.ts.plot.aggr = dt.ts.plot[, 
                             .(xxx = mean(get(l.col$ratioERK))), 
                             by = c(l.col$met.group,
                                    l.col$met.treat, 
                                    l.col$met.rt)]
setnames(dt.ts.plot.aggr, 'xxx', l.col$ratioERK)

l.dummy = lapply(v.grouping, function(x) {
  locP = ggplot(data = dt.ts.plot.aggr[get(l.col$met.group) == x], 
                aes_string(x = l.col$met.rt, 
                           y = l.col$ratioERK)) +
    geom_line(aes_string(color = l.col$met.treat)) +
    geom_segment(data = dt.stim, 
                 aes(x = x, xend = xend, 
                     y = y, yend = yend, 
                     group = group), 
                 color = l.par$plot.col.stim) +
    scale_color_discrete('siRNA:') +
    xlab(l.par$plot.xlab) +
    ylab(l.par$plot.ylab) +
    coord_cartesian(ylim = c(l.par$plot.yaxis.min, 
                             l.par$plot.yaxis.max)) +
    theme_bw()
  
  locFname = file.path(s.dir.plots, sprintf('averages_group%02d.pdf', x))
  ggsave(filename = locFname, plot = locP, width = 6, height = 4)
  
  return(NULL)
})


# Plot distributions ----

l.dummy = Map(function (x, i) {
  
  locFname = file.path(s.dir.plots, sprintf('densities_%s.pdf', i))
  ggsave(filename = locFname, plot = x, width = 5, height = 4)
}, l.p, names(l.p))


cat(sprintf('Analysis finished successfully!\n'))

# Save data for single-cell plots ----
s.fname.tmp = file.path(l.par$dir.root, l.par$dir.data, paste0(gsub('.csv', '', l.par$f.out.ts), '_CNerkWithMeta.csv'))
fwrite(x = LOCsignif_dt(dt.ts.plot, 6), file = s.fname.tmp, row.names = F)
gzip(s.fname.tmp, overwrite = T)

