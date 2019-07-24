# Process data from Coralie's siRNA screen:
# 1. load single-cell trajectories and receptor data, 
# 2. remove data below 5% and above 95% filter (5-95%) based on receptor level and the avergae of first 5 timepoints of nuclear ERK-KTR
# 3. plot single-cell time series
#
# Relies on data produced by trajFromSinglePulse.R

require(data.table, quietly = T)
require(R.utils, quietly = T)
require(optparse, quietly = T)
require(readxl, quietly = T)
require(ggplot2)


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

# keep inner percentage of data; defined by min,max % bounds
LOCkeepMid = function(inDT, inCol, inLR = c(0.05, 0.95), inType = 3) {
  # find position of the left boundary
  locQuantLR = quantile(inDT[[inCol]], inLR, type = inType)
  
  # select only points inside LR boundaries
  return(inDT[get(l.col$imrec.rec.meanint) > locQuantLR[1] & get(l.col$imrec.rec.meanint) < locQuantLR[2]])
}

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

# lower and upper threshold for data trimming (percentage/100)
l.par$trim.low = 0.05
l.par$trim.high = 0.95

# maximum number of frames to take for averaging. Used for cleaning data based on ERK
l.par$max.frame = 5

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
l.col$met.trackiduni = 'track_id_uni'
l.col$ratioERK = 'ratioERK'

# initialsie plot list
l.p = list()


# Load data ----
# single-cell trajectories from trajFromSinglePusle.R script
dt.ts = data.table::fread(file.path(l.par$dir.root, l.par$dir.data, paste0(l.par$f.out.ts, '.gz')))
if (opt$debug)
  cat(sprintf('Reading time series data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.ts)))

# receptor snapshot from trajFromSinglePusle.R script
dt.rec = data.table::fread(file.path(l.par$dir.root, l.par$dir.data, paste0(l.par$f.out.rec, '.gz')))
if (opt$debug)
  cat(sprintf('Reading receptor data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.rec)))

# platemap
dt.pm = as.data.table(readxl::read_xlsx(file.path(l.par$dir.root, l.par$f.pm), sheet = 1, skip = 2, col_names = T))
if (opt$debug)
  cat(sprintf('Reading receptor data from:\n%s\n\n', file.path(l.par$dir.root, l.par$f.pm)))


# Process & clean timeseries data ----
cat(sprintf('\nCleaning data\n'))

# add unique track ids to receptor data
dt.rec[, (l.col$met.trackiduni) := paste0(get(l.col$met.fov), '_', get(l.col$met.trackid))]

## Trim time-series data based on ERK-KTR level
# calculate the mean of first 4 time points
dt.ts.base = dt.ts[get(l.col$met.frame) < l.par$max.frame, .(xxx = mean(get(l.col$imerk.cell.meanint))), by = c(l.col$met.trackiduni)]

l.p$dens.ts = ggplot(data = dt.ts.base, aes_string(x = paste0('log10(xxx)'))) +
  geom_density() +
  geom_vline(xintercept = log10(l.par$imerk.cell.meanint.thresh), color = 'red', linetype = 2) +
  xlab(paste0('log10 ', l.col$imerk.cell.meanint)) +
  ggtitle('Distribution of ERK-KTR mean fl.int. from the whole\narea bounded by the outer edge of the ring') +
  theme_bw()

#l.p$dens.ts


# vector with cell ids to remove due to low levels of ERK-KTR
v.cellID.tmp1 = dt.ts.base[xxx < l.par$imerk.cell.meanint.thresh][[l.col$met.trackiduni]]


## trim receptor data
# based on receptor level, taking inner region from 5-95%

v.quantiles = quantile(x = dt.rec[[l.col$imrec.rec.meanint]], probs = c(l.par$trim.low, l.par$trim.high))

l.p$dens.rec = ggplot(data = dt.rec, aes_string(x = paste0('log10(', l.col$imrec.rec.meanint, ')'))) +
  geom_density() +
  geom_vline(xintercept = log10(v.quantiles), colour = 'red', linetype = 2) +
  xlab(paste0('log10 ', l.col$imrec.rec.meanint)) +
  ggtitle('Distribution of Receptor mean fl.int. from\nthe whole-cell area') +
  theme_bw()

#l.p$dens.rec

# vector with cell ids to remove due to low levels of ERK-KTR
v.cellID.tmp2 = dt.rec[(get(l.col$imrec.rec.meanint) < v.quantiles[1]) | (get(l.col$imrec.rec.meanint) > v.quantiles[2])][[l.col$met.trackiduni]]

# keep only those time series that are NOT in two lists created above
dt.ts.clean = dt.ts[!(get(l.col$met.trackiduni) %in% unique(c(v.cellID.tmp1, v.cellID.tmp2)))]

## Plot single-cell trime series ----
# Add metadata

dt.ts.clean[, (l.col$ratioERK) := get(l.col$imerk.cyto.ring.meanint) / get(l.col$imerk.nuc.meanint)]

dt.ts.plot = dt.ts.clean[, c(l.col$met.fov, l.col$met.frame, l.col$met.trackiduni, l.col$ratioERK), with = F]
dt.ts.plot = merge(dt.ts.plot, 
                   dt.pm[, c('Position', 'Well', 'Stimulation_treatment', 'Grouping'), with = F], 
                   by.x = l.col$met.fov, 
                   by.y = 'Position')


# create & save plots per functional group
v.grouping = unique(dt.ts.plot[['Grouping']])
s.dir.plots = file.path(l.par$dir.root, 'output-plots')

if(!dir.exists(s.dir.plots))
  dir.create(s.dir.plots)

lapply(v.grouping, function(x) {
    locP = ggplot(data = dt.ts.plot[Grouping == x], aes_string(x = l.col$met.frame, y = l.col$ratioERK, group = l.col$met.trackiduni)) +
      geom_line(alpha = 0.1) +
      facet_wrap(~ Stimulation_treatment + Image_Metadata_Site) +
      stat_summary(fun.y=mean, geom="line", lwd=1, color = 'red', aes(group=1)) +
      xlab('Time (min)') +
      ylab('C/N ERK-KTR') +
      coord_cartesian(ylim = c(0, 1.5))
    
    locFname = file.path(s.dir.plots, sprintf('tcourses_group%02d.pdf', x))
    ggsave(filename = locFname, plot = locP, width = 8, height = 6)
})

# Plot distributions ----

Map(function (x, i) {
  
  locFname = file.path(s.dir.plots, sprintf('densities_%s.pdf', i))
  ggsave(filename = locFname, plot = x, width = 5, height = 4)
}, l.p, names(l.p))

# Save data ----
s.fname.tmp = file.path(l.par$dir.root, l.par$dir.data, 'tCoursesSelected_CNerkWithMeta.csv')
fwrite(x = LOCsignif_dt(dt.ts.plot, 6), file = s.fname.tmp, row.names = F)
gzip(s.fname.tmp, overwrite = T)

