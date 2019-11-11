# Process data from Coralie's siRNA screen
# 1. Read data used for plotting produced in the main script "extractCleanPlotTrajFromSinglePulse.R"
# 2. Plot single-cell time series
# 3. Plot population averages

require(data.table, quietly = T)
require(R.utils, quietly = T)
require(optparse, quietly = T)
require(ggplot2, quietly = T)
require(readxl, quietly = T)

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

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

# Temporary for executing from with RStudio

#opt$rootdir = "~/Projects.tmp/Coralie/20191014_NIH3T3_syst_optoFGFR1_siPOOL_100ms100per"
#opt$plotformat = "~/Projects.tmp/Coralie/20191014_NIH3T3_syst_optoFGFR1_siPOOL_100ms100per/plotFormat-5q50q.csv"
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
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.fov", "Image_Metadata_Site")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.treat", "Stimulation_treatment")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.group", "Grouping")
l.col = LOCcheckAndAddListElem(l.par, l.col, "met.rt", "RealTime")

# col names assigned in the main script
l.col$met.trackiduni = 'track_id_uni'
l.col$ratioERK = 'ratioERK'

# checking other parameters; assign defaults if absent from the plotFormt file
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "ts.stimt", 9)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.stim.h", 0.25)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.yaxis.min", 0.0)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.yaxis.max", 1.5)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.col.summary", "#D65252")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.col.stim", "#4879EF")
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.xlab", 'Time (min)')
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "plot.ylab", 'C/N ERK-KTR')
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "rec.int.min", 0.05)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "rec.int.max",  0.5)
l.par = LOCcheckAndAddListElem(l.par, inLto = NULL, "f.out.ts", "tCoursesSelected_cleaned.csv")


# initialsie plot list
l.p = list()

# Load data ----


s.fname.tmp = file.path(l.par$dir.root, l.par$dir.data, paste0(gsub('.csv', '', l.par$f.out.ts), '_CNerkWithMeta.csv.gz'))
dt.ts.plot = fread(file = s.fname.tmp)

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
  
  locFname = file.path(s.dir.plots, sprintf('tcourses2ndPass_group%02d.pdf', x))
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
  
  locFname = file.path(s.dir.plots, sprintf('averages2ndPass_group%02d.pdf', x))
  ggsave(filename = locFname, plot = locP, width = 6, height = 4)
  
  return(NULL)
})


cat(sprintf('Analysis finished successfully!\n'))
