The scripts process data from siRNA screen.

The script `extractCleanPlotTrajFromSinglePulse.R` processes raw data:

1. Single-cell time-series from `Merged/cp.out/output.mer/objNuc_1line_clean_tracks.csv.gz`
2. Receptor snapshots from the end of the experiment from `Receptor/cp.out/output.mer/objNuc.csv.gz`
3. Plate map from `Receptor/cp.out/output.mer/objNuc.csv.gz`


The processing involves the following steps:

1. Remove cells without cytosol identified
2. Remove dark cells in single-cell time series based on the average of first 5 timepoints of nuclear ERK-KTR.
3. Identifies long trajectories from single-cell time-series data from LAP track (merged, 1-line-header csv.gz) 
4. Remove dark and bright cells (below 5% and above 95% filter (5-95%)) based on receptor level.
5. Assigns track-ids to cells in receptor snaphots (via fuzzy merge)
6. Plot single-cell time series
7. Plot population averages
8. Plot densities of ERK-KTR and receptor with bounds for data removal

Saves data in `output-data`:

- `receptorMatchedToTracks_cleaned.csv.gz` with receptor data matched with nuclei assigned track ID from time-series data
- `tCoursesSelected_cleaned.csv.gz` with time-series data after all cleaning steps; no metadata
- `tCoursesSelected_cleaned_CNerkWithMeta.csv.gz` with reduced time-series dataset only with C/N ERK=KTR and metadata from platemap

Saves plots in `output-plots`:

- population averages per treatment; one plot per group
- single-cell time-series faceted per treatment; one plot per group
- probability densities for ERK-KTR and receptor from the area encompassed by the outer edge of the ring

Requires following R packages:

- `data.table`
- `R.utils`
- `fuzzyjoin`
- `optparse`
- `ggplot2`
- `readxl`