# GaBSe
Galaxy Bundle Searcher. Looks for galaxy clusers using coumts in cells.
Designed to use .fits data from DES SVA1 GOLD data release

Files and purposes:
CLUSTERS REPORT.pdf - The report written about the results of this project.

colourGraphing.py - plots colour magnitude diagrams of data exported from datasetToNumpyConverter.py

Counts in Cells.py - The main set of code that runs counts in cells on the DES data. Capable of running on a few regions of the sky, modifiable in the first 40 lines of code to run on any given region of sky at any given cell size. Exports graphs to the graphs folder, graph file names, titles and axes labels are not changed and may not be appropriate for what is being shown.

datasetToNumpyConverter.py - converts .fits file that DES released into 3 numpy files which are used by nearly every other python script here.

GaBSe candidates.npy-  The candidates we used for the project. Import file with numpy, contains a [4,60] array of format [:,n]=[position in RA, position in DEC, Significance, Matched to redMaPPer cluster (1 if it is matched 0 if not)].

README.md - this file.

reddata.npy - the redmapper data containing positions of galaxy clusters, used in Counts in Cells.py
Redmapper data.py - Opens .dat data from redMaPPer and converts it to .npy with relevant positional data

redshift.py - Calculates the expected angular size of galaxy clusters dependent on redshift. Needs no other files to run

Sample generator.py - Generated samples that we used to test the algorithm on in the early stages of the project.
