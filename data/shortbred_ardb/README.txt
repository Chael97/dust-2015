2015-11-30 RH

Three files copied from: /Users/rhickey/Dropbox (BioBE)/BioBE Team Folder/Gerlinger/metagenome analysis/shortbred/ardb figure

- ARDB_Data_comparison.xlsx: data from meta-analysis by Tiffany Hsu (Huttenhower lab); this includes 3 datasets (HMP subjects, time point 1; Yatsunenko et al.; Qin et al.). Email message from Tiffany Hsu on 11/11/15:
	“Figure S3 is a combination of Jim's Shortbred paper data, a PLoS One paper (Yooseph et al. 2013), and the MBTA data.
	Jim combined three different datasets: a subset of the HMP dataset (subjects, timepoint 1 only), data from Yatsunenko et al (for other countries), and data from Qin et al (for China subjects). He ran shortbred-identify on the antibiotic resistance database (ARDB) using Uniref50 as the reference to generate a marker database. He than ran shortbred-quantify to get a table of subjects (as columns), and antibiotic resistance markers (as rows).
	I took his marker database and ran shortbred-quantify on the shotgun data from the PLoS One paper against it (without filtering or anything, just downloaded the sequences from the BioProject). I also ran shortbred-quantify on the cleaned up MBTA data (Kneaddata, filtering). 
	I then joined all the tables together (they had the same # of rows because they were run against the same marker database), and then threw the entire thing into RStudio and used ggplot2 to generate the plot.  
The processed table is under Supplemental Tables, Table S8, Tab 3, ARDB_DataComparisons. A key for the columns can be found here. The marker database used is here.”

- ARDB_Data_comparison.txt: txt version of xlsx (header removed; units are RPKM)

- categories_metadata.txt: source categories for all samples

- gerlinger_shortbred_merged_ardb.txt: antibiotic resistance database (ARDB) ShortBred results from Gerlinger samples (performed by EH)