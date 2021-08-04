## README

Code to reproduce most of the figures originally done by Benjamin Lehne in light of the meQTL project.
Also includes code to put together figure 6 done by Matthias Heinig and Johann Hawe.

Data have been gathered manually from ICL/HMGU.

Available scripts:

- plot_chessboard.R: The chessboard plot showing all cosmo associations and the CpG/SNP densities (Fig 1A)
- plot_cis_hist.R: THe histogram of cis association distances (EDF 2A)
- plot_figure6.R: Figure 6 - loads individual parts and puts them together with minor adjustments. The asthma.example.pdf nees to be extracted/saved in a R env with R ver < 4.0.0;
- plot_longrange_boxplots.R: Boxplots of proportions of longrange associations over trans associations.
- plot_manhattan.R: Manhattan plot showing trans associations based on cosmo object. Gene names have been added manually based on the original figure by BL (TODO would be to do this autmatically based on nearest genes from sentinel SNP of specific region).

File chessboard_bp_070617.sh was received from Marie Loh who extracted it
from Benjamin Lehne's original scripts. THe file formed the basis for the 'plot_chessboard.R' script.

