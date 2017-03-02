Mantel test
-----------

Mantel test, correlation between entries of two square matrices containing 'distances' reltive to pairs of individuals, has been vastly used since its introduction in 1967. As stated in the paper 'Dismantling the Mantel tests - Gilles Guillot and Francois Rousset', one of the main problems of this test is its incapability of dealing with matrices that contain some structure and/or autocorrelation. The structure of the data leads to an accrection of type I error, and so far, nobody has shown any good solution for such problem in the genomic world.

Mantel test without correction for population structure
-------------------------------------------------------

In this present analysis, we have decided to apply the Mantel test without any correction for population structure. Three different data sets were analyzed: \* Pairwise correlation of all core genes GRM \* Pairwise correlation of symbiotic genes \* Pairwise correlation of symbiotic genes versus core genes

### Pairwise correlation of symbiotic genes

The results of the pairwise combination of symbiotic genes is shown in the figure below: ![](https://github.com/bvilhjal/nchain/blob/master/mantel_recA_rpoB.pdf)

From this plot, we are missing nodD, nodC, nodB, nodI, nodA and fixN symbiotic genes:. Bjarni has done some filters with SNP data and it seems that those genes did not pass the criteria he's defined. I will try to include them again into our data set.

Mantel test correcting for population structure
-----------------------------------------------

Mantel test, used as a LD approach, relies on the existence of patterns produced exclusively due to genetic variation. However, population structure, local selection or other factors could be confounders in Mantel tests correlations, leading us to wrong conclusions. One way of dealing with this problem is to correct the data matrix for population structure.
