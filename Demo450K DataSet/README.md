
# Demo450K DataSet

@(Diary)[ChAMP]

This 450K data set contains 8 samples, 4 cancer and 4 normal, which is perfect for testing 450K dataset with ChAMP. This data set is lung cancer. We have incoparated this data set into ChAMPdata package for user to use as default Demo data set.

## Code：

``` r
# Set data path for loading, this dataset was contained in ChAMPdata.
testDir=system.file("extdata",package="ChAMPdata")
# Load Data
myLoad <- champ.load(testDir,arraytype="450K")
myLoad_2 <- champ.load(testDir,method="minfi",arraytype="450K")
# Check Distribution of CpGs on chromosome, Island Regions and so on.
CpG.GUI()
# Normlize with BMIQ method.
myNorm <- champ.norm(myLoad$beta) # Tried 4 methods below
# Normlize with PBC method.
myNorm_2 <- champ.norm(myLoad$beta,method="PBC")
# Normlize with SWAN method.
myNorm_3 <- champ.norm(beta=myLoad_2$beta,rgSet=myLoad_2$rgSet,mset=myLoad_2$mset,method="SWAN")
# Normlize with FunctionNormalization method.
myNorm_4 <- champ.norm(myLoad_2$beta,rgSet=myLoad_2$rgSet,method="FunctionalNormalization")
# Check Quality Control of this data set
QC.GUI(myNorm,pheno=myLoad$pd$Sample_Group)
# Do SVD check on data set
champ.SVD()
# Correct batch effect on "Slide" factor.
myCombat <- champ.runCombat(myNorm,pd = myLoad$pd,batchname="Slide")
# Detect DMP
myDMP <- champ.DMP()
# Check DMP result
DMP.GUI()
# Calculate DMR with Bumphunter Method
myDMR <- champ.DMR()
# Calculate DMR with DMRcate Method
myDMR_2 <- champ.DMR(method="DMRcate",cores=1)
# Calculate DMR with ProbeLasso Method
myDMR_3 <- champ.DMR(method="ProbeLasso",minProbes = 3)
# Check DMR result
DMR.GUI()
# Calculate GSEA result with goseq method.
champ.GSEA()
# Calculate Differential Methylation Block with champ.GSEA()
myBlock <- champ.Block()
# Check DMB result.
Block.GUI()
# Calculate EpiMod with FEM package.
champ.EpiMod()
# Do reference base cell correction.
myrefbase <- champ.refbase()
# Calculate Copy Number Variance.
myCNA <- champ.CNA()
```

---

Firstly we use "ChAMP" method to load data. ChAMP method would NOT return rgSet and mset.

### champ.load() result
``` r
> testDir=system.file("extdata",package="ChAMPdata")
> myLoad <- champ.load(testDir,arraytype="450K")
[===========================]
[<<<< ChAMP.LOAD START >>>>>]
-----------------------------

[ Loading Data with ChAMP Method ]
----------------------------------
Note that ChAMP method will NOT return rgSet or mset, they object defined by minfi. Which means, if you use ChAMP method to load data, you can not use SWAN or FunctionNormliazation method in champ.norm() (you can use BMIQ or PBC still). But All other function should not be influenced.

[===========================]
[<<<< ChAMP.IMPORT START >>>>>]
-----------------------------

[ Section 1: Read PD Files Start ]
  CSV Directory: /home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/lung_test_set.csv
  Find CSV Success
  Reading CSV File
  Replace Sentrix_Position into Array
  Replace Sentrix_ID into Slide
[ Section 1: Read PD file Done ]


[ Section 2: Read IDAT files Start ]
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R03C02_Grn.idat ---- (1/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R05C02_Grn.idat ---- (2/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/9247377086_R01C01_Grn.idat ---- (3/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/9247377086_R02C01_Grn.idat ---- (4/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7766130112_R06C01_Grn.idat ---- (5/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7766130112_R01C02_Grn.idat ---- (6/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R01C01_Grn.idat ---- (7/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R01C02_Grn.idat ---- (8/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R03C02_Red.idat ---- (1/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R05C02_Red.idat ---- (2/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/9247377086_R01C01_Red.idat ---- (3/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/9247377086_R02C01_Red.idat ---- (4/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7766130112_R06C01_Red.idat ---- (5/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7766130112_R01C02_Red.idat ---- (6/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R01C01_Red.idat ---- (7/8)
  Loading:/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/7990895118_R01C02_Red.idat ---- (8/8)

  Extract Mean value for Green and Red Channel Success
    Your Red Green Channel contains 622399 probes.
[ Section 2: Read IDAT Files Done ]


[ Section 3: Use Annotation Start ]

  Reading 450K Annotation >>

  Fetching NEGATIVE ControlProbe.
    Totally, there are 613 control probes in Annotation.
    Your data set contains 613 control probes.

  Generating Meth and UnMeth Matrix
    Extracting Meth Matrix...
      Totally there are 485512 Meth probes in 450K Annotation.
      Your data set contains 485512 Meth probes.
    Extracting UnMeth Matrix...
      Totally there are 485512 UnMeth probes in 450K Annotation.
      Your data set contains 485512 UnMeth probes.

  Generating beta Matrix
  Generating M Matrix
  Generating intensity Matrix
  Calculating Detect P value
  Counting Beads
[ Section 3: Use Annotation Done ]

[<<<<< ChAMP.IMPORT END >>>>>>]
[===========================]
[You may want to process champ.filter() next.]

[===========================]
[<<<< ChAMP.FILTER START >>>>>]
-----------------------------

In New version ChAMP, champ.filter() function has been set to do filtering on the result of champ.import(). You can use champ.import() + champ.filter() to do Data Loading, or set "method" parameter in champ.load() as "ChAMP" to get the same effect.

This function is provided for user need to do filtering on some beta (or M) matrix, which contained most filtering system in champ.load except beadcount. User need to input beta matrix, pd file themselves. If you want to do filterintg on detP matrix and Bead Count, you also need to input a detected P matrix and Bead Count information.

Note that if you want to filter more data matrix, say beta, M, intensity... please make sure they have exactly the same rownames and colnames.


[ Section 1:  Check Input Start ]
  You have inputed beta,intensity for Analysis.

  pd file provided, checking if it's in accord with Data Matrix...
    pd file check success.

  Parameter filterDetP is TRUE, checking if detP in accord with Data Matrix...
    detP check success.

  Parameter filterBeads is TRUE, checking if beadcount in accord with Data Matrix...
    beadcount check success.

  parameter autoimpute is TRUE. Checking if the conditions are fulfilled...
    !!! ProbeCutoff is 0, which means you have no needs to do imputation. autoimpute has been reset FALSE.

  Checking Finished :filterDetP,filterBeads,filterMultiHit,filterSNPs,filterNoCG,filterXY would be done on beta,intensity.
  You also provided :detP,beadcount .
[ Section 1: Check Input Done ]


[ Section 2: Filtering Start >>

  Filtering Detect P value Start
    The fraction of failed positions per sample
    You may need to delete samples with high proportion of failed probes:

   Failed CpG Fraction.
C1         0.0013429122
C2         0.0022162171
C3         0.0003563249
C4         0.0002842360
T1         0.0003831007
T2         0.0011946152
T3         0.0014953286
T4         0.0015447610

    Filtering probes with a detection p-value above 0.01.
    Removing 2728 probes.
    If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples

  Filtering BeadCount Start
    Filtering probes with a beadcount <3 in at least 5% of samples.
    Removing 9291 probes

  Filtering NoCG Start
    Only Keep CpGs, removing 2959 probes from the analysis.

  Filtering SNPs Start
    Using general 450K SNP list for filtering.
    Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
    Removing 49231 probes from the analysis.

  Filtering MultiHit Start
    Filtering probes that align to multiple locations as identified in Nordlund et al
    Removing 7003 probes from the analysis.

  Filtering XY Start
    Filtering probes located on X,Y chromosome, removing 9917 probes from the analysis.

  Updating PD file

  Fixing Outliers Start
    Replacing all value smaller/equal to 0 with smallest positive value.
    Replacing all value greater/equal to 1 with largest value below 1..
[ Section 2: Filtering Done ]

 All filterings are Done, now you have 404383 probes and 8 samples.

[<<<<< ChAMP.FILTER END >>>>>>]
[===========================]
[You may want to process champ.QC() next.]

[<<<<< ChAMP.LOAD END >>>>>>]
[===========================]
[You may want to process champ.QC() next.]

>
```

Then we used "minfi" method to load data, which should return the same result as "ChAMP".

``` r
> myLoad_2 <- champ.load(testDir,arraytype="450K",method="minfi")
[===========================]
[<<<< ChAMP.LOAD START >>>>>]
-----------------------------

[ Loading Data with Minfi Method ]
----------------------------------
Loading data from /home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata
[read.metharray.sheet] Found the following CSV files:
[1] "/home/tianyuan/R/x86_64-redhat-linux-gnu-library/3.4/ChAMPdata/extdata/lung_test_set.csv"
Loading required package: IlluminaHumanMethylation450kmanifest
<< Read DataSet Success. >>

The fraction of failed positions per sample

            (You may need to delete samples with high proportion of failed probes
):
   Failed CpG Fraction.
C1         0.0013429122
C2         0.0022162171
C3         0.0003563249
C4         0.0002842360
T1         0.0003831007
T2         0.0011946152
T3         0.0014953286
T4         0.0015447610
Filtering probes with a detection p-value above 0.01 in one or more samples has removed 2728 probes from the analysis. If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples.
<< Filter DetP Done. >>


There is no NA values in your matrix, there is no need to imputation.

Filtering probes with a beadcount <3 in at least 5% of samples, has removed 9291 from the analysis.
<< Filter Beads Done. >>

Filtering non-cg probes, has removed 2959 from the analysis.
<< Filter NoCG Done. >>

Using general 450K SNP list for filtering.
Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed 49231 from the analysis.
<< Filter SNP Done. >>

Filtering probes that align to multiple locations as identified in Nordlund et al, has removed 7003 from the analysis.
<< Filter MultiHit Done. >>

Filtering probes on the X or Y chromosome has removed 9917 from the analysis.
<< Filter XY chromosome Done. >>

[Beta value is selected as output.]

Zeros in your dataset have been replaced with smallest positive value.

One in your dataset have been replaced with largest value below 1.

The analysis will be proceed with 404383 probes and 8 samples.

Current Data Set contains 0 NA in [Beta] Matrix.

[<<<<< ChAMP.LOAD END >>>>>>]
[===========================]
[You may want to process champ.QC() next.]

>
```

### CpG.GUI() Result

Then we can use CpG.GUI() to check the distribution of these CpGs on chromosomes.

``` r
CpG.GUI()
```
![Alt text](./1490873132787.png)


### champ.norm() Result

Then we demonstrate how to use 4 method to do normalization. The firstly one is "BMIQ" method.

``` r
> myNorm <- champ.norm(myLoad$beta)
[===========================]
[>>>>> ChAMP.NORM START <<<<<<]
-----------------------------
champ.norm Results will be saved in ./CHAMP_Normalization/
[ SWAN method call for BOTH rgSet and mset input, FunctionalNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. ]

<< Normalizing data with BMIQ Method >>
Note that,BMIQ function may fail for bad quality samples (Samples did not even show beta distribution).
3 cores will be used to do parallel BMIQ computing.
[>>>>> ChAMP.NORM END <<<<<<]
[===========================]
[You may want to process champ.SVD() next.]

>
```

The second normalization method we test here is "PBC".

``` r
> myNorm_2 <- champ.norm(myLoad$beta,method="PBC")
[===========================]
[>>>>> ChAMP.NORM START <<<<<<]
-----------------------------
champ.norm Results will be saved in ./CHAMP_Normalization/
[ SWAN method call for BOTH rgSet and mset input, FunctionalNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. ]

<< Normalizing data with PBC Method >>
[1] "Done for sample 1"
[1] "Done for sample 2"
[1] "Done for sample 3"
[1] "Done for sample 4"
[1] "Done for sample 5"
[1] "Done for sample 6"
[1] "Done for sample 7"
[1] "Done for sample 8"
[>>>>> ChAMP.NORM END <<<<<<]
[===========================]
[You may want to process champ.SVD() next.]

>
```

The third normalization method we test here is "SWAN". **This method calls for rgSet and mset, so we have to use myLoad_2 from "minfi" method here.**

``` r
> myNorm_3 <- champ.norm(beta=myLoad_2$beta,rgSet=myLoad_2$rgSet,mset=myLoad_2$mset,method="SWAN")
[===========================]
[>>>>> ChAMP.NORM START <<<<<<]
-----------------------------
champ.norm Results will be saved in ./CHAMP_Normalization/
[ SWAN method call for BOTH rgSet and mset input, FunctionalNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. ]

[>>>>> ChAMP.NORM END <<<<<<]
[===========================]
[You may want to process champ.SVD() next.]

>
```

``` r
> myNorm_4 <- champ.norm(myLoad_2$beta,rgSet=myLoad_2$rgSet,method="FunctionalNormalization")
[===========================]
[>>>>> ChAMP.NORM START <<<<<<]
-----------------------------
champ.norm Results will be saved in ./CHAMP_Normalization/
[ SWAN method call for BOTH rgSet and mset input, FunctionalNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. ]

[preprocessFunnorm] Background and dye bias correction with noob
Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
[preprocessNoob] Applying R/G ratio flip to fix dye bias...
[preprocessFunnorm] Mapping to genome
[preprocessFunnorm] Quantile extraction
[preprocessFunnorm] Normalization
[>>>>> ChAMP.NORM END <<<<<<]
[===========================]
[You may want to process champ.SVD() next.]

>
```


### QC.GUI() Result
``` r
QC.GUI(myNorm,pheno=myLoad$pd$Sample_Group)
```
![Alt text](./1490874689557.png)



### champ.SVD() Result
``` r

> champ.SVD(myNorm)
[===========================]
[<<<<< ChAMP.SVD START >>>>>]
-----------------------------
champ.SVD Results will be saved in ./CHAMP_SVDimages/ .

[SVD analysis will be proceed with 404383 probes and 8 samples.]


[ champ.SVD() will only check the dimensions between data and pd, instead if checking if Sample_Names are correctly matched (because some user may have no Sample_Names in their pd file),thus please make sure your pd file is in accord with your data sets (beta) and (rgSet).]

<< Following Factors in your pd(sample_sheet.csv) will be analysised: >>
<Sample_Group>(character):C, T
<Sample_Well>(character):E09, G09, E02, F02, B09, C09, E08
<Slide>(character):7990895118, 9247377086, 7766130112
<Array>(factor):R03C02, R05C02, R01C01, R02C01, R06C01, R01C02
[champ.SVD have automatically select ALL factors contain at least two different values from your pd(sample_sheet.csv), if you don't want to analysis some of them, please remove them manually from your pd variable then retry champ.SVD().]

<< Following Factors in your pd(sample_sheet.csv) will not be analysis: >>
<Sample_Name>
<Sample_Plate>
<Pool_ID>
<Project>
[Factors are ignored because they only indicate Name or Project, or they contain ONLY ONE value across all Samples.]

<< PhenoTypes.lv generated successfully. >>
<< Calculate SVD matrix successfully. >>
<< Plot SVD matrix successfully. >>
[<<<<<< ChAMP.SVD END >>>>>>]
[===========================]
[If the batch effect is not significant, you may want to process champ.DMP() or champ.DMR() or champ.BlockFinder() next, otherwise, you may want to run champ.runCombat() to eliminat batch effect, then rerun champ.SVD() to check corrected result.]

>
```
![Alt text](./1490874805462.png)

### champ.runCombat Result
``` r
> myCombat <- champ.runCombat(myNorm,pd = myLoad$pd,batchname="Slide")
[===========================]
[<< CHAMP.RUNCOMBAT START >>]
-----------------------------
<< Preparing files for ComBat >>
[Combat correction will be proceed with 404383 probes and 8 samples.]

<< Following Factors in your pd(sample_sheet.csv) could be applied to Combat: >>
<Slide>(factor):7990895118, 9247377086, 7766130112
[champ.runCombat have automatically select ALL factors contain at least two different values from your pd(sample_sheet.csv).]

<< Following Factors in your pd(sample_sheet.csv) can not be corrected: >>
<Sample_Name>
<Sample_Plate>
<Sample_Group>
<Pool_ID>
<Project>
<Sample_Well>
<Array>
[Factors are ignored because they are conflict with variablename, or they contain ONLY ONE value across all Samples, or some phenotype contains less than 2 Samples.]
As your assigned in batchname parameter: Slide will be corrected by Combat function.

<< Checking confounded status between Slide and Sample_Group >>
--------------------------
Model for Correction is:
~Sample_Group
<environment: 0x402e6a48>
Combat can adjust for  1  covariate(s) or covariate level(s)
--------------------------
<< Rank Check Complete, you data is good to proceed. >> ^_^

<< Start Correcting Slide >>
~Sample_Group
<environment: 0x402e6a48>
Generate mod success. Started to run ComBat, which is quite slow...
Found 3 batches
Adjusting for 1 covariate(s) or covariate level(s)
Standardizing Data across genes
Fitting L/S model and finding priors
Finding parametric adjustments
Adjusting the Data
champ.runCombat success. Corrected dataset will be returned.
>
```


### champ.DMP() Result
``` r
> myDMP <- champ.DMP()
[===========================]
[<<<<< ChAMP.DMP START >>>>>]
-----------------------------
!!! Important !!! New Modification has been made on champ.DMP():

    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.

    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted "pheno" parameter is "numeric" type.

--------------------------------

[ Section 1:  Check Input Pheno Start ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.
    [The power of statistics analysis on groups contain very few samples may not strong.]
    pheno contains only 2 phenotypes
    compare.group parameter is NULL, two pheno types will be added into Compare List.
    C_to_T compare group : C, T

[ Section 1:  Check Input Pheno Done ]


[ Section 2:  Find Differential Methylated CpGs Start ]

  -----------------------------
  Start to Compare : C, T
  Contrast Matrix
      Contrasts
Levels pT-pC
    pC    -1
    pT     1
  You have found 3879 significant MVPs with a BH adjusted P-value below 0.05.
  Calculate DMP for C and T done.

[ Section 2:  Find Numeric Vector Related CpGs Done ]


[ Section 3:  Match Annotation Start ]


[ Section 3:  Match Annotation Done ]

[<<<<<< ChAMP.DMP END >>>>>>]
[===========================]
[You may want to process DMP.GUI() or champ.GSEA() next.]

>
```

### DMP.GUI() Result
``` r
DMP.GUI()
```
![Alt text](./1490877011533.png)



### champ.DMR() Result

There are three method calculate DMR, Bumphunter, DMRcate and ProbeLasso. Here we firstly show the performance of Bumphunter.

``` r
> myDMR <- champ.DMR()
[===========================]
[<<<<< ChAMP.DMR START >>>>>]
-----------------------------
!!! important !!! We just upgrate champ.DMR() function, since now champ.DMP() could works on multiple phenotypes, but ProbeLasso can only works on one DMP result, so if your pheno parameter contains more than 2 phenotypes, and you want to use ProbeLasso function, you MUST specify compare.group=c("A","B"). Bumphunter and DMRcate should not be influenced.

[ Section 1:  Check Input Pheno Start ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.

[ Section 1:  Check Input Pheno Done ]


[ Section 2:  Run DMR Algorithm Start ]

Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
<< Find DMR with Bumphunter Method >>
3 cores will be used to do parallel Bumphunter computing.
According to your data set, champ.DMR() detected 7372 clusters contains MORE THAN 7 probes within300 maxGap. These clusters will be used to find DMR.

[bumphunterEngine] Parallelizing using 3 workers/cores (backend: doParallelMC, version: 1.0.10).
[bumphunterEngine] Computing coefficients.
[bumphunterEngine] Smoothing coefficients.
Loading required package: rngtools
Loading required package: pkgmaker
Loading required package: registry

Attaching package: 'pkgmaker'

The following object is masked from 'package:S4Vectors':

    new2

The following object is masked from 'package:base':

    isNamespaceLoaded

[bumphunterEngine] Performing 250 bootstraps.
[bumphunterEngine] Computing marginal bootstrap p-values.
[bumphunterEngine] Smoothing bootstrap coefficients.
[bumphunterEngine] cutoff: 1.616
[bumphunterEngine] Finding regions.
[bumphunterEngine] Found 1010 bumps.
[bumphunterEngine] Computing regions for each bootstrap.
[bumphunterEngine] Estimating p-values and FWER.
<< Calculate DMR success. >>
Bumphunter detected 184 DMRs with P value <= 0.05.

[ Section 2:  Run DMR Algorithm Done ]

[<<<<<< ChAMP.DMR END >>>>>>]
[===========================]
[You may want to process DMR.GUI() or champ.GSEA() next.]

>
```

Then we show DMRcate here.

``` r
> myDMR_2 <- champ.DMR(method="DMRcate",cores=1)
[===========================]
[<<<<< ChAMP.DMR START >>>>>]
-----------------------------
!!! important !!! We just upgrate champ.DMR() function, since now champ.DMP() could works on multiple phenotypes, but ProbeLasso can only works on one DMP result, so if your pheno parameter contains more than 2 phenotypes, and you want to use ProbeLasso function, you MUST specify compare.group=c("A","B"). Bumphunter and DMRcate should not be influenced.

[ Section 1:  Check Input Pheno Start ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.

[ Section 1:  Check Input Pheno Done ]


[ Section 2:  Run DMR Algorithm Start ]

1 cores will be used to do parallel DMRcate computing.
<< Find DMR with DMRcate Method >>
Loading required package: IlluminaHumanMethylationEPICanno.ilm10b2.hg19

Attaching package: 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19'

The following objects are masked from 'package:IlluminaHumanMethylation450kanno.ilmn12.hg19':

    Islands.UCSC, Locations, Manifest, Other, SNPs.132CommonSingle, SNPs.135CommonSingle, SNPs.137CommonSingle,
    SNPs.138CommonSingle, SNPs.141CommonSingle, SNPs.142CommonSingle, SNPs.144CommonSingle, SNPs.146CommonSingle,
    SNPs.147CommonSingle, SNPs.Illumina

Your contrast returned 4944 individually significant probes. We recommend the default setting of pcutoff in dmrcate().
Fitting chr1...
Fitting chr10...
Fitting chr11...
Fitting chr12...
Fitting chr13...
Fitting chr14...
Fitting chr15...
Fitting chr16...
Fitting chr17...
Fitting chr18...
Fitting chr19...
Fitting chr2...
Fitting chr20...
Fitting chr21...
Fitting chr22...
Fitting chr3...
Fitting chr4...
Fitting chr5...
Fitting chr6...
Fitting chr7...
Fitting chr8...
Fitting chr9...
Demarcating regions...
Done!
DMRcate detected 260 DMRs with mafcut as= 0.05.

[ Section 2:  Run DMR Algorithm Done ]

[<<<<<< ChAMP.DMR END >>>>>>]
[===========================]
[You may want to process DMR.GUI() or champ.GSEA() next.]

Warning messages:
1: In if (arraytype == "450K") { :
  the condition has length > 1 and only the first element will be used
2: In if (arraytype == "EPIC") { :
  the condition has length > 1 and only the first element will be used
>
```

Finally, we show the performance of ProbeLasso here.

``` r
> myDMR_3 <- champ.DMR(method="ProbeLasso",minProbes = 3)
[===========================]
[<<<<< ChAMP.DMR START >>>>>]
-----------------------------
!!! important !!! We just upgrate champ.DMR() function, since now champ.DMP() could works on multiple phenotypes, but ProbeLasso can only works on one DMP result, so if your pheno parameter contains more than 2 phenotypes, and you want to use ProbeLasso function, you MUST specify compare.group=c("A","B"). Bumphunter and DMRcate should not be influenced.

[ Section 1:  Check Input Pheno Start ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.
  ProbeLasso Method can only be done between two phenotypes. So we need to do more check here...
    Your pheno parameter contains extactly two phenotypes, which is good and compare.group is not needed, champ.DMR() would proceed with your whole data set.

[ Section 1:  Check Input Pheno Done ]


[ Section 2:  Run DMR Algorithm Start ]

champ.DMR Results will be saved in ./CHAMP_ProbeLasso/
<< Find DMR with ProbeLasso Method >>
[===========================]
[<<<<< ChAMP.DMP START >>>>>]
-----------------------------
!!! Important !!! New Modification has been made on champ.DMP():

    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.

    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted "pheno" parameter is "numeric" type.

--------------------------------

[ Section 1:  Check Input Pheno Start ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.
    [The power of statistics analysis on groups contain very few samples may not strong.]
    pheno contains only 2 phenotypes
    compare.group parameter is NULL, two pheno types will be added into Compare List.
    C_to_T compare group : C, T

[ Section 1:  Check Input Pheno Done ]


[ Section 2:  Find Differential Methylated CpGs Start ]

  -----------------------------
  Start to Compare : C, T
  Contrast Matrix
      Contrasts
Levels pT-pC
    pC    -1
    pT     1
  You have found 404383 significant MVPs with a BH adjusted P-value below 1.
  Calculate DMP for C and T done.

[ Section 2:  Find Numeric Vector Related CpGs Done ]


[ Section 3:  Match Annotation Start ]


[ Section 3:  Match Annotation Done ]

[<<<<<< ChAMP.DMP END >>>>>>]
[===========================]
[You may want to process DMP.GUI() or champ.GSEA() next.]

<< Get closestProbe for each Probe >>
<< Get lassoQuantileThreshold for each featureCgi >>
<< Get expend ranges for each probe >>
<< Get DMR from overlapped probes >>
<< Get adjusted P value for DMR >>
<< Get Start-End Ranges for each DMR >>
<< Calculate Methylation Scores for each DMR >>
<< Generate Probe-level Data >>
<< Generate DMR metadata >>
ProbeLasso detected 56 DMRs with P value <= 0.05.

[ Section 2:  Run DMR Algorithm Done ]

[<<<<<< ChAMP.DMR END >>>>>>]
[===========================]
[You may want to process DMR.GUI() or champ.GSEA() next.]

>
```


### DMR.GUI() Result
``` r
> DMR.GUI()
!!! important !!! Since we just upgrated champ.DMP() function, which is now can support multiple phenotypes. Here in DMR.GUI() function, if you want to use "runDMP" parameter, and your pheno contains more than two groups of phenotypes, you MUST specify compare.group parameter as compare.group=c("A","B") to get DMP value between group A and group B.

[ Section 1: Calculate DMP Start  ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.
  Your pheno contains EXACTLY two phenotypes, which is good, compare.group is not needed.
Calculating DMP
[===========================]
[<<<<< ChAMP.DMP START >>>>>]
-----------------------------
!!! Important !!! New Modification has been made on champ.DMP():

    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.

    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted "pheno" parameter is "numeric" type.

--------------------------------

[ Section 1:  Check Input Pheno Start ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.
    [The power of statistics analysis on groups contain very few samples may not strong.]
    pheno contains only 2 phenotypes
    compare.group parameter is NULL, two pheno types will be added into Compare List.
    C_to_T compare group : C, T

[ Section 1:  Check Input Pheno Done ]


[ Section 2:  Find Differential Methylated CpGs Start ]

  -----------------------------
  Start to Compare : C, T
  Contrast Matrix
      Contrasts
Levels pT-pC
    pC    -1
    pT     1
  You have found 404383 significant MVPs with a BH adjusted P-value below 1.
  Calculate DMP for C and T done.

[ Section 2:  Find Numeric Vector Related CpGs Done ]


[ Section 3:  Match Annotation Start ]


[ Section 3:  Match Annotation Done ]

[<<<<<< ChAMP.DMP END >>>>>>]
[===========================]
[You may want to process DMP.GUI() or champ.GSEA() next.]


[ Section 1: Calculate DMP Done  ]


[ Section 2: Mapping DMR to annotation Start  ]

  Generating Annotation File
  Generating Annotation File Success

[ Section 2: Mapping DMR to annotation Done  ]


Listening on http://127.0.0.1:5873
```
![Alt text](./1490880175417.png)
![Alt text](./1490887690292.png)
![Alt text](./1490888021028.png)


### champ.GSEA() Result

Here we have two ways to do GSEA, one is troditional Fisher Exact Test, and another is "gometh" function from missMethyl package.

``` r
> myGSEA <- champ.GSEA()
[===========================]
[<<<< ChAMP.GSEA START >>>>>]
-----------------------------
<< Prepare Gene List Ready  >>
<< Start Do GSEA on each Gene List  >>
<< Do GSEA on Gene list DMP>>
<< Pale Fisher Exact Test will be used to do GSEA >>
 << The category information is downloaded from MsigDB, and only simple Fisher Exact Test will be used to calculate GSEA. This method is suitable if your genes has equalivalent probability to be enriched. If you are using CpGs mapping genes, gometh method is recommended.>>
<< Done for Gene list DMP >>
<< Do GSEA on Gene list DMR>>
<< Pale Fisher Exact Test will be used to do GSEA >>
 << The category information is downloaded from MsigDB, and only simple Fisher Exact Test will be used to calculate GSEA. This method is suitable if your genes has equalivalent probability to be enriched. If you are using CpGs mapping genes, gometh method is recommended.>>
<< Done for Gene list DMR >>
[<<<<< ChAMP.GSEA END >>>>>>]
[===========================]
>
```
Then we show the performance of gometh method

``` r
> myGSEA <- champ.GSEA(method="gometh")
[===========================]
[<<<< ChAMP.GSEA START >>>>>]
-----------------------------
<< Prepare CpG List Ready  >>
  Calculating GSEA with gometh method on DMP CpG list
  Note that gometh method would count numbers of CpGs in each genes and correct this bias.
  Calculating GSEA with gometh method on DMR CpG list
  Note that gometh method would count numbers of CpGs in each genes and correct this bias.
[<<<<< ChAMP.GSEA END >>>>>>]
[===========================]
Warning messages:
1: In alias2SymbolTable(flat$symbol) :
  Multiple symbols ignored for one or more aliases
2: In alias2SymbolTable(flat$symbol) :
  Multiple symbols ignored for one or more aliases
>
```

### champ.Block() Result
``` r
> myBlock <- champ.Block()
[===========================]
[<<<< ChAMP.Block START >>>>]
-----------------------------
<< Load Annotation Successfully >>
<< Get Clusters by cgi-info Successfully >>
<< Calculate Average Beta Value Successfully >>
<< Generate Block Position Successfully >>
<< New Clusters are generated for blocks >>
<< Generate information for New Clusters >>
[bumphunterEngine] Parallelizing using 3 workers/cores (backend: doParallelMC, version: 1.0.10).
[bumphunterEngine] Computing coefficients.
[bumphunterEngine] Smoothing coefficients.
[bumphunterEngine] Performing 500 permutations.
[bumphunterEngine] Computing marginal permutation p-values.
[bumphunterEngine] Smoothing permutation coefficients.
[bumphunterEngine] cutoff: 0.046
[bumphunterEngine] Finding regions.
[bumphunterEngine] Found 2165 bumps.
[bumphunterEngine] Computing regions for each permutation.
[bumphunterEngine] Estimating p-values and FWER.
<< Run Bumphunter Successfully >>
[<<<<< ChAMP.BLOCK END >>>>>]
[===========================]
[You may want to process Block.GUI() next.]

>
```


### Block.GUI() Result
``` r
> Block.GUI()
!!! important !!! Since we just upgrated champ.DMP() function, which is now can support multiple phenotypes. Here in Block.GUI() function, if you want to use "runDMP" parameter, and your pheno contains more than two groups of phenotypes, you MUST specify compare.group parameter as compare.group=c("A","B") to get DMP value between group A and group B.
Generation CpG information
  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.
  Your pheno contains EXACTLY two phenotypes, which is good, compare.group is not needed.
Calculating DMP
[===========================]
[<<<<< ChAMP.DMP START >>>>>]
-----------------------------
!!! Important !!! New Modification has been made on champ.DMP():

    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.

    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted "pheno" parameter is "numeric" type.

--------------------------------

[ Section 1:  Check Input Pheno Start ]

  You pheno is character type.
    Your pheno information contains following groups. >>
    <C>:4 samples.
    <T>:4 samples.
    [The power of statistics analysis on groups contain very few samples may not strong.]
    pheno contains only 2 phenotypes
    compare.group parameter is NULL, two pheno types will be added into Compare List.
    C_to_T compare group : C, T

[ Section 1:  Check Input Pheno Done ]


[ Section 2:  Find Differential Methylated CpGs Start ]

  -----------------------------
  Start to Compare : C, T
  Contrast Matrix
      Contrasts
Levels pT-pC
    pC    -1
    pT     1
  You have found 404383 significant MVPs with a BH adjusted P-value below 1.
  Calculate DMP for C and T done.

[ Section 2:  Find Numeric Vector Related CpGs Done ]


[ Section 3:  Match Annotation Start ]


[ Section 3:  Match Annotation Done ]

[<<<<<< ChAMP.DMP END >>>>>>]
[===========================]
[You may want to process DMP.GUI() or champ.GSEA() next.]


Listening on http://127.0.0.1:5873
```
 ![Alt text](./1490882069062.png)


### champ.EpiMod() Result
``` r
 myEpiMod <- champ.EpiMod()
```
![Alt text](./1490886365556.png)


### champ.refase() Result
``` r
> myRefebase <- champ.refbase()
[===========================]
[<<< ChAMP.REFBASE START >>>]
-----------------------------
<< Load projectWBC function success. >>
Mean value for each estimated Cell Proportion:
        CD8T         CD4T           NK        Bcell         Mono         Gran
5.347103e-17 3.643681e-01 1.019458e-01 2.007267e-01 1.912713e-01 1.284261e-01
CD8T has smallest cell proportion, all other cell proportions will be corrected by linear regression method.
All cell proportion influence except the one with least cell proportion get corrected.

>
```


### champ.CNA() Result
``` r
> myCNA <- champ.CNA()
[===========================]
[<<<<< ChAMP.CNA START >>>>>]
-----------------------------
champ.CNA Results will be saved in ./CHAMP_CNA .

ChaMP.CNA does not provide batch Correct on intensity data now, but you can use champ.runCombat to correct slides batch yourself.
<< Create Control Data >>
<< Combining champ bloodCtl dataset into your intensity dataset as control >>
champ bloodCtl dataset contains only two samples, they will be used as control groups.
<< Calculate mean value difference between each sample to mean control samples >>
<< Generate CHR and MAPINFO information >>
<< Processing Samples >>
Analyzing: C1.qn
Analyzing: C2.qn
Analyzing: C3.qn
Analyzing: C4.qn
Analyzing: T1.qn
Analyzing: T2.qn
Analyzing: T3.qn
Analyzing: T4.qn
<< Processing Groups >>
Analyzing: C1.C.qn
Analyzing: C2.C.qn
Analyzing: C3.C.qn
Analyzing: C4.C.qn
Analyzing: T1.T.qn
Analyzing: T2.T.qn
Analyzing: T3.T.qn
Analyzing: T4.T.qn
[<<<<<< ChAMP.CNA END >>>>>>]
[===========================]
>
```
![Alt text](./1491367196104.png)
