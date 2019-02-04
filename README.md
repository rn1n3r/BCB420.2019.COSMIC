# `BCB420.2019.COSMIC`

#### (COSMIC data annotation of human genes)

###### Edward Ho &lt;emc.ho@mail.utoronto.ca&gt;

----

****

----

## 1 About this package:

This package describes the pipeline to download cancer somatic mutation from the [COSMIC database](https://cancer.sanger.ac.uk/cosmic), how to ensure the mutation entries have a proper HGNC symbol (where applicable), and calculation of basic statistics. This package will also provide an annotation of the example gene set (provided [here](https://github.com/hyginn/BCB420-2019-resources/blob/master/exampleGeneSet.md)), which is based on the phagosome/lysosome fusion system as reviewed by [Corona & Jackson (2018)](https://www.sciencedirect.com/science/article/pii/S0962892418301223).


#### In this project ...

```text
 --BCB420.2019.COSMIC/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.COSMIC.Rproj
   |__DESCRIPTION
   |__dev/
      |__rptTwee.R
      |__toBrowser.R               # display .md files in your browser
   |__inst/
      |__extdata/
         |__ensp2sym.RData         # ENSP ID to HGNC symbol mapping tool
         |__xSetEdges.tsv          # annotated example edges
      |__img/
         |__[...]                  # image sources for .md document
      |__scripts/
         |__recoverIDs.R           # utility to use biomaRt for ID mapping
   |__LICENSE
   |__NAMESPACE
   |__R/
      |__zzz.R
   |__README.md                    # this file

```

&nbsp;

----

## 2 COSMIC Data

COSMIC (Catalogue Of Somatic Mutations in Cancer) is an expert-curated collection of data describing somatic gene mutations found in human cancer. It presents phenotypic data (primary site and histology) associated with the gene mutations that is collected from the literature by experts from the literature, and also utilizes data from the Cancer Gene Census.

All COSMIC data is free for academic users, and requires an email address from an academic institution to register and download data. More information can be found [here](https://cancer.sanger.ac.uk/cosmic/license+&cd=1&hl=en&ct=clnk&gl=ca).


&nbsp;https://academic.oup.com/nar/article/47/D1/D941/5146192

#### 2.1 Data semantics

Association evidence in STRING is compiled in seven "channels" (Szclarczyk _et al._ 2019):

1. **Genomic context I**: neighbourhood
2. **Genomic context II**: fusion
3. **Genomic context III**: phylogenetic profiles
4. **Co-expression**: correlations across a large number of mRNA and proteome data sets
5. **Text-mining**: statistical co-citation analysis across PubMed and OMIM.
6. **Experiments**: these are the _classical_ protein-protein interactions. These scores are derived from data imported from all [IMEX consortium databases](https://www.imexconsortium.org/)
7. **Curated pathway- and protein-complex databases**: pathways from KEGG, Reactome, BioCyc and the GO consortium. All associations derived from this channel are scored with _p_ = 0.9.

For each "channel", STRING distinguishes evidence that is provided from the organism itself, and annotation transfer by homology.

Finally, a composite score is computed as the sum of the probabilities from each channel, subtracting a "prior" (that indicates the probability of a false positive) and applying a homology correction to co-occurrence and text-mining scores. Source code for this operation is available [here (python)](http://string-gamma.org/download/combine_subscores.py).

In summary: STRING interactions represent the probability of "information flow" between two nodes, since: "biologically meaningful interfaces have evolved to allow the flow of information through the cell, and they are ultimately essential for implementing a functional system."(Szclarczyk _et al._ 2019)

This approach integrates biological data on the largest possible scale, and this is successful: a recent benchmark study (Huang _et al._ 2018) has shown STRING to have the highest network recovery scores of disease-associated gene sets from the DisGeNET database, when compared with 20 other gene network resources. 

&nbsp;

## 3 Data download and cleanup

To download the source data from COSMIC ... :

1. Register on the [**COSMIC** website](
    https://cancer.sanger.ac.uk/cosmic/register), log in and navigate to to the [download section](https://cancer.sanger.ac.uk/cosmic/download).
2. Download the COSMIC Mutation Data (`CosmicMutantExport.tsv.gz`).
3. Uncompress the file and place it in <code>../data</code> (sister directory of working directory). <code>CosmicMutantExport.tsv</code> is **1.7 Gb**.

Additionally, we require data from the HGNC (HUGO Gene Nomenclature Committee) to map the HGNC IDs provided by COSMIC to the approved symbol, and also to validate that gene names are valid symbols.

1. Navigate to https://www.genenames.org/download/custom/
2. Keep the default settings and click Submit.
3. Save the text file as <code>hgnc.tsv</code> in the same <code>../data</code> directory as above.

&nbsp;

## 4 Mapping to HGNC symbols

COSMIC data includes the gene name for each observed mutation sample, and this name should be the associated HGNC symbol in most cases. Additionally, the tabulated data has a column for the HGNC ID, allowing the exact symbol to be found in the HGNC database. To ensure that the symbols are accurate, a map of the HGNC ID to the associated will be constructed for the unique HGNC IDs in the COSMIC database.

Some entries have multiple HGNC associated with them. For now, the HGNC symbol that corresponds with the entry's gene name will be used.

Furthermore, there are some cases where there is no HGNC ID, so the provided gene name will be used after validating that it is a valid HGNC symbol. Lastly, there are some gene names that are provided in other formats, such as Ensembl gene ID. These will be mapped (where possible) to the HGNC symbol using the biomaRt package.

&nbsp;

#### Preparations: packages, functions, files

**Required packages:**

**`readr`** to quickly read in the .tsv files provided by COSMIC.

```R
if (! requireNamespace("readr")) {
  install.packages("readr")
}
```



**`data.table`** to make subsetting tables by key [fast](https://cran.r-project.org/web/packages/data.table/vignettes/datatable). This is important as the COSMIC data is quite large (6,581,004 rows).

```R
if (! requireNamespace("data.table")) {
  install.packages("data.table")
}
```



**`biomaRt`** biomaRt is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```



&nbsp;

#### 4.1 Loading the data into R

```R
# Load the HGNC data
hgnc <- readr::read_tsv(file.path("../data", "HGNC_data_all.tsv"), col_names=TRUE)

# Load the COSMIC data with defined columned types:
col_spec <- cols(
  `Gene name` = col_character(),
  `Accession Number` = col_character(),
  `Gene CDS length` = col_integer(),
  `HGNC ID` = col_character(), # This requires col_character because the "/" is used
  `Sample name` = col_character(),
  ID_sample = col_integer(),
  ID_tumour = col_integer(),
  `Primary site` = col_character(),
  `Site subtype 1` = col_character(),
  `Site subtype 2` = col_character(),
  `Site subtype 3` = col_character(),
  `Primary histology` = col_character(),
  `Histology subtype 1` = col_character(),
  `Histology subtype 2` = col_character(),
  `Histology subtype 3` = col_character(),
  `Genome-wide screen` = col_character(),
  `Mutation ID` = col_character(),
  `Mutation CDS` = col_character(),
  `Mutation AA` = col_character(),
  `Mutation Description` = col_character(),
  `Mutation zygosity` = col_character(),
  LOH = col_character(),
  GRCh = col_integer(),
  `Mutation genome position` = col_character(),
  `Mutation strand` = col_character(),
  SNP = col_character(),
  `Resistance Mutation` = col_character(),
  `FATHMM prediction` = col_character(),
  `FATHMM score` = col_double(),
  `Mutation somatic status` = col_character(),
  Pubmed_PMID = col_integer(),
  ID_STUDY = col_integer(),
  `Sample Type` = col_character(),
  `Tumour origin` = col_character(),
  Age = col_double() # Needed to use col_double as certain entries had decimal values
)

# Load data
cosmicData <- readr::read_tsv(file.path("../data", "CosmicMutantExport.tsv"), col_names = TRUE, col_types = col_spec)

# Data preview
head(cosmicData)
# A tibble: 6 x 36
#  `Gene name` `Accession Numb~ `Gene CDS lengt~ `HGNC ID` `Sample name` ID_sample ID_tumour `Primary site`
#  <chr>       <chr>                       <int> <chr>     <chr>             <int>     <int> <chr>         
#1 GRK6        ENST00000355472              1731 4545      PD1403a          898155    815769 stomach       
#2 TMEM108     ENST00000321871              1728 28451     TCGA-13-0760~   1474815   #1398514 ovary         
#3 HRH1        ENST00000397056              1464 5182      TCGA-13-0900~   1474861   #1398560 ovary         
#4 SLC26A8     ENST00000355574              2913 14468     TCGA-13-1509~   1474839   #1398538 ovary         
#5 OR51E1      ENST00000396952               957 15194     GC1_T           1645677   #1560711 stomach       
#6 CMYA5       ENST00000238522             11259 NA        TCGA-AG-3892~   1651564   #1566351 large_intesti~
# ... with 28 more variables: `Site subtype 1` <chr>, `Site subtype 2` <chr>, `Site subtype 3` <chr>,
```

&nbsp;

#### 4.2  Mapping the HGNC ID to HGNC symbol

Many of the entries in the COSMIC database already have an HGNC ID associated with it. In this case, we use this value and look up the corresponding symbol in the imported HGNC data.

###### 4.2.1 Cleaning the ID entries


As noted earlier, sometimes the mutation corresponds to two genes. In this case, we will use the first ID, which usually corresponds with the given gene name. Additionally, we add the HGNC prefix ("HGNC:") to the IDs so that they match the imported data.

```R
uHGNCID <- unique(cosmicData$`HGNC ID`)

# Take the first ID if a / is found
uHGNCID <- sapply(strsplit(uHGNCID, "/"), head, 1)
uHGNCID <- unique(cosmicData$`HGNC ID`)
uHGNCIDPrefix <- paste("HGNC:", uHGNCID, sep="") 

```


###### 4.2.2 Mapping the ID to the symbol

Using the unique HGNC IDs found in the COSMIC data, we create a data.table (<code>id2sym</code>) with the ID and the approved symbol using the <code>hgnc</code> table imported earlier. The data.table will use the HGNC ID as a key, and this allows for fast lookup through [binary search](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-keys-fast-subset.html). Then, using this data.table, a new column <code>HGNCsym</code> is added to the <code>cosmicData</code> table which stores the approved symbol. 

```R
# Create id2sym that maps HGNC ID to HGNC Symbol
# Use data.table in order to utilize key search to speed up computation
index <- match(uHGNCIDPrefix, hgnc$`HGNC ID`)
id2sym <- data.table::data.table(`HGNC ID` = uHGNCID, `Gene name` = hgnc$`Approved symbol`[index])
data.table::setkey(id2sym, `HGNC ID`)

# Use id2sym to map the ID to the symbol
# allow.cartesian=TRUE to allow for duplicate IDs
cosmicData$HGNCsym <- id2sym[cosmicData$`HGNC ID`, allow.cartesian=TRUE]$`Gene name`

# See what our current coverage is for symbol mapping
sum(!is.na(cosmicData$HGNCsym)) / length(cosmicData$HGNCsym) # 0.8840686
```

So currently 88% of the dataset has an associated HGNC symbol. Not bad, but let's see if we can improve that. A quick observation was that some entries which had <code>NA</code> for its HGNC ID had gene names that resembled HGNC approved symbols. 

```R
head(cosmicData$`Gene name`[is.na(cosmicData$`HGNC ID`)])
#[1] "CMYA5"                    "CMYA5"                    "RNFT1_ENST00000305783"   
#[4] "C12orf41_ENST00000420613" "TFE3_ENST00000336239"     "TIAM1_ENST00000286827"   
# CMYA5 is a HGNC symbol, as is TFE3
```



#### 4.3 Using the gene name as the symbol where appropriate

To do this, first we collect the gene names of the entries with missing symbols (<code>missingNames</code>). From above, we can see that some gene names are HGNC symbols, and some of them also have an Ensemble transcript ID attached to it. We first strip the transcript IDs, and then check to see if the name matches the symbols in the HGNC data we imported. Then, we use the indices of the matches to add onto the <code>HGNCsym</code> column of the COSMIC data.

```R

# Now, need to deal with the entries with HGNC ID = NA
# Check if the provided gene name is found in the HGNC symbol list
missingNames <- cosmicData$`Gene name`[is.na(cosmicData$HGNCsym)]

# Let's also get rid of the _ENSTs
missingNames <- sapply(strsplit(missingNames, "_ENST"), head, 1)
hgncIndex <- match(missingNames, hgnc$`Approved symbol`)

# Add matches to the symbol column
cosmicData$HGNCsym[is.na(cosmicData$HGNCsym)] <- hgnc$`Approved symbol`[hgncIndex]

# Check updated coverage
sum(!is.na(cosmicData$HGNCsym)) / length(cosmicData$HGNCsym) # 0.9821093
sum(is.na(cosmicData$HGNCsym)) # 117739

# What do the missing genes look like?
head(cosmicData$`Gene name`[is.na(cosmicData$HGNCsym)])

#[1] "C12orf41_ENST00000420613" "AC093393.1"               "PCNXL3_ENST00000355703"  
#[4] "ENSG00000167390"          "LOC652153"                "LRRC16A"     
```
&nbsp;

#### 4.4 Using biomaRt to map Ensembl gene IDs to HGNC symbols

&nbsp;

```R
# Update names with missing symbols (could not directly use name as symbol)
missingNames <- missingNames[is.na(hgncIndex)]

# Let's see if we can find any using the Ensembl gene ID
ensgIndex <- grepl("ENSG", missingNames)
ensembl <- missingNames[ensgIndex]

# Make a biomaRt and filter for the Ensemble gene ID
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensemblMap <- biomaRt::getBM(filters="ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = unique(ensembl), mart=myMart)

# Use a data.table for fast lookup using keys
ensemblMap <- data.table::as.data.table(ensemblMap)
data.table::setkey(ensemblMap, "ensembl_gene_id")

# Remove empty matches
ensemblMap <- ensemblMap[hgnc_symbol != ""]

# Verify that the returned symbols are HGNC approved symbols
sum(!(ensemblMap$hgnc_symbol %in% hgnc$`Approved symbol`)) # 0


```



#### 4.5 Final validation

Validation and statistics of our mapping tool:

```R



```

&nbsp;


