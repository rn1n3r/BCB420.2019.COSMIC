# `BCB420.2019.STRING`

#### (STRING data annotatation of human genes)

&nbsp;

###### [Boris Steipe](https://orcid.org/0000-0002-1134-6758), Department of Biochemistry and Department of Molecular Genetics, University of Toronto, Canada. &lt;boris.steipe@utoronto.ca&gt;

----

**If any of this information is ambiguous, inaccurate, outdated, or incomplete, please check the [most recent version](https://github.com/hyginn/BCB420.2019.STRING) of the package on GitHub and [file an issue](https://github.com/hyginn/BCB420.2019.STRING/issues).**

----

# 1 About this package:

...


----

# 2 STRING Data

STRING is a database of functional interactions. Interactions are inferred from a variety of different experimental and computational categories, scored with a confidence score and made available as network edges where the nodes are genes with an ENSEMBL ID.

This document describes work with [STRING version 11 (preview) (2019-01-11)](https://string-db.org/cgi/access.pl?footer_active_subpage=archive) [(Szclarczyk _et al._ 2019)](https://academic.oup.com/nar/article/47/D1/D607/5198476).



&nbsp;

# 3 Data download and cleanup

1. Navigate to the [**STRING** database](https://string-db.org) and follow the link to the [download section](https://string-db.org/cgi/download.pl).
2. Choose "Homo sapiens" as organism.
3. Download data files. Warning: large.

* `9606.protein.links.v11.0.txt.gz` (71.2 Mb)	protein network data (scored links between proteins);
* `9606.protein.aliases.v11.0.txt.gz` (13.3 Mb)	aliases for STRING proteins: locus names, accessions, descriptions....

4. Uncompress the files and place them into a sister directory of your working directory which is called `data`. (It should be reachable with `file.path("..", "data")`). **Warning:**  `../data/9606.protein.aliases.v11.0.txt` is 172.7 Mb; 
`../data/9606.protein.links.v11.0.txt` is 541 Mb.

&nbsp;

# 4 Mapping ENSEMBL IDs to HGNC symbols

STRING network nodes are Ensembl protein IDs. These can usually be mapped to HGNC symbols, but there might be ambiguities because alternatively spliced proteins might different ENSP IDs that  map to the same HGNC symbol. To provide the best possible interpretation, we need to build a map. This requires care, because it is not guaranteed that all ENSP IDs can be mapped uniquely.** However, the usability of the dataset for annotation depends on the quality of this mapping.**

```R
  # Read the alias data: we need that to find which ENSP IDs map to which
  # HGNC symbols
  tmp <- readr::read_tsv(file.path("../data", "9606.protein.aliases.v11.0.txt"),
                  skip = 1,
                  col_names = c("ENSP", "ID", "source"))  # 2,224,813 rows

  # do they all have the right tax id?
  all(grepl("^9606\\.", tmp$ENSP))  # TRUE
  # remove "9606." prefix
  tmp$ENSP <- gsub("^9606\\.", "", tmp$ENSP)

  # how many IDs do we have to map?
  uENSP <- unique(tmp$ENSP)  # 19,566

  # fetch HGNC reference data
  myURL <- paste0("https://github.com/hyginn/",
                  "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
  load(url(myURL))  # loads HGNC data frame


  # map ENSP to HGNC symbols

  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  ensp2sym <- biomaRt::getBM(filters = "ensembl_peptide_id",
                             attributes = c("ensembl_peptide_id",
                                            "hgnc_symbol"),
                             values = uENSP,
                             mart = myMart)

  colnames(ensp2sym) <- c("ensp", "sym")

  # validate
  any(is.na(ensp2sym$ensp))
  any(is.na(ensp2sym$sym))
  nrow(ensp2sym)                 # 19,169
  length(unique(ensp2sym$ensp))  # 19,166  - three duplicates
  length(unique(ensp2sym$sym))   # 18,892  - 277 duplicates


  # we will not be missing important data.

  # Next we address the duplicated ensp IDs

  dupEnsp <- ensp2sym$ensp[duplicated(ensp2sym$ensp)]
  ensp2sym[ensp2sym$ensp %in% dupEnsp, ]

  #                  ensp      sym
  # 11419 ENSP00000380933  PLEKHG7
  # 11420 ENSP00000380933 C12orf74
  # 13253 ENSP00000480558   CCL3L3
  # 13254 ENSP00000480558   CCL3L1
  # 13331 ENSP00000344961  PLEKHG7
  # 13332 ENSP00000344961 C12orf74

  # ENSP00000380933 and ENSP00000344961 should both map to PLEKHG7
  # CCL3L3 and CCL3L3 both have UniProt ID P16619, we map ENSP00000480558
  # (arbitrarily) to CCL3L1

  # validate rows
  ensp2sym[ensp2sym$sym %in% c("C12orf74", "CCL3L3"), ]

  # remove rows
  ensp2sym <- ensp2sym[ ! (ensp2sym$sym %in% c("C12orf74", "CCL3L3")), ]

  # check result
  any(duplicated(ensp2sym$ensp))   # now FALSE

  # now we can use the ensp IDs as rownames
  rownames(ensp2sym) <- ensp2sym$ensp

  # finish up: add the unmapped ENSPs

  x <- data.frame(ensp = unmappedENSP,
                  sym = NA,
                  stringsAsFactors = FALSE)
  rownames(x) <- x$ensp
  ensp2sym <- rbind(ensp2sym, x)

  # do we now have everything mapped?
  all(uENSP %in% ensp2sym$ensp)  # TRUE

```

## 4.1 Data cleanup

...


&nbsp;

# 5 A script for annotating gene sets

 ...



&nbsp;

# 6 Annotating the example set

 ...



&nbsp;

# 7 Further reading

* Szklarczyk, D., Gable, A. L., Lyon, D., Junge, A., Wyder, S., Huerta-Cepas, J., Simonovic, M., Doncheva, N. T., Morris, J. H., Bork, P., Jensen, L. J., & von Mering, C. (2019). STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. _Nucleic acids research_, D1, D607-D613.

&nbsp;

# 8 Acknowledgements

Thanks to Simon Kågedal's [PubMed to APA reference tool](http://helgo.net/simon/pubmed/).

&nbsp;

&nbsp;

<!-- END -->
