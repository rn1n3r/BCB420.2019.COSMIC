# STRINGwork.R
#
# Purpose: Workflow for downloading and annotating STRING data
# Version: 0.1
# Date:    2018-01-22
# Author:  Boris Steipe <boris.steipe@utoronto.ca>
#          ORCID: 0000-0002-1134-6758
# License: see file LICENSE
#
# ToDo:
# Notes:
#
# ==============================================================================

# WARNING: SIDE EFFECTS
# Executing this script will execute code it contains.

# ====  PACKAGES  ==============================================================

if (! requireNamespace("readr")) {
  install.packages("readr")
}

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}


# ====  FUNCTIONS  =============================================================

recoverID <- function(ensp) {
  # Purpose:
  #     Try to recover an ensp to sym mapping from biomart
  # Parameters:
  #     a: ensp an ensemble peptide ID
  # Value:
  #     result: an HGNC symbol or NA

  # code ...
  sym <- NA
  x <- biomaRt::getBM(filters = "ensembl_peptide_id",
                      attributes = c("uniprotswissprot",
                                     "ucsc"),
                      values = ensp,
                      mart = myMart)
  for (iRow in seq_along(nrow(x))) {
    s <- x[iRow, "uniprotswissprot"]
    s <- HGNC$sym[which(HGNC$UniProtId == s)]
    if (length(s) == 1) {
      sym <- s
    } else {
      s <- x[iRow, "ucsc"]
      s <- HGNC$sym[which(HGNC$UCSCID == s)]
      if (length(s) == 1) {
        sym <- s
      }
    }
  }
  print(ensp)
  if (length(sym) == 1 && is.na(sym) && sym == "") {
    sym <- NA
  }

  return(sym)
}










# ====  PROCESS  ===============================================================
if (FALSE) {

  # STRING is a database of functional interactions. Interactions are scored
  # and made available as network edges.

  # Source data was downloaded from STRING database via organism specific
  # download.
  #
  # https://string-db.org/
  #
  #   ../data/9606.protein.links.v11.0.txt (541) Mb  - contains relationships
  #   ../data/9606.protein.aliases.v11.0.txt (172.7) - contains IDs
  #
  #

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

  # what are these duplicated symbols?
  head(sort(table(ensp2sym$sym), decreasing = TRUE))
  # 227 empty strings among a few "normal" duplicated ID's
  # set the empties to NA
  sel <- ensp2sym$sym == ""
  sum(sel)  # 227
  ensp2sym$sym[sel] <- NA

  #
  # can we map the missing ensp IDs to other cross-references?
  unmappedENSP <- uENSP[ ! (uENSP %in% ensp2sym$ensp)]
  unmappedENSP <- c(unmappedENSP, ensp2sym$ensp[is.na(ensp2sym$sym)])

  x <- unmappedENSP[1:10]
  for (id in x) {
    sym <- recoverID(id)
    print(sprintf("%s %s", sym, x))
  }


  x <- unique(ensp2x$ucsc[ensp2x$ucsc %in%
                                        HGNC$UCSCID])

  sym <- character()
  for (id in x) {
    sym <- c(sym, HGNC$sym[which(id == HGNC$UCSCID)])
  }



  # None found ... nothing we can do. What is this anyway?
  unmappedENSP[1]  # Google -> genecards: Caution: Product of a dubious
  #                               CDS prediction
  unmappedENSP[2]  # Google -> HGNC -> UniProt: Putative uncharacterized
                   #                            protein C20orf78

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

  # final validation
  x <- which(is.na(ensp2sym$sym))[1]
  test <- rbind(head(ensp2sym), ensp2sym[(x - 5):x, ])

  x <- biomaRt::getBM(filters = "ensembl_peptide_id",
                      attributes = c("ensembl_peptide_id",
                                     "hgnc_symbol"),
                      values = test,
                      mart = myMart)

  ensp2sym[x-5,2]




  attributes[grep("ucsc", attributes$description, ignore.case = TRUE),]

  # How many HGNC symbols of the proteins do we find in the STRING
  # ID column?
  unique(HGNC$type)
  # subset
  HGNC <- HGNC[grepl("(protein)|(immunoglobulin)|(TCR)", HGNC$type), ]
  # 19,655 rows

  # how many of these are in STRING?
  x <- HGNC$sym %in% tmp$ID
  sum(x)  # 18,839
  sum(x) * 100 / nrow(HGNC) # 95.8 %

  # Are those symbols that we _can_ map unique?
  tmp2 <- tmp[tmp$ID %in% HGNC$sym, ] # 20,742 is > 19,655: some are not unique

  # let's see some examples:

  x <- tmp2$ID[head(which(duplicated(tmp2$ID)))]
  # [1] "AARSD1" "ACAA2"  "ACAT1"  "ACAT2"  "ACE"    "ACP1"
  tmp2$ENSP[tmp2$ID == x[1]]
  # [1] "ENSP00000400870" "ENSP00000409924"


  head(tmp[x, ])


  sel <- grep("BioMart_HUGO", tmp$source)  # 19,149 rows
  tmp[head(sel), ]
  symMapData <- tmp[sel, ]

  # how many unique?
  x <- unique(symMapData$ID)



  tmp <- tmp[grep("BioMart_HUGO", tmp$source), ]  # 19,119 rows

  tmp <- tmp[tmp$ID %in% Chr20GeneData$sym, ] # 497 of 521 symbols mapped
  ENS2symMap <- tmp$ID                        # extract symbols...
  names(ENS2symMap) <- tmp$ENSP               # ... and use ENSP IDs as names


  # Read the interaction graph data: this is a weighted graph defined as an
  # edge list with gene a, gene b, confidence score (0, 999).
  tmp <- read_delim(paste0(STRINGDIR, "9606.protein.links.v10.5.txt"),
                    delim = " ",
                    skip = 1,
                    col_names = c("a", "b", "score"))  # 11,353,056 rows

  tmp <- tmp[tmp$score >= 900, ]  # 547,620 rows of high-confidence edges

  # Extract edges where both genes are in the mapped Chr 20 genes
  Chr20FuncIntx <- tmp[(tmp$a %in% names(ENS2symMap)) &
                         (tmp$b %in% names(ENS2symMap)), ]  # 564 rows

  # Use ENS2symMap to translate ENSP IDs to HGNC symbols
  Chr20FuncIntx$a <- ENS2symMap[Chr20funcIntx$a]
  Chr20FuncIntx$b <- ENS2symMap[Chr20funcIntx$b]

  # We treat this as an undirected graph, thus we remove duplicates.

  # Add a sorted key column

  Chr20FuncIntx$keys <- character(nrow(Chr20FuncIntx))
  for (i in 1:nrow(Chr20FuncIntx)) {
    Chr20FuncIntx$keys[i] <- paste(sort(c(Chr20FuncIntx$a[i],
                                          Chr20FuncIntx$b[i])),
                                   collapse = ":")
  }

  # remove rows with duplicated keys
  Chr20FuncIntx <- Chr20FuncIntx[! duplicated(Chr20FuncIntx$keys), ]

  # remove unneeded column
  Chr20FuncIntx <- Chr20FuncIntx[ , c("a", "b", "score")]

  # Done
  write_tsv(Chr20FuncIntx, path = "Chr20FuncIntx.tsv")



}

# ====  TESTS  =================================================================
if (FALSE) {
  # Enter your function tests here...

}


# [END]
