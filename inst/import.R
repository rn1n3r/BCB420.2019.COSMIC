if (! requireNamespace("readr")) {
  install.packages("readr")
}

if (! requireNamespace("data.table")) {
  install.packages("data.table")
}

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
library(readr)

col_spec <- cols(
  `Gene name` = col_character(),
  `Accession Number` = col_character(),
  `Gene CDS length` = col_integer(),
  `HGNC ID` = col_character(), # "/" character?
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
  Age = col_double() # Change this
)



# col_spec <- cols(
#   `Gene name` = col_character(),
#   `Accession Number` = col_character(),
#   `Gene CDS length` = col_integer(),
#   `HGNC ID` = col_character(), # "/" character?
#   `Sample name` = col_character(),
#   ID_sample = '_', 
#   ID_tumour = '_',
#   `Primary site` = col_character(),
#   `Site subtype 1` = '_',
#   `Site subtype 2` = '_',
#   `Site subtype 3` = '_',
#   `Primary histology` = col_character(),
#   `Histology subtype 1` = '_',
#   `Histology subtype 2` = '_',
#   `Histology subtype 3` = '_',
#   `Genome-wide screen` = '_',
#   `Mutation ID` = '_',
#   `Mutation CDS` = col_character(),
#   `Mutation AA` = col_character(),
#   `Mutation Description` = col_character(),
#   `Mutation zygosity` = col_character(),
#   LOH = '_',
#   GRCh = '_',
#   `Mutation genome position` = '_',
#   `Mutation strand` = '_',
#   SNP = '_',
#   `Resistance Mutation` = '_',
#   `FATHMM prediction` = '_',
#   `FATHMM score` = '_',
#   `Mutation somatic status` = '_',
#   Pubmed_PMID = '_',
#   ID_STUDY = '_',
#   `Sample Type` = '_',
#   `Tumour origin` = '_',
#   Age = '_' # Change this
# )

print("Reading tsv files")
cosmicData <- readr::read_tsv(file.path("../data", "CosmicMutantExport.tsv"), col_names = TRUE, col_types = col_spec)
cosmicData2 <- cosmicData
hgnc <- readr::read_tsv(file.path("../data", "HGNC_data_all.tsv"), col_names=TRUE)
#sum(is.na(tmp$`HGNC ID`))/length(tmp$`HGNC ID`)

uHGNCID <- unique(cosmicData$`HGNC ID`)

# Take the first ID if a / is found
uHGNCID <- sapply(strsplit(uHGNCID, "/"), head, 1)
uHGNCID <- unique(cosmicData$`HGNC ID`)
uHGNCIDPrefix <- paste("HGNC:", uHGNCID, sep="")

#sum(!(uHGNCID) %in% hgnc_table$hgnc_id)

# Create id2sym that maps HGNC ID to HGNC Symbol
# Use data.table in order to utilize key search to speed up computation
index <- match(uHGNCIDPrefix, hgnc$`HGNC ID`)
id2sym <- data.table::data.table(`HGNC ID` = uHGNCID, `Gene name` = hgnc$`Approved symbol`[index])
data.table::setkey(id2sym, `HGNC ID`)

print("Mapping HGNC IDs to symbols")
cosmicData$HGNCsym <- id2sym[cosmicData$`HGNC ID`, allow.cartesian=TRUE]$`Gene name`

# Now, need to deal with the entries with HGNC ID = NA
# Check if the provided gene name is found in the HGNC symbol list
print("Using provided gene name as symbol (validate if found in HGNC list)")
missingNames <- cosmicData$`Gene name`[is.na(cosmicData$HGNCsym)]

# let's also get rid of the _ENSTs
missingNames <- sapply(strsplit(missingNames, "_ENST"), head, 1)
hgncIndex <- match(missingNames, hgnc$`Approved symbol`)

cosmicData$HGNCsym[is.na(cosmicData$HGNCsym)] <- hgnc$`Approved symbol`[hgncIndex]


# Still missing some (perhaps synonyms and outdated symbols)
missingNames <- missingNames[is.na(hgncIndex)] 

# Let's see if we can find any using the Ensembl gene ID
print("Using biomaRt to look for remaining missing entries by using the Ensembl gene ID")
ensgIndex <- grepl("ENSG", missingNames)
ensembl <- missingNames[ensgIndex]
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensemblMap <- biomaRt::getBM(filters="ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = unique(ensembl), mart=myMart)

ensemblMap <- data.table::as.data.table(ensemblMap)
data.table::setkey(ensemblMap, "ensembl_gene_id")
ensemblMap <- ensemblMap[hgnc_symbol != ""]

cosmicData$HGNCsym[is.na(cosmicData$HGNCsym)][ensgIndex] <- ensemblMap[ensembl, allow.cartesian=TRUE]$hgnc_symbol

# couldn't get all of them, let's move on though
missingNames <- cosmicData$`Gene name`[is.na(cosmicData$HGNCsym)]



# print("Searching through previous names and synonyms")
# foundNames = vector(mode="integer", length=length(missingNames))
# i <- 0
# for (name in missingNames) {
#   i <- i+1
#   print(i)
#   index <- grep(paste(name, "($|,)", sep=""), hgnc$`Previous symbols`)
#   if (length(index) >= 1) {
#     foundNames[which(name == missingNames)] <- index[1] # Just take the first match for now
#   }
#   else {
#     index <- grep(paste(name, "($|,)", sep=""), hgnc$`Synonyms`)
#     if (length(index) == 1) {
#       foundNames[which(name == missingNames)] <- index
#     }
#     
#   }
#   
# }
# 
# cosmicData$HGNCsym[is.na(cosmicData$HGNCsym)][foundNames != 0] <- hgnc$`Approved symbol`[foundNames]
# missingNames <- missingNames[foundNames==0]
# 
