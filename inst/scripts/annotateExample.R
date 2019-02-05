# annotateExample.R
# Describes the annotation of the example gene set provided using the cosmicData
# imported through BCB420.2019.COSMIC::importCosmic()
# Saves the annotated data to xSetCosmic.tsv in extdata/, and plots graphs

# Import if not found in environment
if (!exists("cosmicData")) {
  BCB420.2019.COSMIC::importCosmic()
}

genes <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")

notInCosmic <- genes[!(genes %in% cosmicData$HGNCsym)]
hgnc$`Approved name`[match(notInCosmic, hgnc$`Approved symbol`)]

geneSetData <- cosmicData[cosmicData$HGNCsym %in% genes,]

# Save annotated set
write.table(geneSetData, file='inst/extdata/xSetCosmic.tsv', quote=FALSE,
            sep='\t', col.names = TRUE, row.names = FALSE, na="")

# Get the number of times each gene symbol appears in geneSetData
counts <- table(geneSetData$HGNCsym)
counts <- sort(counts, decreasing=TRUE)
barplot(counts[1:10], main="Count of HGNC Symbol in Gene Set")

# Most common sites
counts <- table(geneSetData$`Primary site`)
counts <- sort(counts, decreasing=TRUE)
barplot(counts[1:10], main="Count of Primary Site in Gene Set")
