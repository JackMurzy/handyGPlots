#' Look up gene information from Ensembl
#'
#' Queries the Ensembl REST API to return gene metadata and canonical transcript
#' exon structure. Accepts a gene symbol, Ensembl gene ID (ENSG), or Ensembl
#' transcript ID (ENST).
#'
#' @param input_id A character string. Either a gene symbol (e.g. `"POMC"`),
#'   an Ensembl gene ID (`"ENSG00000115138"`), or an Ensembl transcript ID
#'   (`"ENST00000395826"`).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{SYMBOL}{HGNC gene symbol}
#'     \item{ENSG}{Ensembl gene ID}
#'     \item{ENST}{Canonical transcript ID}
#'     \item{ENSP}{Ensembl protein ID for the canonical transcript}
#'     \item{CHROM}{Chromosome}
#'     \item{STRAND}{Strand (1 or -1)}
#'     \item{GENE_START}{Gene start position (GRCh38)}
#'     \item{GENE_END}{Gene end position (GRCh38)}
#'     \item{EXONS}{data.frame of exon coordinates for the canonical transcript}
#'     \item{CODING_INFO}{Translation start/end for the canonical transcript}
#'   }
#'
#' @examples
#' \dontrun{
#'   ensembl_gene_lookup("POMC")
#'   ensembl_gene_lookup("ENSG00000115138")
#'   ensembl_gene_lookup("ENST00000395826")
#' }
#'
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
#' @export
ensembl_gene_lookup <- function(input_id){
  # ENST -> ENSG via Ensembl REST
  if(grepl("^ENST",input_id)){
    input_id <- sub("\\..*","",input_id)
    url <- paste0("https://rest.ensembl.org/lookup/id/",input_id,"?content-type=application/json&expand=0")
    response <- GET(url)
    if (status_code(response)!=200) stop("Failed to map ENST to ENSG via Ensembl")
    tx_data <- fromJSON(httr::content(response,"text",encoding="UTF-8"))
    input_id <- tx_data$Parent
    if (is.null(input_id) || length(input_id)==0) stop("No ENST-ENSG mapping found")
  }
  # Build endpoint (ENSG or symbol)
  endpoint <- if(grepl("^ENSG",input_id)){
    paste0("/lookup/id/",input_id,"?expand=1")
  } else {
    paste0("/lookup/symbol/human/",input_id,"?expand=1")
  }
  url <- paste0("https://rest.ensembl.org",endpoint,"&content-type=application/json")
  response <- GET(url)
  if(status_code(response)!=200) stop("Failed to retrieve data from Ensembl API")
  data <- fromJSON(httr::content(response,"text",encoding="UTF-8"))
  transcripts <- data$Transcript
  ensp_id <- enst_id <- NA; exons <- data.frame(NULL); translation_info <- NULL;
  if(!is.null(transcripts)){
    canonical_tx <- transcripts[transcripts$is_canonical==1, ]
    enst_id <- canonical_tx$id[1]
    exons <- canonical_tx$Exon[[1]]
    translation_info <- canonical_tx$Translation
    ensp_id <- translation_info$id[1]  # ENSP from canonical translation
  }
  return(list(
    SYMBOL = data$display_name,
    ENSG = data$id,
    ENST = enst_id,
    ENSP = ensp_id,
    CHROM = data$seq_region_name,
    STRAND = data$strand,
    GENE_START = data$start,
    GENE_END = data$end,
    EXONS = exons,
    CODING_INFO = translation_info
  ))
}


#' Look up between UniProt accession and gene name
#'
#' Searches the UniProt REST API and resolves either a gene symbol to its
#' UniProt accession, or a UniProt accession to its gene symbol, depending
#' on the input provided.
#'
#' @param input_id Character. Either an HGNC gene symbol (e.g. `"POMC"`) or a
#'   UniProt accession (e.g. `"P01189"`).
#'
#' @return A character string. If `input_id` is a gene symbol, returns the
#'   best-matching UniProt accession. If `input_id` is a UniProt accession,
#'   returns the corresponding gene symbol.
#'
#' @examples
#' \dontrun{
#'   uniprot_lookup("POMC")    # returns "P01189"
#'   uniprot_lookup("P01189")  # returns "POMC"
#' }
#'
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
#' @importFrom data.table data.table
#' @export
uniprot_lookup <- function(input_id){
  url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=",input_id,"&format=json")
  response <- GET(url)
  if(status_code(response)!=200) stop("Failed to retrieve data from UniProt. Check gene name")
  data <- fromJSON(httr::content(response,"text",encoding="UTF-8"),simplifyDataFrame=FALSE,flatten=FALSE)
  results <- data$results
  accessions <- NULL
  for(i in seq_along(results)){
    if(results[[i]]$organism$taxonId==9606){
      accessions <- rbind(accessions, data.table(
        accession = results[[i]]$primaryAccession,
        gene = results[[i]]$genes[[1]]$geneName$value,
        status = results[[i]]$entryType
      ))
    }
  }
  if(is.null(accessions)) stop("No human results found for the provided gene name")
  is_uniprot <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",input_id)
  if(is_uniprot){
    # Input was an accession so return the gene name
    accession_match <- accessions[accession == input_id]
    if(nrow(accession_match)!=0) return(accession_match$gene[1])
    return(accessions$gene[1])
  } else {
    # Input was a gene name so return the accession (existing logic)
    accession_gene <- accessions[gene==input_id]
    accession_status <- accessions[grepl("UniProtKB reviewed",status)]
    overlap <- merge(accession_status,accession_gene,by=c("accession","status","gene"))
    if(nrow(overlap)!= 0) return(overlap$accession)
    if(nrow(accession_gene)!=0) return(accession_gene$accession)
    if(nrow(accession_status)!=0) return(accession_status$accession)
    return(accessions$accession[1])
  }
}



#' Map protein feature coordinates to genomic positions
#'
#' Uses the Ensembl REST API to convert amino acid positions on a given
#' protein (ENSP) to genomic coordinates.
#'
#' @param ensp_id Character. Ensembl protein ID, e.g. `"ENSP00000360512"`.
#' @param aa_start Integer. Start amino acid position.
#' @param aa_end Integer. End amino acid position.
#'
#' @return A `data.table` with columns `genomic_start`, `genomic_end`,
#'   `strand`, `aa_start`, `aa_end`, or `NULL` if mapping fails.
#'
#' @examples
#' \dontrun{
#'   map_protein_feature("ENSP00000379170", 1, 10)
#' }
#'
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
#' @importFrom data.table data.table
#' @export
map_protein_feature <- function(ensp_id, aa_start, aa_end){
  url <- paste0("https://rest.ensembl.org/map/translation/",
                ensp_id,"/", aa_start,"..",aa_end,
                "?content-type=application/json")
  resp <- GET(url)
  if(status_code(resp)!=200) return(NULL)
  res <- fromJSON(httr::content(resp,"text",encoding="UTF-8"))$mappings
  if(length(res)==0) return(NULL)
  data.table(
    genomic_start = res$start,
    genomic_end = res$end,
    strand = res$strand,
    aa_start = aa_start,
    aa_end = aa_end
  )
}


#' Fetch protein features from UniProt
#'
#' Retrieves annotated features (domains, active sites, etc.) for a
#' protein. The input can be a gene symbol, Ensembl gene/transcript ID, or a
#' UniProt accession - the function resolves all IDs automatically.
#'
#' @param input_id Character. One of: a gene symbol (`"POMC"`), Ensembl gene ID
#'   (`"ENSG00000115138"`), Ensembl transcript ID (`"ENST00000395826"`), or a
#'   UniProt accession (`"P01189"`).
#' @param map_to_genomic Logical. If `TRUE`, amino acid feature coordinates are
#'   mapped to GRCh38 genomic positions via the Ensembl REST API using the ENSP
#'   ID from the canonical transcript. Default `FALSE`.
#' @param type_limit Character vector or NULL. Restrict to specific feature
#'   types (e.g. `c("Domain", "Region")`). Default `NULL`.
#' @param desc_regex Character or NULL. Regex to filter features by description.
#'   Default `NULL`.
#' @param full_data Logical. If `TRUE`, returns a list with both the features
#'   table and the full parsed UniProt JSON. Default `FALSE`.
#'
#' @return A `data.table` of features with columns `type`, `description`,
#'   `start`, `end` (and `genomic_start`, `genomic_end`, `strand` if
#'   `map_to_genomic = TRUE`). If `full_data = TRUE`, a named list with
#'   elements `Features` and `Protein_data`.
#'
#' @examples
#' \dontrun{
#'   # All inputs below resolve to the same protein:
#'   get_protein_features("POMC")
#'   get_protein_features("ENSG00000115138")
#'   get_protein_features("ENST00000395826")
#'   get_protein_features("P01189")
#'
#'   # With genomic coordinate mapping:
#'   get_protein_features("POMC", map_to_genomic = TRUE)
#'
#'   # Filter to peptides only:
#'   get_protein_features("POMC", type_limit = "Peptide", map_to_genomic = TRUE)
#'   
#'   # Filter to helical regions only
#'   get_protein_features("COL1A1", desc_regex = "helical", map_to_genomic = TRUE)
#'   
#'   # Get full data
#'   get_protein_features("POMC",full_data=T)
#' }
#'
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
#' @importFrom data.table as.data.table setorder fifelse rbindlist
#' @export
get_protein_features <- function(input_id, map_to_genomic = FALSE,
                                 type_limit = NULL, desc_regex = NULL,
                                 full_data = FALSE) {
  # Resolve accession or ENSP
  ensp_id <- accession <- NULL
  is_uniprot <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",input_id)
  if(is_uniprot){
    accession <- input_id
    if(map_to_genomic){
      # Resolve gene symbol from accession, then look up ENSP via Ensembl
      gene_symbol <- uniprot_lookup(input_id)
      gene_info <- ensembl_gene_lookup(gene_symbol)
      ensp_id <- gene_info$ENSP
      if(is.null(ensp_id) || is.na(ensp_id)) stop("Could not resolve ENSP ID for genomic mapping")
    }
  } else {
    gene_info <- ensembl_gene_lookup(input_id)
    accession <- uniprot_lookup(gene_info$SYMBOL)
    if(map_to_genomic){
      ensp_id <- gene_info$ENSP
      if(is.null(ensp_id) || is.na(ensp_id)) stop("Could not resolve ENSP ID for genomic mapping")
    }
  }
  
  # Fetch from UniProt
  url <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".json")
  response <- GET(url)
  if(status_code(response)!=200) stop("Failed to retrieve data from UniProt. Check your accession ID")
  protein_data <- fromJSON(httr::content(response,"text",encoding="UTF-8"),simplifyDataFrame=TRUE,flatten=TRUE)
  features <- tryCatch({
    ft <- as.data.table(protein_data$features)[, .(
      type, description,
      start = location.start.value,
      end = location.end.value
    )]
    ft[, description := fifelse(description=="", NA_character_, gsub(",", "", description))]
    setorder(ft, start, end)
  }, error = function(e) { warning("Error extracting features: ", e$message); NULL })
  if(is.null(features)) stop("No features found for accession: ", accession)
  
  # Optional filtering
  if(!is.null(type_limit)){
    features <- features[type %chin% type_limit]
    if(nrow(features)==0){ warning("type_limit retains no rows"); return(NULL)}
  }
  if(!is.null(desc_regex)){
    features <- features[grepl(desc_regex, description, ignore.case = TRUE)]
    if(nrow(features)==0){ warning("desc_regex retains no rows"); return(NULL)}
  }
  
  # Optional mapping
  if(!is.null(ensp_id)){
    mapping_list <- lapply(seq_len(nrow(features)),function(i){
      cat("\nMapping feature", i, "/", nrow(features), "\n")
      row <- features[i, ]
      mapping <- map_protein_feature(ensp_id,row$start,row$end)
      Sys.sleep(0.5)
      if(!is.null(mapping)){
        cbind(row[rep(1,nrow(mapping)), ,drop=FALSE],mapping)
      } else {
        cbind(row,genomic_start=NA,genomic_end=NA,strand=NA)
      }
    })
    features <- rbindlist(mapping_list,fill=TRUE)[,c("aa_start","aa_end"):=NULL]
  }
  if(full_data) return(list(Features = features, Protein_data = protein_data))
  return(features[])
}
