#' Split an exon into CDS and UTR regions
#'
#' Classifies an exon or splits it into multiple segments based on its
#' overlap with the coding sequence (CDS) boundaries. Returns one row per
#' region segment with a label of either `"CDS"` or `"UTR"`.
#'
#' @param start Integer. Genomic start position of the exon.
#' @param end Integer. Genomic end position of the exon.
#' @param coding_start Integer. Genomic start position of the CDS.
#' @param coding_end Integer. Genomic end position of the CDS.
#'
#' @return A `data.table` with columns `region` (`"CDS"` or `"UTR"`),
#'   `start`, and `end`. May contain 1–3 rows depending on how the exon
#'   overlaps the coding region.
#'
#' @noRd
#' 
split_exon <- function(start, end, coding_start, coding_end) {
  if (end < coding_start | start > coding_end) {
    # Fully outside coding region
    return(data.table(region = "UTR", start = start, end = end))
  } else if (start >= coding_start & end <= coding_end) {
    # Fully inside coding region
    return(data.table(region = "CDS", start = start, end = end))
  } else if (start < coding_start & end > coding_end) {
    # Exon completely spans coding region so split into three
    return(rbind(
      data.table(region = "UTR", start = start, end = coding_start - 1),
      data.table(region = "CDS", start = coding_start, end = coding_end),
      data.table(region = "UTR", start = coding_end + 1, end = end)
    ))
  } else if (start < coding_start & end <= coding_end) {
    # Overlaps coding start
    return(rbind(
      data.table(region = "UTR", start = start, end = coding_start - 1),
      data.table(region = "CDS", start = coding_start, end = end)
    ))
  } else if (start >= coding_start & start <= coding_end & end > coding_end) {
    # Overlaps coding end
    return(rbind(
      data.table(region = "CDS", start = start, end = coding_end),
      data.table(region = "UTR", start = coding_end + 1, end = end)
    ))
  }
}


#' Axis transformation for compressing genomic intervals
#'
#' Creates a custom `scales` transformation that squashes specified genomic
#' intervals (long introns) to a shorter display width, while keeping
#' all other coordinates proportional. Intended for use with
#' `ggplot2::scale_x_continuous(trans = ...)`.
#'
#' @param from Numeric vector. Start positions of the intervals to compress.
#' @param to Numeric vector. End positions of the intervals to compress.
#' @param factor Numeric vector. Compression factor for each interval — higher
#'   values squash more aggressively.
#' @param min_squashed_length Numeric. Minimum display width (in original
#'   coordinate units) for any compressed interval, preventing regions from
#'   collapsing to nothing. Default `100`.
#'
#' @return A `scales` transformation object suitable for passing to
#'   `scale_x_continuous(trans = ...)`.
#'
#' @noRd
improved_squash_axis <- function(from, to, factor, min_squashed_length = 100) {
  from <- as.vector(from)
  to <- as.vector(to)
  factor <- as.vector(factor)
  if (length(from) != length(to) || length(from) != length(factor)) {
    stop("from, to, and factor must be vectors of the same length")
  }
  # Ensure intervals are ordered left to right and don't overlap
  if (length(from) > 1) {
    ord <- order(from)
    from <- from[ord]
    to <- to[ord]
    factor <- factor[ord]
    for (i in 1:(length(from) - 1)) {
      if (!is.na(from[i]) && !is.na(to[i]) && !is.na(from[i+1]) && 
          to[i] >= from[i+1]) {
        stop("Regions cannot overlap")
      }
    }
  }
  # Pre-compute how many display units each squashed region saves.
  # Cumulated so that coordinates after multiple squashed regions
  # can all be shifted left by the correct total amount.
  reductions <- numeric(length(from))
  for (i in 1:length(from)) {
    if (!is.na(from[i]) && !is.na(to[i]) && !is.na(factor[i])) {
      original_width <- to[i] - from[i]
      squashed_width <- max(original_width / factor[i], min_squashed_length)
      reductions[i] <- original_width - squashed_width
    }
  }
  cumulative_reductions <- cumsum(reductions)
  
  # Forward transform: genomic -> display coordinates
  # Points inside a squashed region are rescaled to fit the shorter width;
  # points after it are shifted left by the cumulative reduction so far.
  trans <- function(x) {
    if (length(x) == 0 || all(is.na(x))) return(x)
    result <- x
    for (i in 1:length(from)) {
      if (is.na(from[i]) || is.na(to[i]) || is.na(factor[i])) next
      offset <- if (i > 1) cumulative_reductions[i-1] else 0
      idx_inside <- which(!is.na(x) & x >= from[i] & x <= to[i])
      if (length(idx_inside) > 0) {
        original_width <- to[i] - from[i]
        squashed_width <- max(original_width / factor[i], min_squashed_length)
        scaling_factor <- original_width / squashed_width
        rel_pos <- (x[idx_inside] - from[i]) / scaling_factor
        result[idx_inside] <- from[i] - offset + rel_pos
      }
      idx_after <- which(!is.na(x) & x > to[i])
      if (length(idx_after) > 0) {
        result[idx_after] <- x[idx_after] - cumulative_reductions[i]
      }
    }
    return(result)
  }
  # Inverse transform: display -> genomic coordinates
  # Required by ggplot2 to place axis tick labels at the correct genomic
  # positions. Processes regions in reverse order to correctly undo
  # cumulative shifts.
  inv <- function(x) {
    if (length(x) == 0 || all(is.na(x))) return(x)
    result <- x
    for (i in length(from):1) {
      if (is.na(from[i]) || is.na(to[i]) || is.na(factor[i])) next
      offset <- if (i > 1) cumulative_reductions[i-1] else 0
      original_width <- to[i] - from[i]
      squashed_width <- max(original_width / factor[i], min_squashed_length)
      scaling_factor <- original_width / squashed_width
      squashed_from <- from[i] - offset
      squashed_to <- squashed_from + squashed_width
      idx_after <- which(!is.na(x) & x >= squashed_to)
      if (length(idx_after) > 0) {
        result[idx_after] <- x[idx_after] + reductions[i]
      }
      idx_inside <- which(!is.na(x) & x >= squashed_from & x <= squashed_to)
      if (length(idx_inside) > 0) {
        rel_pos <- x[idx_inside] - squashed_from
        result[idx_inside] <- from[i] + rel_pos * scaling_factor
      }
    }
    return(result)
  }
  return(scales::trans_new("improved_squash_axis", trans, inv, domain = c(-Inf, Inf)))
}

#' Generate x-axis break positions for a gene plot
#'
#' Computes sensible genomic tick positions for the x-axis depending on the
#' requested verbosity mode. In `"all"` mode, break density scales with exon
#' size so that larger exons get more ticks while tiny exons are unlabelled.
#'
#' @param mode Character. One of:
#'   \describe{
#'     \item{`"none"`}{No tick marks.}
#'     \item{`"minimal"`}{Only the coding region start and end.}
#'     \item{`"all"`}{Gene boundaries, intron edges, and exon boundaries for
#'       exons large enough to be visible, with extra interior ticks for
#'       particularly wide exons.}
#'   }
#'   Default `"all"`.
#' @param coding_start Integer. Genomic start of the coding region.
#' @param coding_end Integer. Genomic end of the coding region.
#' @param cds_exons_dt A `data.table` of CDS exon coordinates with columns
#'   `start` and `end`.
#' @param long_introns A `data.table` of long introns with columns
#'   `intron_start` and `intron_end`, as produced inside `create_geneplot`.
#'
#' @return A sorted numeric vector of genomic positions, or `NULL` if
#'   `mode = "none"`.
#'
#' @noRd
generate_breaks <- function(mode = "all", coding_start, coding_end,
                            cds_exons_dt, long_introns) {
  if (mode == "none")    return(NULL)
  if (mode == "minimal") return(sort(unique(c(coding_start, coding_end))))
  
  # Start with gene boundaries and the edges of any compressed introns
  # (so ticks always appear either side of a "//" marker)
  positions <- c(coding_start, coding_end,
                 long_introns$intron_start, long_introns$intron_end)
  
  # Add exon boundary ticks, scaling with exon size:
  # - exons < 1/20th of total CDS length are too small to label
  # - exons > 1/5th of total CDS length get additional interior ticks
  cds_exons_dt[, exon_length := end - start]
  total_exon_length <- sum(cds_exons_dt$exon_length)
  for (i in 1:nrow(cds_exons_dt)) {
    exon_start <- cds_exons_dt$start[i]
    exon_end <- cds_exons_dt$end[i]
    exon_length <- exon_end - exon_start
    if (exon_length > total_exon_length / 20) {
      positions <- c(positions, exon_start, exon_end)
      if (exon_length > total_exon_length / 5) {
        n_ticks <- min(5, floor(exon_length / (total_exon_length / 10)))
        tick_positions <- exon_start + (1:n_ticks) * exon_length / (n_ticks + 1)
        positions <- c(positions, tick_positions)
      }
    }
  }
  return(sort(unique(round(positions))))
}



#' Create a gene lolliplot
#'
#' Plots variant-level (SNP) or domain-level association results as a lolliplot
#' overlaid on the canonical transcript structure of a gene, with optional
#' intron compression and protein domain annotations.
#'
#' @param gene_name Character. Gene symbol, ENSG, or ENST ID.
#' @param markers A `data.table` of variant or domain results, or `NULL`.
#'   For SNP mode must contain columns: `SYMBOL`, `POS`, `BETA`, `AC`, and the
#'   column named by `value_col`. For domain/burden mode must contain `SYMBOL`,
#'   `start`, `end`, `name`, `BETA`, and `value_col`.
#' @param marker_class Character. One of `"snp"`, `"burden"`, or `"domain"`.
#'   Default `"snp"`.
#' @param domains A `data.table` with columns `start`, `end`, `domain` (and
#'   optionally `colour`) for annotating protein domains. Default `NULL`.
#' @param value_col Character. Column name in `markers` to use as the y-axis
#'   value. Default `"P_BOLT_LMM_INF"`.
#' @param pvalue Logical. If `TRUE`, values are treated as p-values and
#'   transformed to signed -log10(p). Default `TRUE`.
#' @param title_use Character. Plot title. Default `"Gene Lolliplot"`.
#' @param subtitle_use Character or NULL. Plot subtitle.
#' @param colour_title Character. Legend title for colour mapping.
#' @param shape_title Character. Legned title for shape mapping.
#' @param min_exon_width Numeric or NULL. Minimum exon box half-height in data
#'   units. Auto-calculated if `NULL`.
#' @param intron_atten_cuttoff Numeric. Introns longer than this (bp) are
#'   compressed. Default `500`.
#' @param max_squashed_length Numeric or NULL. Maximum display width (in bp)
#'   for compressed introns. Auto-calculated if `NULL`.
#' @param marker_gene_buffer Numeric. Gap between gene body and lollipop stems.
#'   Default `0.4`.
#' @param exon_fill Character. Fill colour for CDS exon boxes. Default
#'   `"steelblue"`.
#' @param x_axis_labels Character. One of `"all"`, `"minimal"`, or `"none"`.
#'   Default `"all"`.
#' @param major_sig_thresh Numeric. Y-value for the major significance line
#'   (default genome-wide: `-log10(5e-8)`).
#' @param major_sig_colour Character. Colour for the major significance line.
#'   Default `"blue"`.
#' @param minor_sig_thresh Numeric. Y-value for the minor significance line
#'   (default: `-log10(0.05)`).
#' @param minor_sig_colour Character. Colour for the minor significance line.
#'   Default `"red"`.
#' @param plot_sig_lines Logical. Whether to draw significance threshold lines.
#'   Default `TRUE`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#'   # Basic gene structure only
#'   create_geneplot("POMC",intron_atten_cuttoff=NULL)
#'   
#'   # Compressed introns
#'   create_geneplot("POMC",intron_atten_cuttoff=500)
#'
#'   # With SNP results
#'   set.seed(42)
#'   pomc_snps <- data.table(
#'     SYMBOL = "POMC",
#'     POS    = c(sample(25161081:25162000, 40),sample(25164500:25164772, 10)),
#'     BETA   = rnorm(50, mean = 0.5, sd = 0.3),
#'     P_BOLT_LMM_INF = c(
#'       runif(45, 0.01, 1),      # background variants
#'       runif(4,  1e-6, 1e-4),   # suggestive hits
#'       5e-9                      # one genome-wide significant hit
#'     ),
#'     AC = sample(c(1, 5, 20, 100, 500, 2000), 50, replace = TRUE)
#'   )
#'   create_geneplot("POMC", markers = pomc_snps, value_col = "P_BOLT_LMM_INF", 
#'   title_use = "POMC burden")
#'   # Add on different colours and shapes
#'   pomc_snps[,class:=sample(c("ptv","missense"), 50, replace = TRUE)]
#'   create_geneplot("POMC", markers = pomc_snps, value_col = "P_BOLT_LMM_INF", 
#'   title_use = "POMC burden")
#'   pomc_snps[,`:=`(
#'     colour=fifelse(class=="ptv","red","orange"),
#'     colour_key=class
#'   )]
#'   create_geneplot("POMC", markers = pomc_snps, value_col = "P_BOLT_LMM_INF", 
#'   title_use = "POMC burden",colour_title="Class")
#'   pomc_snps[,`:=`(
#'     shape=fifelse(P_BOLT_LMM_INF < 0.05,16,1),
#'     shape_key=fifelse(P_BOLT_LMM_INF < 0.05,"Sig","Insig")
#'   )]
#'   create_geneplot("POMC", markers = pomc_snps, value_col = "P_BOLT_LMM_INF", 
#'   title_use = "POMC burden",colour_title="Class",shape_title="Significance")
#'
#'   # Overlaying protein domains on SNP results
#'   pomc_protein_domains <- get_protein_features("POMC", type_limit = "Peptide", 
#'   map_to_genomic = TRUE)[,
#'     .(domain=description,start=genomic_start,end=genomic_end)][!is.na(start)]
#'   ]
#'   create_geneplot("POMC", markers = pomc_snps, value_col = "P_BOLT_LMM_INF",
#'                   domains = pomc_protein_domains, 
#'                   title_use = "POMC SNPs + domains")
#'                   
#'   # With domain results
#'   set.seed(42)
#'   pomc_domains <- data.table(
#'     SYMBOL = "POMC",
#'     name   = c("Signal peptide", "NPP", "ACTH", "alpha-MSH", "beta-MSH"),
#'     start  = c(25161081, 25161200, 25161500, 25161500, 25163000),
#'     end    = c(25161150, 25161450, 25161800, 25161600, 25163500),
#'     BETA   = c(1, 1, -1, 1, -1),
#'     P_BOLT_LMM_INF = c(1e-9, 0.02, 1e-4, 0.5, 0.03)
#'   )
#'   create_geneplot("POMC", markers = pomc_domains, marker_class = "burden",
#'                   value_col = "P_BOLT_LMM_INF", 
#'                   title_use = "POMC domain burden")
#'   # With custom domain colours
#'   pomc_domains[, colour := c("red", "orange", "steelblue", "pink", "green")]
#'   create_geneplot("POMC", markers = pomc_domains, marker_class = "burden",
#'                   value_col = "P_BOLT_LMM_INF", 
#'                   title_use = "POMC domain burden")
#'
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_rect geom_point geom_hline
#'   scale_x_continuous scale_y_continuous scale_colour_identity
#'   scale_colour_manual scale_colour_discrete scale_fill_manual
#'   scale_fill_discrete scale_size_manual guide_legend guides labs
#'   theme_minimal theme element_line element_text element_blank annotate
#' @importFrom data.table as.data.table data.table setorder rbindlist
#' @importFrom scales trans_new
#' @export
create_geneplot <- function(gene_name=NULL,markers = NULL,
                            marker_class = "snp", domains = NULL,
                            value_col = "P_BOLT_LMM_INF",
                            pvalue = T,
                            title_use = "Gene Lolliplot",
                            subtitle_use = NULL,
                            colour_title = NULL,
                            shape_title = NULL,
                            min_exon_width=NULL,
                            intron_atten_cuttoff=500,
                            max_squashed_length=NULL,
                            marker_gene_buffer = 0.4,
                            exon_fill="steelblue",
                            x_axis_labels = "all",
                            major_sig_thresh = -log10(5e-08),
                            major_sig_colour = "blue",
                            minor_sig_thresh = -log10(0.05),
                            minor_sig_colour = "red",
                            plot_sig_lines = TRUE){
  ##### Gene info #####
  # Query the base coords for a genes exons and length of coding sequence, calculating introns based on gaps
  if(is.null(gene_name)) stop("Please provide a gene")
  # Get raw coords for introns and exons
  gene_invest_info <- ensembl_gene_lookup(input_id=gene_name)
  coding_start <- gene_invest_info$CODING_INFO$start
  coding_end <- gene_invest_info$CODING_INFO$end
  strand <- gene_invest_info$STRAND
  exon_dt <- as.data.table(gene_invest_info$EXONS)
  # Make the coding line
  coding_line <- data.frame(
    start = coding_start,
    end = coding_end
  )
  # Classify exons based on coding sequence length
  split_exons_list <- lapply(1:nrow(exon_dt), function(i) {
    split_exon(exon_dt$start[i], exon_dt$end[i], coding_start, coding_end)
  })
  split_exons_dt <- rbindlist(split_exons_list, idcol = "exon_number")
  # Filter for just CDS's
  cds_exons_dt <- split_exons_dt[region == "CDS"]
  setorder(cds_exons_dt, start)
  # Compute intron start/end between exons (remove last row as no intron after last exon)
  intron_dt <- cds_exons_dt[, .(
    intron_start = end,
    intron_end = data.table::shift(start, type = "lead")
  )][!is.na(intron_start) & !is.na(intron_end)]
  
  
  ##### Intron compression setup #####
  # Long introns are visually problematic - they push exons far apart and waste plot space.
  # We identify which introns exceed the threshold and will need compressing on the axis.
  # Calculate intron length
  intron_dt[, intron_length := intron_end - intron_start]
  # Flag if we need to collapse any introns
  long_introns <- intron_dt[intron_length > intron_atten_cuttoff]
  from_vec <- long_introns$intron_start
  to_vec <- long_introns$intron_end
  apply_squash <- nrow(long_introns) > 0
  short_introns <- intron_dt[intron_length < intron_atten_cuttoff]
  # The compressed display width for long introns is set to match the shortest
  # real intron (so nothing looks out of place), or a small fraction of the
  # total gene width if there are no short introns to reference.
  # NOTE - this can be overridden directly via max_squashed_length. 
  if(nrow(short_introns)!=0){
    smallest_short_intron_length <- min(intron_dt[intron_length < 500]$intron_length)
  } else {
    # Set to 5% total width
    smallest_short_intron_length <- (coding_end-coding_start)/200
  }
  # Override with max_squashed_length if availible
  if(!is.null(max_squashed_length)) smallest_short_intron_length <- max_squashed_length
  # Set factor for squashing introns
  if (apply_squash) {
    if (!is.null(max_squashed_length)) {
      factors <- long_introns$intron_length / max_squashed_length
    } else {
      factors <- pmin(long_introns$intron_length / 100, 20)
      factors <- pmax(factors, 5)
    }
  }
  
  ##### Misc helpers and saftey setups #####
  # Sort out min exon width
  if (is.null(min_exon_width)){
    min_exon_width <- 1
  }
  round_up_to <- function(x, multiple = 10) ceiling(x / multiple) * multiple
  round_down_to <- function(x, base = 10) base * floor(x / base)
  
  
  ##### Plotting Code #####
  # Initial empty plot 
  p <- ggplot()
  
  #### Plot variants/domains ####
  # If variant data was provided, add it to the plot
  if (!is.null(markers)) {
    ### Handle domains ###
    if(marker_class %in% c("burden","domain")){
      cat("\nDomain level results detected\n")
      # Curtail any domain regions that extend beyond the coding sequence boundaries
      marker_df <- markers[
        SYMBOL == gene_invest_info$SYMBOL &
          !is.na(get(value_col)) &
          end >= gene_invest_info$CODING_INFO$start &
          start <= gene_invest_info$CODING_INFO$end
      ][, `:=`(
        start = pmax(start, gene_invest_info$CODING_INFO$start),
        end = pmin(end, gene_invest_info$CODING_INFO$end)
      )]
      # Convert p-values to signed -log10 scale, or use raw value
      if(pvalue){
        marker_df[, plot_val := (-log10(get(value_col))) * BETA]
      } else {
        marker_df[, plot_val := get(value_col)]
      }
      # Scale exon box height to 5% of the data range so domains sit visibly above the gene body
      min_exon_width <-  max(round_down_to(max(marker_df$plot_val, na.rm = TRUE), 10) * 0.05,1)
    } else {
      ### Handle variants ###
      cat("\nSNP level results detected/assumed\n")
      # Keep only variants that fall within the coding region of the target gene
      marker_df <- markers[SYMBOL == gene_invest_info$SYMBOL & !is.na(get(value_col)) & POS > coding_start & POS < coding_end]
      # Convert p-values to signed -log10 scale, or use raw value
      if(pvalue){
        marker_df[, plot_val := (-log10(get(value_col))) * sign(BETA)]
      } else {
        marker_df[, plot_val := get(value_col)]
      }
      # Bucket allele count into display bands for point sizing
      ac_levels <- c("1", "2+", "10+", "50+", "100+", "500+", "1000+")
      marker_df[, AC_BIN := factor(cut(
        AC,
        breaks = c(-Inf, 1, 10, 50, 100, 500, 1000, Inf),
        labels = ac_levels,
        right = TRUE
      ), levels = ac_levels)]
    }
    # Shift each marker value away from zero by min_offset so it clears the gene body
    min_offset <- min_exon_width + marker_gene_buffer
    marker_df[, buffered_val:= fifelse(plot_val<0,plot_val-min_offset,plot_val + min_offset)]
    
    ### Y-axis labels ###
    # Compute pretty axis limits from the data range, then nudge them outward slightly if the data sits too close to a rounded boundary (avoids clipping points)
    # Calculate y-labels based on original axis:
    min_p <- min(marker_df$plot_val, na.rm = TRUE)
    max_p <- max(marker_df$plot_val, na.rm = TRUE)
    # Round up/down to nearest 5 and check they aren't too close to true max
    y_limit_top <- round_up_to(max_p, 5)
    if ((y_limit_top - max_p) < 1) y_limit_top <- y_limit_top + 5
    y_limit_bottom <- round_down_to(min_p, 5)
    if ((min_p - y_limit_bottom) < 1) y_limit_bottom <- y_limit_bottom - 5
    labels_pretty <- pretty(c(y_limit_bottom, y_limit_top), n = 7)
    # Handle one-sided plots (all effects in same direction) by dropping the
    # irrelevant half of the axis and inserting a blank spacer at the centre
    # where the gene body sits
    if (all(sign(marker_df$BETA) > 0)) {
      labels_pretty <- labels_pretty[labels_pretty >= 0]
      labels_full <- c(0,"",labels_pretty)
    } else if (all(sign(marker_df$BETA) < 0)) {
      labels_pretty <- labels_pretty[labels_pretty <= 0]
      labels_full <- c(labels_pretty,"",0)
    } else {
      labels_full <- append(labels_pretty, c(0, ""), after = which(labels_pretty == 0)[1] - 1)
    }
    # Remap label positions to account for the gene body gap (min_offset).
    # Labels below centre are shifted down, labels above centre are shifted up,
    # and the blank spacer maps to exactly zero (the gene body baseline).
    center_idx <- which(labels_full == "")
    breaks_needed <- sapply(seq_along(labels_full), function(i) {
      val <- labels_full[i]
      if (val == "") return(0)
      adj <- as.numeric(val)
      if (i < center_idx) adj - min_offset
      else if (i > center_idx) adj + min_offset
      else adj  # exact center stays the same (already covered by "" -> 0)
    })
    
    ### Significance lines ##
    # Only draw lines if the threshold falls within the plotted y range (within 25% headroom).
    # This prevents lines appearing on plots where no variants come close to significance
    if(plot_sig_lines){
      cat("\nChecking if significance lines are suitable\n")
      # Check if either is within our pretty labels + 25% (if not then don't plot)
      if(!is.null(minor_sig_thresh)){
        if(minor_sig_thresh <= max(labels_pretty)*1.25){
          p <- p +  geom_hline(yintercept = minor_sig_thresh+min_offset, linetype = "dashed", colour = minor_sig_colour)
        }
        if(-minor_sig_thresh >= min(labels_pretty)*1.25){
          p <- p +  geom_hline(yintercept = -minor_sig_thresh-min_offset, linetype = "dashed", colour = minor_sig_colour)
        }
      }
      if(!is.null(major_sig_thresh)){
        if(major_sig_thresh <= max(labels_pretty)*1.25){
          p <- p +  geom_hline(yintercept = major_sig_thresh+min_offset, linetype = "dashed", colour = major_sig_colour)
        }
        if(-major_sig_thresh >= min(labels_pretty)*1.25){
          p <- p +  geom_hline(yintercept = -major_sig_thresh-min_offset, linetype = "dashed", colour = major_sig_colour)
        }
      }
    }
    
    ### Add on markers ###
    # Burden / domain regions
    if (marker_class %in% c("burden","domain")) {
      # Domain/burden results: filled rectangles spanning each region
      p <- p +
        geom_rect(
          data = marker_df[plot_val != 0],
          aes(
            xmin = start, xmax = end,
            ymin = min_offset * sign(plot_val),
            ymax = buffered_val,
            fill = name
          ),
          alpha = 0.7
        )
      if ("colour" %in% colnames(marker_df)) {
        colour_map <- setNames(marker_df$colour, marker_df$name)
        p <- p + scale_fill_manual(name = "Region", values = colour_map)
      } else {
        p <- p + scale_fill_discrete(name = "Region")
      }
    } else {
      # SNP results: lollipops
      # Determine colour aesthetic source in order of preference: user-supplied hex > class column > none
      if ("colour" %in% colnames(marker_df)) {
        colour_aes    <- quote(colour)
        colour_source <- "colour"
      } else if ("class" %in% colnames(marker_df)) {
        colour_aes    <- quote(class)
        colour_source <- "class"
        if(is.null(colour_title)) colour_title <- "Class"
      } else {
        colour_aes    <- "none"
        colour_source <- "none"
      }
      # Shape aesthetic - use if column present
      shape_present <- "shape" %in% colnames(marker_df)
      
      p <- p +
        geom_segment(
          data = marker_df[plot_val != 0],
          aes(
            x = POS, xend = POS,
            y    = ifelse(plot_val >= 0, min_offset, -min_offset),
            yend = buffered_val,
            colour = if (colour_aes != "none") !!colour_aes else NULL
          ),
          size = 0.3
        ) +
        geom_point(
          data = marker_df[plot_val != 0],
          aes(
            x = POS, y = buffered_val,
            size   = AC_BIN,
            colour = if (colour_aes != "none") !!colour_aes else NULL,
            shape  = if (shape_present) shape else NULL
          ),
          stroke = 0.2, alpha = 0.6
        ) +
        scale_size_manual(
          values = c("1" = 1, "2+" = 1.5, "10+" = 2, "50+" = 2.5, "100+" = 3, "500+" = 3.5, "1000+" = 4),
          guide  = guide_legend(title = "AC")
        )
      
      # Shape scale - if colour_key column present show legend, otherwise suppress it
      if (shape_present) {
        if ("shape_key" %in% colnames(marker_df)) {
          shape_key_map <- unique(marker_df[!is.na(shape) & !is.na(shape_key), .(shape, shape_key)])
          p <- p + scale_shape_identity(
            guide  = guide_legend(
              title = shape_title,
              override.aes = list(shape = shape_key_map$shape)
            ),
            labels = shape_key_map$shape_key,
            breaks = shape_key_map$shape
          )
        } else {
          # No key - shapes render correctly but no legend
          p <- p + scale_shape_identity(guide = "none")
        }
      }
      
      # Colour scale
      if (colour_source == "colour") {
        if ("colour_key" %in% colnames(marker_df)) {
          # colour_key provided - build legend using same approach as original class mapping
          marker_df$colour_key <- as.character(marker_df$colour_key)
          colour_key_map <- unique(marker_df[!is.na(colour_key) & !is.na(colour), .(colour_key, colour)])
          colour_mapping <- setNames(colour_key_map$colour_key, colour_key_map$colour)
          p <- p + scale_colour_manual(
            values     = names(colour_mapping),
            labels     = colour_mapping,
            breaks     = names(colour_mapping),
            aesthetics = "colour",
            guide      = guide_legend(title = colour_title)
          )
        } else {
          # No key - colours render correctly but no legend
          p <- p + scale_colour_identity(guide = "none")
        }
      } else if (colour_source == "class") {
        p <- p + scale_colour_discrete(name = colour_title)
      } else {
        p <- p + guides(colour = "none")
      }
    }
    # Add on horizontal lines marking top and bottom of gene-body to help seperate markers
    p <- p +
      geom_hline(yintercept = c(-min_offset, min_offset), 
                 linetype = "solid", colour = "black", linewidth = 0.4) + 
      scale_y_continuous(
        limits = c(min(breaks_needed),max(breaks_needed)),
        breaks = breaks_needed,
        labels = labels_full
      )
  }
  
  
  #### Plot gene body ####
  p <- p +
    # Thin grey baseline spanning the full coding region
    geom_segment(data = coding_line, aes(x = start, xend = end, y = 0, yend = 0),
                 size = 1, colour = "grey40") + 
    # Thick boxes for CDS exons
    geom_rect(data = cds_exons_dt,
              aes(xmin = start, xmax = end, ymin = -min_exon_width, ymax = min_exon_width),
              fill = exon_fill, colour = "black", linewidth = 0.3)
  
  ### Domains ###
  # Add domains if specified
  if (!is.null(domains) & all(c("start","end","domain") %in% colnames(domains))){
    if(marker_class %in% c("burden","domain")){
      warning("Skipping domain annotations since markers are for domain results already !")
    } else {
      # Clip domains to coding region boundaries
      domains_filtered <- domains[
        end >= gene_invest_info$CODING_INFO$start &
          start <= gene_invest_info$CODING_INFO$end
      ][, `:=`(
        start = pmax(start, gene_invest_info$CODING_INFO$start),
        end = pmin(end, gene_invest_info$CODING_INFO$end)
      )]
      # Compute display height for each domain. Domains are sorted by total width
      # (widest drawn tallest) and overlapping domains are stepped down in height
      # so they remain distinguishable rather than being drawn on top of each other
      domains_filtered[, seg_width := end - start]
      domain_sizes <- domains_filtered[, .(
        total_width = sum(seg_width),
        span_start = min(start),
        span_end = max(end)
      ), by = domain]
      setorder(domain_sizes, -total_width)
      domain_sizes[, height_factor := NA_real_]
      placed_spans <- list()
      for (i in seq_len(nrow(domain_sizes))) {
        this <- domain_sizes[i]
        hf <- 0.9 # start at 90% of exon height to distinguish domain from exons themselves
        for (prev in placed_spans) {
          # Check overlap by full span
          if (this$span_start <= prev$span_end && this$span_end >= prev$span_start) {
            hf <- min(hf, prev$height_factor - 0.1) # step down if overlapping
          }
        }
        domain_sizes[i, height_factor := max(hf, 0.3)] # floor at 30% to stay visible
        placed_spans <- append(placed_spans, list(domain_sizes[i]))
      }
      domains_expanded <- merge(domains_filtered, domain_sizes[, .(domain, height_factor)], by = "domain", all.x = TRUE)
      # Preserve original domain order for the fill legend
      original_order <- unique(domains_filtered$domain)
      domains_expanded[, domain := factor(domain, levels = domain_sizes$domain)]
      domains_expanded[, domain_fill := factor(domain, levels = original_order)]
      setorder(domains_expanded, domain)
      
      # Add domains to plot
      p <- p + geom_rect(
        data = domains_expanded,
        aes(xmin = start, xmax = end, 
            ymin = -min_exon_width * height_factor, ymax = min_exon_width * height_factor,  # Slightly smaller than exons
            fill = domain_fill),
        linewidth = 2,
        alpha = 1
      ) 
      # Custom colours if specified
      if ("colour" %in% colnames(domains_filtered)) {
        colour_map <- setNames(domains_filtered$colour, domains_filtered$domain)
        p <- p + scale_fill_manual(name = "Domains", values = colour_map)
      } else {
        p <- p + scale_fill_discrete(name = "Domains")
      }
    }
  }
  
  #### X-axis scale ####
  if (apply_squash) {
    # Apply the custom compression transformation if needed
    p <- p +
      scale_x_continuous(
        trans = improved_squash_axis(from_vec, to_vec, factors, min_squashed_length = ifelse(is.null(max_squashed_length), smallest_short_intron_length, max_squashed_length)),
        breaks = generate_breaks(x_axis_labels, coding_start, coding_end, cds_exons_dt, long_introns),
        labels = if (x_axis_labels == "none") NULL else generate_breaks(x_axis_labels, coding_start, coding_end, cds_exons_dt, long_introns)
      )
    # Add break markers where introns are compressed
    for (i in seq_along(from_vec)) {
      p <- p + annotate("text",
                        x = (from_vec[i] + to_vec[i]) / 2,
                        y  = 0,
                        label = "//",
                        size  = 3,
                        color = "black")
    }
  } else {
    # Simply scale as is
    p <- p +
      scale_x_continuous(
        breaks = generate_breaks(x_axis_labels, coding_start, coding_end, cds_exons_dt, long_introns),
        labels = if (x_axis_labels == "none") NULL else generate_breaks(x_axis_labels, coding_start, coding_end, cds_exons_dt, long_introns)
      )
  }
  
  #### Labels and theme ####
  x_axis_label <- if (x_axis_labels == "none") NULL else if (apply_squash) "Genomic position (introns compressed)" else "Genomic position"
  y_axis_label <- if (pvalue) expression(Signed ~ -log[10](italic(P))) else value_col
  p <- p +
    labs(
      x = x_axis_label,
      y = y_axis_label,
      title = title_use,
      subtitle = subtitle_use
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.4),
      panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = if (x_axis_labels == "none") element_blank() else element_text(angle = 60, hjust = 1),
      axis.ticks.x = if (x_axis_labels == "none") element_blank() else element_line(),
      panel.grid.major.x = element_line(color = "lightgray", linetype = "dotted", size = 0.3)
    ) 
  return(p)
} 