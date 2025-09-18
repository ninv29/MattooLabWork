# ================================================================
# Final Combined SILAF + SILAC Analysis
# Robust POI Detection, Fixed Shapes, White Background, Comparative Plots
# ================================================================

library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)

# ---------------------------
# Load datasets
# ---------------------------
# SILAF_originaldata: Drosophila
# SILAC_final: Human

# ---------------------------
# Convert numeric columns
# ---------------------------
SILAF_originaldata$halflife    <- as.numeric(as.character(SILAF_originaldata$halflife))
SILAF_originaldata$MedianAbund <- as.numeric(as.character(SILAF_originaldata$MedianAbund))
SILAC_final$halflife <- as.numeric(as.character(SILAC_final$halflife))
SILAC_final$MedianAbund <- as.numeric(as.character(SILAC_final$MedianAbund))

# ---------------------------
# Robust POI detection
# ---------------------------
poi_patterns <- list(
  HYPE       = c("HYPE", "FicD"),
  BiP        = c("BiP", "GRP78", "HSPA5", "HspA5"),
  SCNA       = c("SCNA", "alpha\\-?Synuclein"),
  Lysosomal  = c("Lysos"),
  HTT        = c("HTT", "Huntingtin"),
  HSP70      = c("HSP70", "Hsp70"),
  PRPH2      = c("PRPH2"),
  RDH12      = c("RDH12"),
  ATF4       = c("ATF4"),
  Cytokine   = c("IL6", "IL-6", "TGFB", "TGF", "TGF-beta")
)

poi_shapes <- c(
  HYPE=17, BiP=15, SCNA=16, Lysosomal=18, HTT=8,
  HSP70=3, PRPH2=4, RDH12=25, ATF4=9, Cytokine=1
)

find_POI <- function(df, poi_patterns) {
  df$POI <- FALSE
  df$POI_category <- NA
  for (poi in names(poi_patterns)) {
    patterns <- poi_patterns[[poi]]
    combined_pattern <- paste(patterns, collapse="|")
    hits <- grepl(combined_pattern, df$Gene, ignore.case=TRUE) |
      grepl(combined_pattern, df$Protein.Name, ignore.case=TRUE) |
      grepl(combined_pattern, df$Description, ignore.case=TRUE)
    df$POI[hits] <- TRUE
    df$POI_category[hits] <- poi
  }
  return(df)
}

# ---------------------------
# Data cleaning and k calculations
# ---------------------------
clean_sil_data <- function(df, poi_patterns) {
  hl_cut <- quantile(df$halflife, 0.99, na.rm = TRUE)
  ab_cut <- quantile(df$MedianAbund, 0.99, na.rm = TRUE)
  
  df_clean <- df %>%
    filter((is.na(halflife) | halflife <= hl_cut)) %>%
    filter((is.na(MedianAbund) | MedianAbund <= ab_cut)) %>%
    mutate(
      k_deg = log(2)/halflife,
      k0    = MedianAbund * k_deg
    ) %>%
    find_POI(poi_patterns)
  
  return(df_clean)
}

SILAF_clean <- clean_sil_data(SILAF_originaldata, poi_patterns)
SILAC_clean <- clean_sil_data(SILAC_final, poi_patterns)

# ---------------------------
# Gene subsets & export
# ---------------------------
export_gene_subsets <- function(df, prefix) {
  fast <- df %>% filter(!is.na(halflife) & halflife < 5)
  write.csv(unique(fast$Gene), paste0("mattoolabwork/", prefix, "_fast_turnover_genes.csv"), row.names=FALSE)
  
  no <- df %>% filter((!is.na(k0) & k0 <= 0.02) | (!is.na(halflife) & halflife > 30))
  write.csv(unique(no$Gene), paste0("mattoolabwork/", prefix, "_no_turnover_genes.csv"), row.names=FALSE)
  
  low_ab <- quantile(df$MedianAbund, 0.25, na.rm=TRUE)
  late <- df %>% filter(!is.na(MedianAbund) & MedianAbund < low_ab) %>%
    filter((!is.na(k0) & k0 > 0.05) | (!is.na(halflife) & halflife < 10))
  write.csv(unique(late$Gene), paste0("mattoolabwork/", prefix, "_late_expressed_genes.csv"), row.names=FALSE)
  
  cat(paste0("Exported gene lists for ", prefix, ".\n"))
  return(list(fast=unique(fast$Gene), no=unique(no$Gene), late=unique(late$Gene)))
}

SILAF_subsets <- export_gene_subsets(SILAF_clean, "SILAF")
SILAC_subsets <- export_gene_subsets(SILAC_clean, "SILAC")

# ---------------------------
# GO enrichment
# ---------------------------
convert_to_entrez <- function(gene_list, orgdb) {
  if (identical(orgdb, org.Dm.eg.db)) {
    gene_list <- gsub("Dmel\\\\", "", gene_list)
    id_types <- c("SYMBOL","FLYBASECG")
    mapped_all <- data.frame()
    for (type in id_types) {
      unmapped <- setdiff(gene_list, mapped_all$INPUT)
      if(length(unmapped)==0) break
      mapped <- suppressWarnings(tryCatch(
        bitr(unmapped, fromType=type, toType="ENTREZID", OrgDb=orgdb), 
        error=function(e) NULL))
      if(!is.null(mapped) && nrow(mapped)>0){
        colnames(mapped)[1] <- "INPUT"
        mapped_all <- rbind(mapped_all, mapped)
      }
    }
    return(unique(mapped_all$ENTREZID))
  } else {
    return(bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)$ENTREZID)
  }
}

run_go <- function(fast, no, late, prefix, orgdb) {
  fast_entrez <- convert_to_entrez(fast, orgdb)
  no_entrez <- convert_to_entrez(no, orgdb)
  late_entrez <- convert_to_entrez(late, orgdb)
  
  go_fast <- enrichGO(fast_entrez, OrgDb=orgdb, keyType="ENTREZID",
                      ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
  go_no <- enrichGO(no_entrez, OrgDb=orgdb, keyType="ENTREZID",
                    ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
  go_late <- enrichGO(late_entrez, OrgDb=orgdb, keyType="ENTREZID",
                      ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
  
  write.csv(as.data.frame(go_fast), paste0("mattoolabwork/", prefix, "_GO_fast_turnover.csv"), row.names=FALSE)
  write.csv(as.data.frame(go_no), paste0("mattoolabwork/", prefix, "_GO_no_turnover.csv"), row.names=FALSE)
  write.csv(as.data.frame(go_late), paste0("mattoolabwork/", prefix, "_GO_late_expressed.csv"), row.names=FALSE)
  
  cat(paste0("GO enrichment saved for ", prefix, ".\n"))
}

run_go(SILAF_subsets$fast, SILAF_subsets$no, SILAF_subsets$late, "SILAF", org.Dm.eg.db)
run_go(SILAC_subsets$fast, SILAC_subsets$no, SILAC_subsets$late, "SILAC", org.Hs.eg.db)

# ---------------------------
# Volcano plot function with fixed POI shapes & white background
# ---------------------------
plot_volcano <- function(df, prefix) {
  volc_k0 <- df %>% dplyr::select(Gene, MedianAbund, k0, POI, POI_category)
  volc_kdeg <- df %>% dplyr::select(Gene, halflife, k_deg, POI, POI_category)
  
  theme_white <- theme_minimal(base_size=12) + 
    theme(panel.background=element_rect(fill="white", color="black"),
          plot.background=element_rect(fill="white", color="black"))
  
  # k0 volcano
  p1 <- ggplot(volc_k0, aes(x=log10(MedianAbund+1), y=log10(k0+1))) +
    geom_point(alpha=0.5, color="grey70") +
    geom_point(data=subset(volc_k0, POI==TRUE), aes(shape=POI_category, color=POI_category), size=3) +
    scale_shape_manual(values=poi_shapes) +
    geom_text_repel(data=subset(volc_k0, POI==TRUE), aes(label=Gene), size=3) +
    theme_white +
    labs(title=paste0("Volcano Plot: Synthesis (k0) - ", prefix),
         x="log10(Median Abundance+1)", y="log10(Synthesis rate k0+1)")
  
  ggsave(paste0("mattoolabwork/", prefix, "_volcano_k0.png"), p1, width=7, height=6, dpi=300, bg="white")
  write.csv(volc_k0, paste0("mattoolabwork/", prefix, "_volcano_k0_data.csv"), row.names=FALSE)
  
  # k_deg volcano
  p2 <- ggplot(volc_kdeg, aes(x=log10(halflife+1), y=log10(k_deg+1))) +
    geom_point(alpha=0.5, color="grey70") +
    geom_point(data=subset(volc_kdeg, POI==TRUE), aes(shape=POI_category, color=POI_category), size=3) +
    scale_shape_manual(values=poi_shapes) +
    geom_text_repel(data=subset(volc_kdeg, POI==TRUE), aes(label=Gene), size=3) +
    theme_white +
    labs(title=paste0("Volcano Plot: Degradation (k_deg) - ", prefix),
         x="log10(Halflife+1)", y="log10(Degradation rate k_deg+1)")
  
  ggsave(paste0("mattoolabwork/", prefix, "_volcano_kdeg.png"), p2, width=7, height=6, dpi=300, bg="white")
  write.csv(volc_kdeg, paste0("mattoolabwork/", prefix, "_volcano_kdeg_data.csv"), row.names=FALSE)
  
  cat(paste0("Volcano plots and data exported for ", prefix, ".\n"))
}

plot_volcano(SILAF_clean, "SILAF")
plot_volcano(SILAC_clean, "SILAC")

# ---------------------------
# Comparative analysis with combined legend
# ---------------------------
SILAF_POI <- SILAF_clean %>% filter(POI==TRUE) %>% mutate(Dataset="SILAF")
SILAC_POI <- SILAC_clean %>% filter(POI==TRUE) %>% mutate(Dataset="SILAC")

common_cols <- intersect(colnames(SILAF_POI), colnames(SILAC_POI))
combined_POI <- rbind(SILAF_POI[, common_cols], SILAC_POI[, common_cols])

theme_white <- theme_minimal(base_size=12) + 
  theme(panel.background=element_rect(fill="white", color="black"),
        plot.background=element_rect(fill="white", color="black"))

# Comparative volcano: k0 vs MedianAbund
p_comp_k0 <- ggplot(combined_POI, aes(x=log10(MedianAbund+1), y=log10(k0+1), color=Dataset, shape=POI_category)) +
  geom_point(size=3) +
  scale_shape_manual(values=poi_shapes) +
  geom_text_repel(aes(label=Gene), size=3) +
  theme_white +
  labs(title="Comparative Volcano: Synthesis (k0)",
       x="log10(Median Abundance+1)", y="log10(Synthesis rate k0+1)")

ggsave("mattoolabwork/Comparative_volcano_k0.png", p_comp_k0, width=8, height=6, dpi=300, bg="white")

# Comparative volcano: k_deg vs halflife
p_comp_kdeg <- ggplot(combined_POI, aes(x=log10(halflife+1), y=log10(k_deg+1), color=Dataset, shape=POI_category)) +
  geom_point(size=3) +
  scale_shape_manual(values=poi_shapes) +
  geom_text_repel(aes(label=Gene), size=3) +
  theme_white +
  labs(title="Comparative Volcano: Degradation (k_deg)",
       x="log10(Halflife+1)", y="log10(Degradation rate k_deg+1)")

ggsave("mattoolabwork/Comparative_volcano_kdeg.png", p_comp_kdeg, width=8, height=6, dpi=300, bg="white")

cat("Comparative analysis completed with combined legend and white background.\n")