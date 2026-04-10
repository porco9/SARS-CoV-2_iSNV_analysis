################################################################################
# Combined Figure A + B + C
#
# Figure A  (상단): 개체별 iSNV rate (callable sites normalized), NS+S stacked
# Figure B  (하단 좌): Dose × Tissue unique iSNV, NS+S stacked
# Figure C  (하단 우): Dose × Tissue unique iSNV sharing 분포
#
# 모든 strip text / axis text / legend 크기 통일
# patchwork: A / (B | C) 레이아웃
#
# Input 1: 10__최종정리본_per_callable_sites.xlsx  (Figure A)
# Input 2: 7__iSNV_summary...Syn분류.xlsx          (Figure B, C)
################################################################################

required_packages <- c("readxl","dplyr","purrr","ggplot2","stringr","patchwork")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 공통 설정 — 모든 figure 동일 적용
# ══════════════════════════════════════════════════════════════════════════════
STRIP_SIZE  <- 14   # Gene 패널 제목 (Spike, ORF1a 등)
AXIS_TEXT   <- 13    # y축 숫자
LEGEND_TEXT <- 13    # legend 글씨
AXIS_TITLE  <- 14    # y축 제목
DOSE_LABEL  <- 4.2  # 패널 내 dose 레이블 (geom_text size)

muttype_colors <- c("NS" = "#E74C3C", "S" = "#85C1E9")
sharing_colors <- c(
  "1 individual"   = "#FADBD8",
  "2 individuals"  = "#F1948A",
  "3+ individuals" = "#E74C3C"
)
genes_of_interest <- c("Spike","ORF1a","ORF1b","ORF3a","ORF8","ORF10")
group_levels      <- c("0.01\nURT","0.01\nLRT","0.1\nURT","1\nURT")

# 공통 테마 함수
theme_unified <- function(x_blank = TRUE) {
  base <- theme_classic(base_size = 12) +
    theme(text = element_text(family = "Times New Roman")) +
    theme(
      axis.line          = element_line(color = "black", linewidth = 0.5),
      axis.text.y        = element_text(size = AXIS_TEXT),
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(size = AXIS_TITLE, face = "bold"),
      strip.text         = element_text(face = "bold.italic", size = STRIP_SIZE),
      strip.background   = element_rect(fill = "gray95", color = NA),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size = LEGEND_TEXT),
      legend.key.size    = unit(0.35, "cm"),
      legend.spacing.x   = unit(0.15, "cm"),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.margin        = margin(6, 8, 6, 8),
      panel.spacing      = unit(0.7, "lines")
    )
  if (x_blank) {
    base <- base + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  } else {
    base <- base + theme(
      axis.text.x  = element_text(size = AXIS_TEXT, lineheight = 0.85),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5)
    )
  }
  base
}

# ══════════════════════════════════════════════════════════════════════════════
# 파일 선택
# ══════════════════════════════════════════════════════════════════════════════
cat("[1/2] Select callable sites Excel (10__최종정리본...)\n")
callable_file <- file.choose()
cat("[2/2] Select iSNV summary Excel (7__iSNV_summary...Syn분류.xlsx)\n")
isnv_file     <- file.choose()

out_dir <- dirname(callable_file)
out_pfx <- file.path(out_dir, "Combined_Figure_A_B_C")
cat("\nFiles loaded.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE A 데이터
# ══════════════════════════════════════════════════════════════════════════════
gene_ref_lengths <- c(
  Total=29903, Spike=3822, ORF1a=13218,
  ORF1b=8085,  ORF3a=828,  ORF8=366,  ORF10=117
)

make_colnames <- function(infile, sheet_name) {
  raw_top  <- readxl::read_excel(infile, sheet=sheet_name, col_names=FALSE, n_max=2)
  gene_row <- as.character(unlist(raw_top[1,]))
  for (i in 2:length(gene_row))
    if (is.na(gene_row[i])||gene_row[i]=="NA") gene_row[i] <- gene_row[i-1]
  sub_row   <- as.character(unlist(raw_top[2,]))
  col_names <- character(length(sub_row))
  for (j in seq_along(col_names)) {
    g <- gene_row[j]; s <- sub_row[j]
    if (is.na(s)||s=="NA")                       { col_names[j]<-paste0("Col",j); next }
    if (grepl("^Sample$",s,ignore.case=TRUE))     { col_names[j]<-"Sample";       next }
    if (grepl("^Group$",s,ignore.case=TRUE))      { col_names[j]<-"Group";        next }
    if (grepl("Total reads",s,ignore.case=TRUE))  { col_names[j]<-"Total_reads";  next }
    sub_clean <- s %>%
      stringr::str_replace_all("\\n"," ") %>%
      stringr::str_replace_all("iSNVs count","iSNVs") %>%
      stringr::str_replace_all("callable sites","callable") %>%
      stringr::str_replace_all("Coverage\\s*","coverage") %>%
      stringr::str_replace_all("iSNVs per\\s+callable sites","rate") %>%
      stringr::str_replace_all("iSNVs per callable sites","rate") %>%
      stringr::str_trim()
    col_names[j] <- paste0(g,"_",sub_clean)
  }
  col_names
}

parse_total_meta <- function(infile) {
  cnames   <- make_colnames(infile,"Total")
  raw_data <- readxl::read_excel(infile,sheet="Total",col_names=FALSE,skip=2)
  names(raw_data) <- cnames
  raw_data %>%
    dplyr::mutate(
      Dose   = as.numeric(as.character(Group)),
      Dose   = ifelse(is.na(Dose),0,Dose),
      Tissue = dplyr::case_when(
        stringr::str_detect(tolower(Sample),"lung")  ~ "LRT",
        stringr::str_detect(tolower(Sample),"nasal") ~ "URT",
        TRUE ~ "Unknown")
    ) %>%
    dplyr::filter(!is.na(Sample)) %>%
    dplyr::select(Sample,Dose,Tissue,Total_reads)
}

parse_sheet_long <- function(infile, sheet_name, meta, genes, gene_ref_lengths) {
  cnames   <- make_colnames(infile, sheet_name)
  raw_data <- readxl::read_excel(infile,sheet=sheet_name,col_names=FALSE,skip=2)
  names(raw_data) <- cnames
  raw_data <- raw_data %>%
    dplyr::mutate(Dose=as.numeric(as.character(Group)),
                  Dose=ifelse(is.na(Dose),0,Dose),
                  Total_reads=as.numeric(Total_reads)) %>%
    dplyr::filter(!is.na(Total_reads))
  if (!"Sample" %in% names(raw_data)) {
    n <- min(nrow(raw_data),nrow(meta))
    raw_data <- raw_data[seq_len(n),]
    raw_data$Sample      <- meta$Sample[seq_len(n)]
    raw_data$Total_reads <- meta$Total_reads[seq_len(n)]
  }
  genes_avail <- genes[sapply(genes, function(g)
    paste0(g,"_iSNVs") %in% names(raw_data) & paste0(g,"_callable") %in% names(raw_data))]
  purrr::map_dfr(genes_avail, function(g) {
    ref_len <- gene_ref_lengths[[g]]
    raw_data %>%
      dplyr::transmute(
        Sample=as.character(Sample), Dose=Dose, MutType=sheet_name, Gene=g,
        Count=as.numeric(.data[[paste0(g,"_iSNVs")]]),
        Callable=as.numeric(.data[[paste0(g,"_callable")]]),
        Total_reads=as.numeric(Total_reads),
        Callable_norm=Callable/ref_len,
        Rate=ifelse(Callable>0,
                    as.numeric(.data[[paste0(g,"_iSNVs")]])/(Callable/ref_len)*1000,
                    NA_real_)
      ) %>% dplyr::filter(is.finite(Callable),Callable>0)
  })
}

cat("Parsing callable sites...\n")
meta     <- parse_total_meta(callable_file)
all_long <- purrr::map_dfr(c("NS","S"),
  ~ parse_sheet_long(callable_file,.x,meta,genes_of_interest,gene_ref_lengths)) %>%
  dplyr::left_join(meta %>% dplyr::select(Sample,Tissue), by="Sample")

# Figure A 보조 데이터
sample_order_a <- meta %>% dplyr::arrange(Dose,Sample) %>% dplyr::pull(Sample) %>% unique()

dose_vlines_a <- meta %>%
  dplyr::arrange(Dose,Sample) %>%
  dplyr::mutate(x_idx=seq_len(dplyr::n())) %>%
  dplyr::group_by(Dose) %>%
  dplyr::summarise(xmax=max(x_idx)+0.5,.groups="drop") %>%
  dplyr::slice(-dplyr::n()) %>% dplyr::pull(xmax)

dose_labels_a <- meta %>%
  dplyr::arrange(Dose,Sample) %>%
  dplyr::mutate(x_idx=seq_len(dplyr::n())) %>%
  dplyr::group_by(Dose) %>%
  dplyr::summarise(x_mid=mean(x_idx),.groups="drop") %>%
  dplyr::mutate(label=dplyr::case_when(
    as.numeric(as.character(Dose))==0~"Pos", TRUE~as.character(Dose)))

df_a <- all_long %>%
  dplyr::mutate(
    Sample  = factor(Sample, levels=sample_order_a),
    x_idx   = as.integer(factor(Sample, levels=sample_order_a)),
    MutType = factor(MutType, levels=c("S","NS")),
    Gene    = factor(Gene, levels=genes_of_interest)
  )

p_a <- ggplot() +
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  geom_col(data=df_a, aes(x=x_idx,y=Rate,fill=MutType), width=0.75, color=NA) +
  scale_fill_manual(values=muttype_colors,
                    labels=c("NS"="Non-synonymous","S"="Synonymous"),
                    breaks=c("NS","S")) +
  geom_vline(xintercept=dose_vlines_a, color="black", linewidth=0.5, linetype="dashed") +
  geom_text(data=dose_labels_a,
            aes(x=x_mid, y=Inf, label=label),
            vjust=1.4, size=DOSE_LABEL, color="black", fontface="bold", inherit.aes=FALSE) +
  facet_wrap(~Gene, ncol=3, scales="free_y") +
  scale_x_continuous(breaks=NULL, expand=expansion(add=0.5)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(y="iSNVs per 1,000 normalized callable sites") +
  theme_unified(x_blank=TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE B + C 데이터
# ══════════════════════════════════════════════════════════════════════════════
parse_isnv_sheet <- function(isnv_file, sheet_name) {
  raw <- readxl::read_excel(isnv_file, sheet=sheet_name, col_names=FALSE)
  h1  <- as.character(unlist(raw[1,]))
  h2  <- as.character(unlist(raw[2,]))
  ref_col    <- which(h1=="Reference")[1]
  major_cols <- which(h2=="Major")
  results <- list()
  for (data_row in 3:nrow(raw)) {
    row_vals <- as.character(unlist(raw[data_row,]))
    pos <- suppressWarnings(as.numeric(row_vals[ref_col]))
    if (is.na(pos)) next
    for (mc in major_cols) {
      sample_name <- NA
      for (k in mc:1) {
        v <- h1[k]
        if (!is.na(v)&&v!="NA"&&
            !v %in% c("Reference","Ref_codon","Ref_AA","Codon#","nt ","NA")) {
          sample_name <- trimws(v); break
        }
      }
      if (is.na(sample_name)) next
      minor_val <- row_vals[mc+1]; syn_val <- row_vals[mc+2]
      if (is.na(minor_val)||minor_val=="NA") next
      if (is.na(syn_val)||!syn_val %in% c("NS","S")) next
      results[[length(results)+1]] <- data.frame(
        Gene=sheet_name, Position=pos, Sample=sample_name,
        Syn_Status=syn_val, stringsAsFactors=FALSE)
    }
  }
  if (length(results)==0) return(NULL)
  dplyr::bind_rows(results)
}

cat("Parsing iSNV summary...\n")
all_isnv <- purrr::map_dfr(genes_of_interest, function(g) {
  cat(" ",g,"...")
  df <- tryCatch(parse_isnv_sheet(isnv_file,g), error=function(e) NULL)
  if (is.null(df)) { cat("SKIP\n"); return(NULL) }
  cat(nrow(df),"rows\n"); df
}) %>%
  dplyr::mutate(
    Dose=dplyr::case_when(
      stringr::str_detect(Sample,"_0\\.01_") ~ "0.01",
      stringr::str_detect(Sample,"_0\\.1_")  ~ "0.1",
      stringr::str_detect(Sample,"_1_")      ~ "1",
      TRUE ~ NA_character_),
    Tissue=dplyr::case_when(
      stringr::str_detect(tolower(Sample),"lung")  ~ "LRT",
      stringr::str_detect(tolower(Sample),"nasal") ~ "URT",
      TRUE ~ "Unknown")
  ) %>%
  dplyr::filter(!is.na(Dose), Tissue!="Unknown")

# Figure B data
b_data <- all_isnv %>%
  dplyr::distinct(Dose,Tissue,Gene,Position,Syn_Status) %>%
  dplyr::count(Dose,Tissue,Gene,Syn_Status, name="Unique_iSNV") %>%
  dplyr::mutate(
    Group      = paste0(Dose,"\n",Tissue),
    Group      = factor(Group, levels=group_levels),
    Syn_Status = factor(Syn_Status, levels=c("S","NS")),
    Gene       = factor(Gene, levels=genes_of_interest)
  )

p_b <- ggplot(b_data, aes(x=Group, y=Unique_iSNV, fill=Syn_Status)) +
  geom_col(width=0.65, color=NA) +
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  scale_fill_manual(values=muttype_colors,
                    labels=c("NS"="Non-synonymous","S"="Synonymous"),
                    breaks=c("NS","S")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  facet_wrap(~Gene, ncol=3, scales="free_y") +
  labs(y="Number of unique iSNVs") +
  theme_unified(x_blank=FALSE)

# Figure C data
sharing_raw <- all_isnv %>%
  dplyr::distinct(Dose,Tissue,Gene,Position,Sample) %>%
  dplyr::count(Dose,Tissue,Gene,Position, name="N_individuals") %>%
  dplyr::mutate(
    Sharing=dplyr::case_when(
      N_individuals==1 ~ "1 individual",
      N_individuals==2 ~ "2 individuals",
      N_individuals>=3 ~ "3+ individuals"),
    Sharing=factor(Sharing, levels=c("1 individual","2 individuals","3+ individuals"))
  )

c_data <- sharing_raw %>%
  dplyr::count(Dose,Tissue,Gene,Sharing, name="Unique_iSNV") %>%
  dplyr::mutate(
    Group=paste0(Dose,"\n",Tissue),
    Group=factor(Group, levels=group_levels),
    Gene =factor(Gene, levels=genes_of_interest)
  )

p_c <- ggplot(c_data, aes(x=Group, y=Unique_iSNV, fill=Sharing)) +
  geom_col(width=0.65, color=NA) +
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  scale_fill_manual(values=sharing_colors, name="Shared in") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  facet_wrap(~Gene, ncol=3, scales="free_y") +
  labs(y="Number of unique iSNVs") +
  theme_unified(x_blank=FALSE) +
  theme(legend.title=element_blank())

# ══════════════════════════════════════════════════════════════════════════════
# COMBINED (patchwork)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nCombining with patchwork...\n")

combined <- p_a / (p_b | p_c) +
  patchwork::plot_layout(heights=c(1.5, 1))

fig_out <- paste0(out_pfx, ".png")
ggsave(fig_out, combined, width=18, height=13, dpi=300, bg="white")
cat("Saved:", basename(fig_out), "\n")

# 개별 저장
ggsave(paste0(out_pfx,"_A.png"), p_a, width=18, height=6.5, dpi=300, bg="white")
ggsave(paste0(out_pfx,"_B.png"), p_b, width=10, height=6,   dpi=300, bg="white")
ggsave(paste0(out_pfx,"_C.png"), p_c, width=10, height=6,   dpi=300, bg="white")
cat("Individual figures also saved.\n")

cat("\n✓ Done!\n")
cat("  Strip size :", STRIP_SIZE, "pt\n")
cat("  Axis text  :", AXIS_TEXT,  "pt\n")
cat("  Legend text:", LEGEND_TEXT,"pt\n")