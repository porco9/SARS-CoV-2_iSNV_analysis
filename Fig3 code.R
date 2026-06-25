################################################################################
# Combined Figure A + B + C  [REVISED]
#
# 변경사항:
#   1) Strip text 기울기(italic) 제거 → bold만 유지
#   2) Figure A X축: 각 샘플 바 아래에 "1st URT" / "2nd URT" / "1st LRT" 표기
#   3) Figure B, C X축 레이블: "1st URT", "2nd URT", "1st LRT" 형식
#
# Input 1: 10__최종정리본...xlsx  (Figure A) — A열=Trial, B열=Sites 포함
# Input 2: 7__iSNV_summary...Syn분류.xlsx  (Figure B, C)
################################################################################

required_packages <- c("readxl","dplyr","purrr","ggplot2","stringr","patchwork","tidyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 공통 설정
# ══════════════════════════════════════════════════════════════════════════════
STRIP_SIZE  <- 14
AXIS_TEXT   <- 11
LEGEND_TEXT <- 13
AXIS_TITLE  <- 14
DOSE_LABEL  <- 4.2

muttype_colors <- c("NS" = "#E74C3C", "S" = "#85C1E9")
sharing_colors <- c(
  "1 individual"   = "#FADBD8",
  "2 individuals"  = "#F1948A",
  "3+ individuals" = "#E74C3C"
)
genes_of_interest <- c("Spike","ORF1a","ORF1b","ORF3a","ORF8","ORF10")

# ── [수정 1] strip.text에서 italic 제거 ──────────────────────────────────────
theme_unified <- function(x_blank = TRUE) {
  base <- theme_classic(base_size = 12) +
    theme(text = element_text(family = "Times New Roman")) +
    theme(
      axis.line          = element_line(color = "black", linewidth = 0.5),
      axis.text.y        = element_text(size = AXIS_TEXT),
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(size = AXIS_TITLE, face = "bold"),
      # ★ bold만 유지, italic 제거
      strip.text         = element_text(face = "bold", size = STRIP_SIZE),
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
      axis.text.x  = element_text(size = AXIS_TEXT, angle = 45, hjust = 1, lineheight = 0.85),
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
    if (is.na(s)||s=="NA")                          { col_names[j]<-paste0("Col",j); next }
    if (grepl("^Sample$",s,ignore.case=TRUE))        { col_names[j]<-"Sample";        next }
    if (grepl("^Trials?$",s,ignore.case=TRUE))       { col_names[j]<-"Trial";         next }
    if (grepl("^Sites?$",s,ignore.case=TRUE))        { col_names[j]<-"Sites";         next }
    if (grepl("^Group$",s,ignore.case=TRUE))         { col_names[j]<-"Group";         next }
    if (grepl("Total reads",s,ignore.case=TRUE))     { col_names[j]<-"Total_reads";   next }
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

# ── [수정 2] parse_total_meta: Trial과 Sites 컬럼 추가 파싱 ─────────────────
parse_total_meta <- function(infile) {
  cnames   <- make_colnames(infile,"Total")
  raw_data <- readxl::read_excel(infile,sheet="Total",col_names=FALSE,skip=2)
  names(raw_data) <- cnames

  raw_data <- raw_data %>%
    dplyr::mutate(
      Dose  = as.numeric(as.character(Group)),
      Dose  = ifelse(is.na(Dose), 0, Dose)
    ) %>%
    dplyr::filter(!is.na(Total_reads))

  # Trial 컬럼 처리
  if ("Trial" %in% names(raw_data)) {
    raw_data <- raw_data %>%
      dplyr::mutate(
        Trial_label = dplyr::case_when(
          stringr::str_detect(tolower(Trial), "1st|1 st|first|trial 1|trial1") ~ "1st",
          stringr::str_detect(tolower(Trial), "2nd|2 nd|second|trial 2|trial2") ~ "2nd",
          TRUE ~ stringr::str_extract(Trial, "\\d+")
        ),
        Trial_label = ifelse(is.na(Trial_label), "1st", paste0(Trial_label))
      )
  } else {
    # Trial 컬럼 없으면 행 순서 기반으로 그룹 내 자동 부여
    raw_data <- raw_data %>%
      dplyr::group_by(Dose) %>%
      dplyr::mutate(Trial_label = ifelse(dplyr::row_number() <= ceiling(dplyr::n()/2), "1st", "2nd")) %>%
      dplyr::ungroup()
  }

  # Sites 컬럼 처리 (URT/LRT)
  if ("Sites" %in% names(raw_data)) {
    raw_data <- raw_data %>%
      dplyr::mutate(
        Tissue = dplyr::case_when(
          stringr::str_detect(toupper(Sites), "URT") ~ "URT",
          stringr::str_detect(toupper(Sites), "LRT") ~ "LRT",
          TRUE ~ "Unknown"
        )
      )
  } else {
    raw_data <- raw_data %>%
      dplyr::mutate(
        Tissue = dplyr::case_when(
          stringr::str_detect(tolower(Sample), "lung")  ~ "LRT",
          stringr::str_detect(tolower(Sample), "nasal") ~ "URT",
          TRUE ~ "Unknown"
        )
      )
  }

  # X축 레이블: "1st URT", "2nd LRT" 형식
  raw_data <- raw_data %>%
    dplyr::mutate(
      XLabel = paste0(Trial_label, " ", Tissue)
    )

  if (!"Sample" %in% names(raw_data)) {
    raw_data <- raw_data %>%
      dplyr::mutate(Sample = paste0("S", dplyr::row_number()))
  }

  raw_data %>%
    dplyr::select(Sample, Dose, Tissue, Trial_label, XLabel, Total_reads) %>%
    dplyr::filter(!is.na(Sample))
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
    n <- min(nrow(raw_data), nrow(meta))
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
                    as.numeric(.data[[paste0(g,"_iSNVs")]])/Callable*1000,
                    NA_real_)
      ) %>% dplyr::filter(is.finite(Callable), Callable>0)
  })
}

cat("Parsing callable sites...\n")
meta     <- parse_total_meta(callable_file)
all_long <- purrr::map_dfr(c("NS","S"),
  ~ parse_sheet_long(callable_file,.x,meta,genes_of_interest,gene_ref_lengths)) %>%
  dplyr::left_join(meta %>% dplyr::select(Sample, Tissue, Trial_label, XLabel), by="Sample")

cat("Sample meta preview:\n")
print(meta %>% dplyr::select(Sample, Dose, Tissue, Trial_label, XLabel))

# ── Tissue 배경색 설정 ────────────────────────────────────────────────────────
# Option 1: 파란 계열(URT) / 복숭아 계열(LRT) — 막대색과 간섭 최소
TISSUE_COLORS <- c("URT" = "#D6EAF8", "LRT" = "#FDEBD0")
TISSUE_ALPHA  <- 0.45   # 투명도: 낮을수록 연하게 (0~1)

# ── Figure A 보조 데이터 ─────────────────────────────────────────────────────
# 샘플 정렬: Dose → Trial → Tissue
sample_order_a <- meta %>%
  dplyr::arrange(Dose, Trial_label, Tissue, Sample) %>%
  dplyr::pull(Sample) %>% unique()

sample_xidx <- tibble::tibble(Sample = sample_order_a) %>%
  dplyr::mutate(x_idx = dplyr::row_number()) %>%
  dplyr::left_join(meta %>% dplyr::select(Sample, Dose, Tissue, XLabel), by="Sample")

# dose 구분선 위치
dose_vlines_a <- sample_xidx %>%
  dplyr::group_by(Dose) %>%
  dplyr::summarise(xmax=max(x_idx)+0.5, .groups="drop") %>%
  dplyr::slice(-dplyr::n()) %>% dplyr::pull(xmax)

# dose 레이블 (패널 상단)
dose_labels_a <- sample_xidx %>%
  dplyr::group_by(Dose) %>%
  dplyr::summarise(x_mid=mean(x_idx), .groups="drop") %>%
  dplyr::mutate(label=dplyr::case_when(
    as.numeric(as.character(Dose))==0 ~ "Pos", TRUE ~ as.character(Dose)))

# ── Tissue 배경 rect 데이터 (Figure A용) ────────────────────────────────────
# 연속된 같은 Tissue 블록을 하나의 rect로 합침
tissue_bg_a <- sample_xidx %>%
  dplyr::arrange(x_idx) %>%
  dplyr::mutate(
    block = cumsum(Tissue != dplyr::lag(Tissue, default = dplyr::first(Tissue)))
  ) %>%
  dplyr::group_by(block, Tissue) %>%
  dplyr::summarise(xmin = min(x_idx) - 0.5, xmax = max(x_idx) + 0.5, .groups = "drop") %>%
  dplyr::mutate(fill_color = TISSUE_COLORS[Tissue])

# ── [수정 3] X축 레이블: 각 바 아래 "1st URT" 형식 ─────────────────────────
xaxis_labels <- sample_xidx %>%
  dplyr::arrange(x_idx) %>%
  dplyr::pull(XLabel)

df_a <- all_long %>%
  dplyr::left_join(sample_xidx %>% dplyr::select(Sample, x_idx), by="Sample") %>%
  dplyr::mutate(
    Sample  = factor(Sample, levels=sample_order_a),
    MutType = factor(MutType, levels=c("S","NS")),
    Gene    = factor(Gene, levels=genes_of_interest)
  )

p_a <- ggplot() +
  # ★ Tissue 배경: annotate로 직접 색 지정해 fill scale 충돌 방지
  purrr::pmap(as.list(tissue_bg_a), function(block, Tissue, xmin, xmax, fill_color) {
    annotate("rect", xmin=xmin, xmax=xmax,
             ymin=-Inf, ymax=Inf, fill=fill_color, alpha=TISSUE_ALPHA)
  }) +
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  geom_col(data=df_a, aes(x=x_idx, y=Rate, fill=MutType), width=0.75, color=NA) +
  scale_fill_manual(values=muttype_colors,
                    labels=c("NS"="Non-synonymous","S"="Synonymous"),
                    breaks=c("NS","S")) +
  geom_vline(xintercept=dose_vlines_a, color="black", linewidth=0.5, linetype="dashed") +
  geom_text(data=dose_labels_a,
            aes(x=x_mid, y=Inf, label=label),
            vjust=1.4, size=DOSE_LABEL, color="black", fontface="bold", inherit.aes=FALSE) +
  facet_wrap(~Gene, ncol=3, scales="free_y") +
  scale_x_continuous(
    breaks = seq_along(xaxis_labels),
    labels = xaxis_labels,
    expand = expansion(add=0.5)
  ) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(y="iSNVs per 1,000 callable sites") +
  theme_unified(x_blank=FALSE) +
  theme(
    axis.text.x = element_text(size=8, angle=45, hjust=1, lineheight=0.85)
  )

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
      TRUE ~ "Unknown"),
    # ── [수정 4] Trial 정보 추출 (샘플명 앞 숫자: 1_ = 1st, 2_ = 2nd) ──────
    Trial_label = dplyr::case_when(
      stringr::str_detect(Sample,"^1_") ~ "1st",
      stringr::str_detect(Sample,"^2_") ~ "2nd",
      TRUE ~ "1st"
    ),
    XLabel_bc = paste0(Dose, "\n", Trial_label, " ", Tissue)
  ) %>%
  dplyr::filter(!is.na(Dose), Tissue!="Unknown")

# ── [수정 5] Figure B/C 그룹 레벨: Trial × Tissue 분리 ──────────────────────
# 존재하는 조합만 레벨로 사용
bc_levels_df <- all_isnv %>%
  dplyr::distinct(Dose, Trial_label, Tissue) %>%
  dplyr::mutate(
    Dose_num = as.numeric(Dose),
    Trial_ord = ifelse(Trial_label=="1st", 1, 2),
    Tissue_ord = ifelse(Tissue=="URT", 1, 2)
  ) %>%
  dplyr::arrange(Dose_num, Trial_ord, Tissue_ord) %>%
  dplyr::mutate(Group_label = paste0(Dose, "\n", Trial_label, " ", Tissue))

bc_group_levels <- bc_levels_df$Group_label

# Figure B data
b_data <- all_isnv %>%
  dplyr::distinct(Dose, Trial_label, Tissue, Gene, Position, Syn_Status) %>%
  dplyr::count(Dose, Trial_label, Tissue, Gene, Syn_Status, name="Unique_iSNV") %>%
  dplyr::mutate(
    Group      = factor(paste0(Dose,"\n",Trial_label," ",Tissue), levels=bc_group_levels),
    Tissue     = Tissue,   # 배경색용 유지
    Syn_Status = factor(Syn_Status, levels=c("S","NS")),
    Gene       = factor(Gene, levels=genes_of_interest)
  )

# Figure B 배경 rect 데이터 — fill을 aes 밖에서 직접 지정
b_bg <- b_data %>%
  dplyr::distinct(Group, Tissue) %>%
  dplyr::mutate(
    x_num      = as.integer(Group),
    fill_color = TISSUE_COLORS[Tissue]
  )

p_b <- ggplot(b_data, aes(x=Group, y=Unique_iSNV, fill=Syn_Status)) +
  # ★ Tissue 배경: fill을 aes() 밖으로 꺼내 막대 fill과 충돌 방지
  purrr::pmap(list(b_bg$x_num, b_bg$fill_color), function(xn, fc) {
    annotate("rect", xmin=xn-0.5, xmax=xn+0.5,
             ymin=-Inf, ymax=Inf, fill=fc, alpha=TISSUE_ALPHA)
  }) +
  geom_col(width=0.65, color=NA) +
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  scale_fill_manual(values=muttype_colors,
                    labels=c("NS"="Non-synonymous","S"="Synonymous"),
                    breaks=c("NS","S")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  facet_wrap(~Gene, ncol=3, scales="free_y") +
  labs(y="Number of unique iSNVs") +
  theme_unified(x_blank=FALSE) +
  theme(axis.text.x = element_text(size=8, angle=45, hjust=1, lineheight=0.85))

# Figure C data
sharing_raw <- all_isnv %>%
  dplyr::distinct(Dose, Trial_label, Tissue, Gene, Position, Sample) %>%
  dplyr::count(Dose, Trial_label, Tissue, Gene, Position, name="N_individuals") %>%
  dplyr::mutate(
    Sharing=dplyr::case_when(
      N_individuals==1 ~ "1 individual",
      N_individuals==2 ~ "2 individuals",
      N_individuals>=3 ~ "3+ individuals"),
    Sharing=factor(Sharing, levels=c("1 individual","2 individuals","3+ individuals"))
  )

c_data <- sharing_raw %>%
  dplyr::count(Dose, Trial_label, Tissue, Gene, Sharing, name="Unique_iSNV") %>%
  dplyr::mutate(
    Group  = factor(paste0(Dose,"\n",Trial_label," ",Tissue), levels=bc_group_levels),
    Tissue = Tissue,   # 배경색용 유지
    Gene   = factor(Gene, levels=genes_of_interest)
  )

# Figure C 배경 rect 데이터
c_bg <- c_data %>%
  dplyr::distinct(Group, Tissue) %>%
  dplyr::mutate(
    x_num      = as.integer(Group),
    fill_color = TISSUE_COLORS[Tissue]
  )

p_c <- ggplot(c_data, aes(x=Group, y=Unique_iSNV, fill=Sharing)) +
  # ★ Tissue 배경
  purrr::pmap(list(c_bg$x_num, c_bg$fill_color), function(xn, fc) {
    annotate("rect", xmin=xn-0.5, xmax=xn+0.5,
             ymin=-Inf, ymax=Inf, fill=fc, alpha=TISSUE_ALPHA)
  }) +
  geom_col(width=0.65, color=NA) +
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  scale_fill_manual(values=sharing_colors, name="Shared in") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  facet_wrap(~Gene, ncol=3, scales="free_y") +
  labs(y="Number of unique iSNVs") +
  theme_unified(x_blank=FALSE) +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size=8, angle=45, hjust=1, lineheight=0.85))

# ══════════════════════════════════════════════════════════════════════════════
# COMBINED (patchwork)
# ══════════════════════════════════════════════════════════════════════════════
cat("\nCombining with patchwork...\n")

combined <- p_a / (p_b | p_c) +
  patchwork::plot_layout(heights=c(1.5, 1.2))

fig_out <- paste0(out_pfx, ".png")
ggsave(fig_out, combined, width=18, height=14, dpi=300, bg="white")
cat("Saved:", basename(fig_out), "\n")

# 개별 저장
ggsave(paste0(out_pfx,"_A.png"), p_a, width=18, height=7,   dpi=300, bg="white")
ggsave(paste0(out_pfx,"_B.png"), p_b, width=10, height=6.5, dpi=300, bg="white")
ggsave(paste0(out_pfx,"_C.png"), p_c, width=10, height=6.5, dpi=300, bg="white")
cat("Individual figures also saved.\n")

cat("\n✓ Done!\n")
cat("  Strip size :", STRIP_SIZE, "pt  (italic 제거됨)\n")
cat("  Axis text  :", AXIS_TEXT,  "pt\n")
cat("  Legend text:", LEGEND_TEXT,"pt\n")
cat("  X-axis labels: Trial + Site 형식 (예: 1st URT, 2nd LRT)\n")