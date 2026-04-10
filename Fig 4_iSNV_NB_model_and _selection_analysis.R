################################################################################
# Two-Set dNdS Analysis
#
# 두 세트(A, B)를 각각 분석하여 동일한 figure에 나란히 표기
#
# Figure 1 (combined): IRR forest (Total/NS/S) + IRR ratio bar
#   - 상단: Set A / 하단: Set B
#   - facet by Gene, 각 세트별 IRR + ratio
#
# Figure 4 (combined): pN/pS 평균 ± SD (개체 점 없음, 오른쪽 패널 없음)
#   - Set A (왼쪽 열) / Set B (오른쪽 열)
#   - 배경색: 그룹별 pN/pS 기반
#
# Input: 각 세트당 2개 파일
#   [1] iSNV summary (7__iSNV_summary...Syn분류.xlsx)
#   [2] callable sites (10__최종정리본...xlsx)
################################################################################

required_packages <- c("readxl","dplyr","purrr","ggplot2","stringr",
                        "MASS","openxlsx","tidyr","scales","patchwork")
for (pkg in required_packages) {
  if (!require(pkg, character.only=TRUE, quietly=TRUE)) {
    install.packages(pkg, dependencies=TRUE, quiet=TRUE)
    library(pkg, character.only=TRUE)
  }
}

# ── 유전 코드 + possible sites ────────────────────────────────────────────────
GENETIC_CODE <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
  AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
  TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
  AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G"
)
BASES <- c("A","T","G","C")

calc_possible_sites <- function(codon) {
  codon <- toupper(trimws(codon))
  if (nchar(codon)!=3||grepl("[^ATGC]",codon)) return(list(S=NA,N=NA))
  ref_aa <- GENETIC_CODE[codon]; if (is.na(ref_aa)) return(list(S=NA,N=NA))
  sT=nT=0
  for (pos in 1:3) {
    s=n=0
    for (alt in setdiff(BASES,substr(codon,pos,pos))) {
      mc=codon; substr(mc,pos,pos)=alt; ma=GENETIC_CODE[mc]
      if (is.na(ma)) next
      if (ma=="*"||ma!=ref_aa) n=n+1 else s=s+1
    }
    sT=sT+s/3; nT=nT+n/3
  }
  list(S=sT,N=nT)
}

genes_of_interest <- c("Spike","ORF1a","ORF1b","ORF3a","ORF8","ORF10")
gene_ref_lengths  <- c(Spike=3822,ORF1a=13218,ORF1b=8085,
                        ORF3a=828,ORF8=366,ORF10=117)

# ── 전역 색상/모양/레이블 정의 ─────────────────────────────────────────────────
SET_COLORS  <- c("A"="#c0392b", "B"="#2980b9")   # A=빨강, B=파랑
SET_SHAPES  <- c("A"=16,        "B"=17)           # A=원, B=삼각
SET_LABELS  <- c("A"="Filtered iSNVs (excluding variants present in positive control)",
               "B"="Unfiltered iSNVs")
SET_BASE_SIZE <- 13

# ── 공통 헬퍼 함수 ────────────────────────────────────────────────────────────
make_colnames <- function(infile, sheet_name) {
  raw_top  <- readxl::read_excel(infile,sheet=sheet_name,col_names=FALSE,n_max=2)
  gene_row <- as.character(unlist(raw_top[1,]))
  for (i in 2:length(gene_row))
    if (is.na(gene_row[i])||gene_row[i]=="NA") gene_row[i] <- gene_row[i-1]
  sub_row   <- as.character(unlist(raw_top[2,]))
  col_names <- character(length(sub_row))
  for (j in seq_along(col_names)) {
    g=gene_row[j]; s=sub_row[j]
    if (is.na(s)||s=="NA")                      {col_names[j]=paste0("Col",j);next}
    if (grepl("^Sample$",s,ignore.case=TRUE))    {col_names[j]="Sample";next}
    if (grepl("^Group$",s,ignore.case=TRUE))     {col_names[j]="Group";next}
    if (grepl("Total reads",s,ignore.case=TRUE)) {col_names[j]="Total_reads";next}
    sc <- s %>%
      stringr::str_replace_all("\\n"," ") %>%
      stringr::str_replace_all("iSNVs count","iSNVs") %>%
      stringr::str_replace_all("callable sites","callable") %>%
      stringr::str_replace_all("Coverage\\s*","coverage") %>%
      stringr::str_replace_all("iSNVs per\\s+callable sites","rate") %>%
      stringr::str_replace_all("iSNVs per callable sites","rate") %>%
      stringr::str_trim()
    col_names[j] <- paste0(g,"_",sc)
  }
  col_names
}

parse_meta <- function(infile) {
  cnames <- make_colnames(infile,"Total")
  rd <- readxl::read_excel(infile,sheet="Total",col_names=FALSE,skip=2)
  names(rd) <- cnames
  rd %>%
    dplyr::mutate(Dose=as.numeric(as.character(Group)),
                  Dose=ifelse(is.na(Dose),0,Dose)) %>%
    dplyr::filter(!is.na(Sample)) %>%
    dplyr::select(Sample,Dose)
}

parse_sheet <- function(infile, sheet_name, meta) {
  cnames <- make_colnames(infile,sheet_name)
  rd <- readxl::read_excel(infile,sheet=sheet_name,col_names=FALSE,skip=2)
  names(rd) <- cnames
  rd <- rd %>%
    dplyr::mutate(Dose=as.numeric(as.character(Group)),
                  Dose=ifelse(is.na(Dose),0,Dose)) %>%
    dplyr::filter(!is.na(as.numeric(Total_reads)))
  if (!"Sample" %in% names(rd)) {
    n <- min(nrow(rd),nrow(meta))
    rd <- rd[seq_len(n),]; rd$Sample <- meta$Sample[seq_len(n)]
  }
  genes_avail <- intersect(genes_of_interest,
    names(gene_ref_lengths)[sapply(names(gene_ref_lengths), function(g)
      paste0(g,"_iSNVs") %in% names(rd) & paste0(g,"_callable") %in% names(rd))])
  purrr::map_dfr(genes_avail, function(g) {
    rl <- gene_ref_lengths[[g]]
    rd %>% dplyr::transmute(
      Sample=as.character(Sample), Dose=Dose, MutType=sheet_name, Gene=g,
      Count=as.integer(.data[[paste0(g,"_iSNVs")]]),
      Callable=as.numeric(.data[[paste0(g,"_callable")]]),
      Callable_norm=Callable/rl
    ) %>% dplyr::filter(is.finite(Count),is.finite(Callable),Callable>0)
  })
}

# ── 세트별 분석 함수 ──────────────────────────────────────────────────────────
analyze_set <- function(isnv_file, callable_file, set_name) {
  cat(sprintf("\n=== Analyzing Set %s ===\n", set_name))

  # possible sites
  codon_df <- purrr::map_dfr(genes_of_interest, function(g) {
    raw <- readxl::read_excel(isnv_file,sheet=g,col_names=FALSE)
    h1  <- as.character(unlist(raw[1,]))
    rc  <- which(h1=="Reference")[1]; cc <- which(h1=="Ref_codon")[1]
    res <- list()
    for (i in 3:nrow(raw)) {
      rv  <- as.character(unlist(raw[i,]))
      pos <- suppressWarnings(as.numeric(rv[rc]))
      cv  <- rv[cc]
      if (is.na(pos)||is.na(cv)||nchar(trimws(cv))!=3||grepl("^=",cv)) next
      res[[length(res)+1]] <- data.frame(Gene=g,Position=pos,
                                          Codon=toupper(trimws(cv)),
                                          stringsAsFactors=FALSE)
    }
    if (!length(res)) return(NULL)
    dplyr::bind_rows(res)
  })

  gene_possible <- codon_df %>%
    dplyr::distinct(Gene,Position,Codon) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sites=list(calc_possible_sites(Codon)),
                  poss_S=sites$S, poss_N=sites$N) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(poss_S)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Total_poss_S=sum(poss_S),
                     Total_poss_N=sum(poss_N), .groups="drop")

  meta <- parse_meta(callable_file)

  # NS/S (Dose>0, for NB)
  count_data <- purrr::map_dfr(c("NS","S","Total"),
    ~ parse_sheet(callable_file,.x,meta)) %>%
    dplyr::left_join(
      gene_possible %>% dplyr::mutate(Total_poss_all=Total_poss_N+Total_poss_S),
      by="Gene") %>%
    dplyr::mutate(
      Possible_sites=dplyr::case_when(
        MutType=="NS"~Total_poss_N, MutType=="S"~Total_poss_S,
        MutType=="Total"~Total_poss_all),
      log_offset=log(Possible_sites*Callable_norm)
    ) %>%
    dplyr::filter(Dose>0, is.finite(log_offset))

  # NB regression
  run_nb <- function(df,g,mt) {
    dat <- df %>% dplyr::filter(Gene==g,MutType==mt,is.finite(Count),Count>=0)
    if (nrow(dat)<5) return(NULL)
    tryCatch({
      model <- MASS::glm.nb(Count~log10(Dose)+offset(log_offset),data=dat)
      cs <- summary(model)$coefficients
      cv <- cs["log10(Dose)","Estimate"]; se <- cs["log10(Dose)","Std. Error"]
      pv <- cs["log10(Dose)","Pr(>|z|)"]; irr <- exp(cv)
      data.frame(Gene=g,MutType=mt,IRR=irr,
        IRR_CI_lower=max(exp(cv-1.96*se),0.001),
        IRR_CI_upper=min(exp(cv+1.96*se),10000),
        P_value=pv,
        Significant=ifelse(pv<0.001,"***",ifelse(pv<0.01,"**",
                    ifelse(pv<0.05,"*","ns"))),
        stringsAsFactors=FALSE)
    },error=function(e) NULL)
  }

  m1 <- purrr::map_dfr(genes_of_interest, function(g)
    dplyr::bind_rows(run_nb(count_data,g,"NS"),
                     run_nb(count_data,g,"S"),
                     run_nb(count_data,g,"Total")))

  # Interaction model
  m2 <- purrr::map_dfr(genes_of_interest, function(g) {
    dat <- count_data %>% dplyr::filter(Gene==g,MutType %in% c("NS","S")) %>%
      dplyr::mutate(MutType=factor(MutType,levels=c("NS","S")))
    if (nrow(dat)<10) return(NULL)
    tryCatch({
      model <- MASS::glm.nb(Count~log10(Dose)*MutType+offset(log_offset),data=dat)
      cs <- summary(model)$coefficients
      cd <- cs["log10(Dose)","Estimate"]; sd_d <- cs["log10(Dose)","Std. Error"]
      pd <- cs["log10(Dose)","Pr(>|z|)"]
      ir <- "log10(Dose):MutTypeS"
      if (ir %in% rownames(cs)) {
        ci=cs[ir,"Estimate"]; si=cs[ir,"Std. Error"]; pi=cs[ir,"Pr(>|z|)"]
      } else { ci=NA; si=NA; pi=NA }
      ratio <- exp(ci)
      sig   <- ifelse(is.na(pi),"NA",ifelse(pi<0.001,"***",
               ifelse(pi<0.01,"**",ifelse(pi<0.05,"*","ns"))))
      data.frame(Gene=g,
        IRR_NS=exp(cd), IRR_NS_lo=max(exp(cd-1.96*sd_d),0.001),
        IRR_NS_hi=min(exp(cd+1.96*sd_d),10000), P_NS=pd,
        Sig_NS=ifelse(pd<0.001,"***",ifelse(pd<0.01,"**",ifelse(pd<0.05,"*","ns"))),
        IRR_S=exp(cd+ci),
        IRR_S_lo=max(exp(cd+ci-1.96*sqrt(sd_d^2+si^2)),0.001),
        IRR_S_hi=min(exp(cd+ci+1.96*sqrt(sd_d^2+si^2)),10000),
        Sig_S=ifelse((pd)<0.001,"***",ifelse(pd<0.01,"**",ifelse(pd<0.05,"*","ns"))),
        IRR_ratio=ratio, P_interaction=pi, Sig_interaction=sig,
        stringsAsFactors=FALSE)
    },error=function(e) NULL)
  })

  # pN/pS (Dose=0 포함)
  count_all <- purrr::map_dfr(c("NS","S"),
    ~ parse_sheet(callable_file,.x,meta)) %>%
    dplyr::left_join(gene_possible,by="Gene")

  indiv_pnps <- count_all %>%
    dplyr::filter(MutType %in% c("NS","S")) %>%
    dplyr::select(Sample,Dose,Gene,MutType,Count,Total_poss_N,Total_poss_S) %>%
    tidyr::pivot_wider(names_from=MutType,
                       values_from=c(Count,Total_poss_N,Total_poss_S),
                       names_sep="_") %>%
    dplyr::mutate(
      pN=Count_NS/Total_poss_N_NS, pS=Count_S/Total_poss_S_S,
      pNpS=ifelse(!is.na(pS)&pS>0&!is.na(pN),pN/pS,NA_real_),
      Dose_label=ifelse(Dose==0,"Pos",as.character(Dose))
    ) %>%
    dplyr::filter(!is.na(pNpS),is.finite(pNpS))

  # Dose factor levels (Pos 있을 때만)
  has_pos <- any(indiv_pnps$Dose==0)
  dose_lvs <- if (has_pos) c("Pos","0.01","0.1","1") else c("0.01","0.1","1")
  indiv_pnps <- indiv_pnps %>%
    dplyr::mutate(Dose_f=factor(Dose_label,levels=dose_lvs)) %>%
    dplyr::filter(!is.na(Dose_f))

  list(m1=m1, m2=m2, indiv_pnps=indiv_pnps,
       gene_possible=gene_possible, dose_lvs=dose_lvs,
       set_name=set_name)
}

# ════════════════════════════════════════════════════════════════
# 파일 선택
# ════════════════════════════════════════════════════════════════
cat("=== Set A ===\n")
cat("[A-1] Select iSNV summary (7__iSNV_summary...Syn분류.xlsx)\n")
isnv_A <- file.choose()
cat("[A-2] Select callable sites (10__최종정리본...xlsx)\n")
callable_A <- file.choose()

cat("\n=== Set B ===\n")
cat("[B-1] Select iSNV summary (7__iSNV_summary...Syn분류.xlsx)\n")
isnv_B <- file.choose()
cat("[B-2] Select callable sites (10__최종정리본...xlsx)\n")
callable_B <- file.choose()

out_dir <- dirname(callable_A)
out_pfx <- file.path(out_dir, "TwoSet_Analysis")

# ── 분석 실행 ─────────────────────────────────────────────────────────────────
res_A <- analyze_set(isnv_A, callable_A, "A")
res_B <- analyze_set(isnv_B, callable_B, "B")

# ════════════════════════════════════════════════════════════════
# FIGURE 1: IRR Forest + IRR Ratio
# Set A (상단) / Set B (하단) — patchwork
# ════════════════════════════════════════════════════════════════
cat("\nGenerating Figure 1 (IRR Forest + Ratio, two sets)...\n")

type_colors <- c("NS"="#c0392b","S"="#2980b9","Total"="#7f8c8d")
theme_forest <- theme_bw(base_size=11) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=11,face="bold"),
        axis.text.x=element_text(size=10),
        strip.text=element_text(face="bold",size=11),
        strip.background=element_rect(fill="gray90"),
        legend.position="bottom", legend.title=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title=element_blank(), plot.subtitle=element_blank(),
        plot.background=element_rect(fill="white",color=NA),
        panel.spacing=unit(0.8,"lines"), plot.margin=margin(8,12,8,12))


# ── Fig1: A/B를 하나의 패널에 ──────────────────────────────────────────────────
make_irr_combined <- function(res_a, res_b) {

  make_df <- function(res, set_nm) {
    m1_tot <- res$m1 %>% dplyr::filter(MutType=="Total") %>%
      dplyr::transmute(Gene,MutType="Total",IRR,
        CI_lo=IRR_CI_lower,CI_hi=IRR_CI_upper,Sig=Significant,Set=set_nm)
    m2_ns <- res$m2 %>%
      dplyr::transmute(Gene,MutType="NS",IRR=IRR_NS,
        CI_lo=IRR_NS_lo,CI_hi=IRR_NS_hi,Sig=Sig_NS,Set=set_nm)
    m2_s <- res$m2 %>%
      dplyr::transmute(Gene,MutType="S",IRR=IRR_S,
        CI_lo=IRR_S_lo,CI_hi=IRR_S_hi,Sig=Sig_S,Set=set_nm)
    dplyr::bind_rows(m2_s,m2_ns,m1_tot)
  }

  df <- dplyr::bind_rows(make_df(res_a,"A"),make_df(res_b,"B")) %>%
    dplyr::mutate(
      MutType = factor(MutType,levels=c("S","NS","Total")),
      Gene    = factor(Gene,levels=genes_of_interest),
      Set     = factor(Set,levels=c("A","B")),
      y_num   = as.numeric(MutType),
      y_off   = y_num + ifelse(Set=="A",0.2,-0.2)
    )

  # IRR ratio 텍스트
  ratio_df <- dplyr::bind_rows(
    res_a$m2 %>% dplyr::mutate(Set="A"),
    res_b$m2 %>% dplyr::mutate(Set="B")
  ) %>%
    dplyr::mutate(
      IRR_ratio=ifelse(abs(IRR_ratio)>1000|is.infinite(IRR_ratio),
                       NA_real_,IRR_ratio),
      Sig_interaction=ifelse(is.na(IRR_ratio),"—",Sig_interaction),
      Gene=factor(Gene,levels=genes_of_interest),
      Set =factor(Set,levels=c("A","B")),
      symbol=ifelse(Set=="A","●","▲"),
      label=sprintf("%s IRR(S/NS)=%.2f (%s)",
                    symbol,IRR_ratio,Sig_interaction),
      y_pos=ifelse(Set=="A",3.55,3.25)
    ) %>%
    dplyr::filter(!is.na(IRR_ratio))

  all_ci  <- c(df$CI_lo,df$CI_hi)
  xmn <- max(min(all_ci,na.rm=TRUE)*0.7,0.3)
  xmx <- min(max(all_ci,na.rm=TRUE)*1.3,100)
  brk <- c(0.5,1,2,5,10)[c(0.5,1,2,5,10)>=xmn & c(0.5,1,2,5,10)<=xmx]

  ggplot(df,aes(x=IRR,y=y_off,color=Set,shape=Set)) +
    geom_vline(xintercept=1,linetype="dashed",color="gray50",linewidth=0.7) +
    geom_errorbarh(aes(xmin=pmax(CI_lo,xmn),xmax=pmin(CI_hi,xmx)),
                   height=0.14,linewidth=1.3) +
    geom_point(size=4) +
    geom_text(aes(label=sprintf("%.2f %s",IRR,Sig)),
              hjust=0.5,vjust=-1.2,size=3.5,fontface="bold",
              show.legend=FALSE) +
    # ratio 텍스트
    geom_text(data=ratio_df,
              aes(x=xmx^0.95,y=y_pos,
                  label=label,color=Set),
              hjust=1,size=3.5,fontface="bold",
              inherit.aes=FALSE,show.legend=FALSE) +
    scale_color_manual(values=SET_COLORS,
                       labels=SET_LABELS,name=NULL) +
    scale_shape_manual(values=SET_SHAPES,
                       labels=SET_LABELS,name=NULL) +
    scale_x_log10(limits=c(xmn,xmx),breaks=brk,
                  labels=scales::label_number(accuracy=0.1),
                  name="Incidence Rate Ratio (IRR, log scale)") +
    scale_y_continuous(breaks=c(1,2,3),
                       labels=c("S","NS","Total"),
                       limits=c(0.5,3.8)) +
    facet_wrap(~Gene,ncol=3,scales="free_x") +
    theme_bw(base_size=SET_BASE_SIZE) +
    theme(
      axis.title.y      = element_blank(),
      axis.text.y       = element_text(size=SET_BASE_SIZE,face="bold"),
      axis.text.x       = element_text(size=SET_BASE_SIZE-2),
      strip.text        = element_text(face="bold",size=SET_BASE_SIZE),
      strip.background  = element_rect(fill="gray90"),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(size=SET_BASE_SIZE),
      legend.key.size   = unit(0.8,"cm"),
      panel.grid.minor  = element_blank(),
      panel.grid.major.y= element_line(color="gray90",linewidth=0.3),
      plot.background   = element_rect(fill="white",color=NA),
      plot.title=element_blank(),plot.subtitle=element_blank(),
      panel.spacing     = unit(1,"lines"),
      plot.margin       = margin(10,15,10,15)
    ) +
    guides(color=guide_legend(override.aes=list(size=5)),
           shape=guide_legend(override.aes=list(size=5)))
}


p_fig1 <- make_irr_combined(res_A, res_B)

fig1_out <- paste0(out_pfx,"_Fig1_IRR_TwoSets.png")
ggsave(fig1_out, p_fig1, width=16, height=9, dpi=300, bg="white")
cat("Figure 1 saved:", basename(fig1_out),"\n")

# ════════════════════════════════════════════════════════════════
# FIGURE 4: pN/pS 평균 ± SD (개체 점 없음)
# Set A (왼쪽) / Set B (오른쪽)
# ════════════════════════════════════════════════════════════════
cat("\nGenerating Figure 4 (pN/pS mean±SD, two sets)...\n")

pnps_to_bg <- function(val, alpha=0.40) {
  if (is.na(val)||!is.finite(val)) return(grDevices::adjustcolor("gray85",alpha))
  if (val<1) {
    intensity <- max(0,min(1,1-val))
    col <- grDevices::colorRamp(c("gray90","steelblue"))(intensity)
  } else {
    intensity <- max(0,min(1,(val-1)/1.5))
    col <- grDevices::colorRamp(c("gray90","#c0392b"))(intensity)
  }
  grDevices::adjustcolor(grDevices::rgb(col[1]/255,col[2]/255,col[3]/255),alpha)
}

make_pnps_plot <- function(res, set_label) {
  ip   <- res$indiv_pnps
  dlvs <- res$dose_lvs

  # 그룹별 평균+SD
  mean_sd <- ip %>%
    dplyr::group_by(Gene, Dose_f) %>%
    dplyr::summarise(
      x_pos     = as.numeric(Dose_f[1]),
      mean_pNpS = mean(pNpS,na.rm=TRUE),
      sd_pNpS   = sd(pNpS,na.rm=TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(Gene=factor(Gene,levels=genes_of_interest),
                  x_pos=as.numeric(Dose_f))

  # 배경색 (gene × dose group별 pN/pS 기반)
  bg_df <- mean_sd %>%
    dplyr::mutate(
      xmin   = x_pos-0.5,
      xmax   = x_pos+0.5,
      bg_col = mapply(pnps_to_bg, mean_pNpS)
    )

  # KW test
  kw_df <- purrr::map_dfr(genes_of_interest, function(g) {
    dat <- ip %>% dplyr::filter(Gene==g,is.finite(pNpS))
    if (nrow(dat)<4||length(unique(dat$Dose_f))<2) return(NULL)
    kw <- tryCatch(kruskal.test(pNpS~Dose_f,data=dat),error=function(e) NULL)
    if (is.null(kw)) return(NULL)
    kp <- kw$p.value
    data.frame(Gene=g,
      KW_sig=ifelse(kp<0.001,"***",ifelse(kp<0.01,"**",
             ifelse(kp<0.05,"*","ns"))),
      x_mid=(length(dlvs)+1)/2, stringsAsFactors=FALSE)
  }) %>% dplyr::mutate(Gene=factor(Gene,levels=genes_of_interest),
                        label=paste0("KW: ",KW_sig))

  pnps_max <- max(c(mean_sd$mean_pNpS+mean_sd$sd_pNpS,1),na.rm=TRUE)
  color_lim <- c(0,max(mean_sd$mean_pNpS,na.rm=TRUE)*1.05)

  ggplot() +
    geom_rect(data=bg_df,
              aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=bg_col),
              inherit.aes=FALSE) +
    scale_fill_identity(guide="none") +
    geom_hline(yintercept=1,linetype="dashed",color="gray40",linewidth=0.6) +
    geom_errorbar(data=mean_sd,
                  aes(x=x_pos,
                      ymin=pmax(mean_pNpS-sd_pNpS,0),
                      ymax=mean_pNpS+sd_pNpS),
                  width=0.22,linewidth=0.9,color="gray20",
                  inherit.aes=FALSE) +
    geom_point(data=mean_sd,
               aes(x=x_pos,y=mean_pNpS,color=mean_pNpS),
               size=4.5,shape=23,fill="white",stroke=1.4,
               inherit.aes=FALSE) +
    geom_text(data=mean_sd,
              aes(x=x_pos,y=mean_pNpS,label=sprintf("%.2f",mean_pNpS)),
              vjust=-1.3,size=2.6,fontface="bold",color="gray20",
              show.legend=FALSE,inherit.aes=FALSE) +
    geom_text(data=kw_df,
              aes(x=x_mid,y=Inf,label=label),
              inherit.aes=FALSE,vjust=1.4,size=2.8,
              color="gray20",fontface="bold") +
    scale_color_gradient2(low="steelblue",mid="gray75",high="#c0392b",
                          midpoint=1,name="pN/pS",limits=color_lim,
                          guide=guide_colorbar()) +
    scale_x_continuous(breaks=seq_along(dlvs),labels=dlvs,
                       expand=expansion(add=0.7)) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.22))) +
    facet_wrap(~Gene,ncol=3,scales="free_y") +
    labs(x="Dose group",
         y=sprintf("pN/pS — Set %s",set_label)) +
    theme_classic(base_size=11) +
    theme(
      axis.line          = element_line(color="black",linewidth=0.5),
      axis.text.x        = element_text(size=10,face="bold"),
      axis.text.y        = element_text(size=9),
      axis.title         = element_text(size=11,face="bold"),
      strip.text         = element_text(face="bold.italic",size=11),
      strip.background   = element_rect(fill="gray95",color=NA),
      legend.position    = "right",
      legend.title       = element_text(size=9,face="bold"),
      legend.text        = element_text(size=8),
      panel.grid.major.y = element_line(color="gray90",linewidth=0.3),
      plot.background    = element_rect(fill="white",color=NA),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.margin        = margin(10,15,10,15),
      panel.spacing      = unit(0.8,"lines")
    )
}

# ── Fig4: A/B 한 패널 — 모양으로 구분 ────────────────────────────────────────
make_pnps_combined <- function(res_a, res_b) {

  prep <- function(res, set_nm) {
    ip   <- res$indiv_pnps
    dlvs <- res$dose_lvs
    ms <- ip %>%
      dplyr::group_by(Gene, Dose_f) %>%
      dplyr::summarise(
        x_pos     = as.numeric(Dose_f[1]),
        mean_pNpS = mean(pNpS,na.rm=TRUE),
        sd_pNpS   = sd(pNpS,na.rm=TRUE),
        .groups="drop"
      ) %>%
      dplyr::mutate(
        Gene  = factor(Gene,levels=genes_of_interest),
        x_pos = as.numeric(Dose_f),
        Set   = set_nm,
        # x 살짝 offset
        x_jit = x_pos + ifelse(set_nm=="A",-0.15,0.15)
      )
    list(ms=ms, dlvs=dlvs)
  }

  pa <- prep(res_a,"A"); pb <- prep(res_b,"B")
  ms_all <- dplyr::bind_rows(pa$ms, pb$ms) %>%
    dplyr::mutate(Set=factor(Set,levels=c("A","B")))

  # 배경색: Set A의 pN/pS 기준 (A가 없으면 B)
  bg_base <- pa$ms %>%
    dplyr::mutate(xmin=x_pos-0.5, xmax=x_pos+0.5,
                  bg_col=mapply(pnps_to_bg,mean_pNpS))

  # KW test (set별)
  kw_list <- list()
  for (set_nm in c("A","B")) {
    res <- if (set_nm=="A") res_a else res_b
    ip  <- res$indiv_pnps
    for (g in genes_of_interest) {
      dat <- ip %>% dplyr::filter(Gene==g,is.finite(pNpS))
      if (nrow(dat)<4||length(unique(dat$Dose_f))<2) next
      kw <- tryCatch(kruskal.test(pNpS~Dose_f,data=dat),error=function(e) NULL)
      if (is.null(kw)) next
      kp <- kw$p.value
      kw_list[[length(kw_list)+1]] <- data.frame(
        Gene=g, Set=set_nm,
        KW_sig=ifelse(kp<0.001,"***",ifelse(kp<0.01,"**",
               ifelse(kp<0.05,"*","ns"))),
        stringsAsFactors=FALSE)
    }
  }
  kw_df <- dplyr::bind_rows(kw_list) %>%
    dplyr::mutate(Gene=factor(Gene,levels=genes_of_interest),
                  Set=factor(Set,levels=c("A","B")),
                  label=paste0(Set,": ",KW_sig),
                  # A는 위, B는 약간 아래
                  y_inf=ifelse(Set=="A",Inf,Inf),
                  vjust_v=ifelse(Set=="A",1.4,2.8),
                  x_mid=length(pa$dlvs)/2+0.5)

  dlvs <- pa$dlvs
  color_lim <- c(0,max(ms_all$mean_pNpS,na.rm=TRUE)*1.05)

  # Set 색상: A=진한 파랑/빨강, B=연한 파랑/빨강
  set_colors  <- c("A"="#c0392b","B"="#2980b9")
  set_shapes  <- c("A"=23,"B"=22)  # 다이아몬드 vs 사각형
  set_fills   <- c("A"="white","B"="white")

  ggplot() +
    geom_rect(data=bg_base,
              aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf,fill=bg_col),
              inherit.aes=FALSE) +
    scale_fill_identity(guide="none") +
    geom_hline(yintercept=1,linetype="dashed",color="gray40",linewidth=0.6) +
    # 에러바
    geom_errorbar(data=ms_all,
                  aes(x=x_jit,
                      ymin=pmax(mean_pNpS-sd_pNpS,0),
                      ymax=mean_pNpS+sd_pNpS,
                      color=Set),
                  width=0.18,linewidth=0.9,inherit.aes=FALSE) +
    # 평균 점
    geom_point(data=ms_all,
               aes(x=x_jit,y=mean_pNpS,
                   color=Set,shape=Set,fill=Set),
               size=4,stroke=1.3,inherit.aes=FALSE) +
    # 수치
    geom_text(data=ms_all,
              aes(x=x_jit,y=mean_pNpS,
                  label=sprintf("%.2f",mean_pNpS),color=Set),
              vjust=-1.3,size=2.4,fontface="bold",
              show.legend=FALSE,inherit.aes=FALSE) +
    # KW 레이블 (A, B 각각)
    geom_text(data=kw_df,
              aes(x=x_mid,y=Inf,label=label,color=Set),
              vjust=kw_df$vjust_v,size=2.6,fontface="bold",
              inherit.aes=FALSE,show.legend=FALSE) +
    scale_color_manual(values=set_colors,name="Set") +
    scale_shape_manual(values=set_shapes,name="Set") +
    scale_fill_manual(values=set_fills, guide="none") +
    scale_x_continuous(breaks=seq_along(dlvs),labels=dlvs,
                       expand=expansion(add=0.7)) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.25))) +
    facet_wrap(~Gene,ncol=3,scales="free_y") +
    labs(x="Dose group",y="pN/pS (mean ± SD)") +
    theme_classic(base_size=11) +
    theme(
      axis.line          = element_line(color="black",linewidth=0.5),
      axis.text.x        = element_text(size=10,face="bold"),
      axis.text.y        = element_text(size=9),
      axis.title         = element_text(size=11,face="bold"),
      strip.text         = element_text(face="bold.italic",size=11),
      strip.background   = element_rect(fill="gray95",color=NA),
      legend.position    = "bottom",
      legend.title       = element_text(size=10,face="bold"),
      legend.text        = element_text(size=10),
      panel.grid.major.y = element_line(color="gray90",linewidth=0.3),
      plot.background    = element_rect(fill="white",color=NA),
      plot.title=element_blank(), plot.subtitle=element_blank(),
      plot.margin        = margin(10,15,10,15),
      panel.spacing      = unit(0.8,"lines")
    ) +
    guides(color=guide_legend(override.aes=list(shape=c(23,22),size=4)),
           shape="none")
}

p_fig4 <- make_pnps_combined(res_A,res_B)
fig4_out <- paste0(out_pfx,"_Fig4_pNpS_TwoSets.png")
ggsave(fig4_out, p_fig4, width=14, height=8, dpi=300, bg="white")
cat("Figure 4 saved:", basename(fig4_out),"\n")

# ── Excel 저장 ────────────────────────────────────────────────────────────────
out_xlsx <- paste0(out_pfx,"_RESULTS.xlsx")
wb <- openxlsx::createWorkbook()
for (set_nm in c("A","B")) {
  res <- if (set_nm=="A") res_A else res_B
  openxlsx::addWorksheet(wb,paste0("Set",set_nm,"_IRR"))
  openxlsx::writeData(wb,paste0("Set",set_nm,"_IRR"),res$m1)
  openxlsx::addWorksheet(wb,paste0("Set",set_nm,"_Interaction"))
  openxlsx::writeData(wb,paste0("Set",set_nm,"_Interaction"),res$m2)
  openxlsx::addWorksheet(wb,paste0("Set",set_nm,"_pNpS"))
  pnps_out <- res$indiv_pnps %>%
    dplyr::select(Sample,Dose,Gene,Count_NS,Count_S,pN,pS,pNpS)
  openxlsx::writeData(wb,paste0("Set",set_nm,"_pNpS"),pnps_out)
}
openxlsx::saveWorkbook(wb,out_xlsx,overwrite=TRUE)
cat("Excel saved:", basename(out_xlsx),"\n")

cat("\n✓ Done!\n")
cat(sprintf("  Fig1 IRR (two sets) : %s\n",basename(fig1_out)))
cat(sprintf("  Fig4 pN/pS (two sets): %s\n",basename(fig4_out)))
cat(sprintf("  Excel               : %s\n",basename(out_xlsx)))# ── Fig4: A/B 한 패널 — 흰 배경, Fig1 동일 색/모양 ───────────────────────────
make_pnps_combined <- function(res_a, res_b) {

  prep <- function(res, set_nm) {
    ip <- res$indiv_pnps
    ms <- ip %>%
      dplyr::group_by(Gene,Dose_f) %>%
      dplyr::summarise(
        x_pos     = as.numeric(Dose_f[1]),
        mean_pNpS = mean(pNpS,na.rm=TRUE),
        sd_pNpS   = sd(pNpS,na.rm=TRUE),
        n         = dplyr::n(),
        .groups   = "drop"
      ) %>%
      dplyr::mutate(
        Gene  = factor(Gene,levels=genes_of_interest),
        x_pos = as.numeric(Dose_f),
        x_jit = x_pos + ifelse(set_nm=="A",-0.18,0.18),
        Set   = set_nm
      )
    ms
  }

  ms_all <- dplyr::bind_rows(prep(res_a,"A"),prep(res_b,"B")) %>%
    dplyr::mutate(Set=factor(Set,levels=c("A","B")))

  # KW test per set
  kw_df <- purrr::map_dfr(c("A","B"), function(set_nm) {
    res <- if (set_nm=="A") res_a else res_b
    ip  <- res$indiv_pnps
    purrr::map_dfr(genes_of_interest, function(g) {
      dat <- ip %>% dplyr::filter(Gene==g,is.finite(pNpS))
      if (nrow(dat)<4||length(unique(dat$Dose_f))<2) return(NULL)
      kw  <- tryCatch(kruskal.test(pNpS~Dose_f,data=dat),error=function(e)NULL)
      if (is.null(kw)) return(NULL)
      kp  <- kw$p.value
      data.frame(Gene=g,Set=set_nm,
        KW_sig=ifelse(kp<0.001,"***",ifelse(kp<0.01,"**",
               ifelse(kp<0.05,"*","ns"))),
        stringsAsFactors=FALSE)
    })
  }) %>%
    dplyr::mutate(
      Gene  = factor(Gene,levels=genes_of_interest),
      Set   = factor(Set,levels=c("A","B")),
      symbol = ifelse(Set=="A","●","▲"),
      label  = paste0(symbol," ",KW_sig),
      dlvs_n= dplyr::if_else(Set=="A",
                length(res_a$dose_lvs),length(res_b$dose_lvs)),
      x_mid = (dlvs_n+1)/2,
      vjust_v=ifelse(Set=="A",1.5,3.2)
    )

  # A/B 합집합 dose levels (Pos 포함 여부 자동)
  all_lvs_ordered <- c("Pos","0.01","0.1","1")
  dlvs <- all_lvs_ordered[all_lvs_ordered %in%
            unique(c(res_a$dose_lvs, res_b$dose_lvs))]

  # ms_all x_pos를 합집합 기준으로 재매핑
  ms_all <- ms_all %>%
    dplyr::mutate(
      Dose_label = as.character(Dose_f),
      x_pos_new  = match(Dose_label, dlvs),
      x_jit      = x_pos_new + ifelse(Set=="A",-0.18,0.18)
    ) %>%
    dplyr::filter(!is.na(x_pos_new))

  # KW 레이블: A/B 나란히 같은 x 위치, 한 줄로
  kw_wide <- kw_df %>%
    dplyr::select(Gene,Set,symbol,KW_sig) %>%
    tidyr::pivot_wider(names_from=Set,
                       values_from=c(symbol,KW_sig),
                       names_sep="_") %>%
    dplyr::mutate(
      Gene  = factor(Gene,levels=genes_of_interest),
      label = paste0(symbol_A," ",KW_sig_A,"  ",symbol_B," ",KW_sig_B),
      x_mid = length(dlvs)/2 + 0.5
    )

  ggplot() +
    geom_hline(yintercept=1,linetype="dashed",color="gray40",linewidth=0.7) +
    # 에러바
    geom_errorbar(data=ms_all,
                  aes(x=x_jit,
                      ymin=pmax(mean_pNpS-sd_pNpS,0),
                      ymax=mean_pNpS+sd_pNpS,
                      color=Set),
                  width=0.2,linewidth=1.0,inherit.aes=FALSE) +
    # 평균 점
    geom_point(data=ms_all,
               aes(x=x_jit,y=mean_pNpS,color=Set,shape=Set),
               size=4.5,stroke=1.4,fill="white",inherit.aes=FALSE) +
    # KW: A/B 나란히 한 줄
    geom_text(data=kw_wide,
              aes(x=x_mid,y=Inf,label=label),
              vjust=1.5,size=4,fontface="bold",color="gray20",
              inherit.aes=FALSE,show.legend=FALSE) +
    scale_color_manual(values=SET_COLORS,labels=SET_LABELS,name=NULL) +
    scale_shape_manual(values=SET_SHAPES,labels=SET_LABELS,name=NULL) +
    scale_fill_manual(values=c("A"="white","B"="white"),guide="none") +
    scale_x_continuous(breaks=seq_along(dlvs),labels=dlvs,
                       expand=expansion(add=0.7)) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.28))) +
    facet_wrap(~Gene,ncol=3,scales="free_y") +
    labs(x="Dose group",y="pN/pS (mean ± SD)") +
    theme_bw(base_size=SET_BASE_SIZE) +
    theme(
      axis.text.x       = element_text(size=SET_BASE_SIZE,face="bold"),
      axis.text.y       = element_text(size=SET_BASE_SIZE-2),
      axis.title        = element_text(size=SET_BASE_SIZE,face="bold"),
      strip.text        = element_text(face="bold.italic",size=SET_BASE_SIZE),
      strip.background  = element_rect(fill="gray95",color=NA),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.text       = element_text(size=SET_BASE_SIZE),
      legend.key.size   = unit(0.8,"cm"),
      panel.grid.minor  = element_blank(),
      panel.grid.major.y= element_line(color="gray90",linewidth=0.3),
      panel.background  = element_rect(fill="white",color=NA),
      plot.background   = element_rect(fill="white",color=NA),
      plot.title=element_blank(),plot.subtitle=element_blank(),
      plot.margin       = margin(10,15,10,15),
      panel.spacing     = unit(0.8,"lines")
    ) +
    guides(color=guide_legend(override.aes=list(size=5)),
           shape=guide_legend(override.aes=list(size=5)))
}


p_fig4 <- make_pnps_combined(res_A,res_B)
fig4_out <- paste0(out_pfx,"_Fig4_pNpS_TwoSets.png")
ggsave(fig4_out, p_fig4, width=14, height=8, dpi=300, bg="white")
cat("Figure 4 saved:", basename(fig4_out),"\n")

# ── Excel 저장 ────────────────────────────────────────────────────────────────
out_xlsx <- paste0(out_pfx,"_RESULTS.xlsx")
wb <- openxlsx::createWorkbook()
for (set_nm in c("A","B")) {
  res <- if (set_nm=="A") res_A else res_B
  openxlsx::addWorksheet(wb,paste0("Set",set_nm,"_IRR"))
  openxlsx::writeData(wb,paste0("Set",set_nm,"_IRR"),res$m1)
  openxlsx::addWorksheet(wb,paste0("Set",set_nm,"_Interaction"))
  openxlsx::writeData(wb,paste0("Set",set_nm,"_Interaction"),res$m2)
  openxlsx::addWorksheet(wb,paste0("Set",set_nm,"_pNpS"))
  pnps_out <- res$indiv_pnps %>%
    dplyr::select(Sample,Dose,Gene,Count_NS,Count_S,pN,pS,pNpS)
  openxlsx::writeData(wb,paste0("Set",set_nm,"_pNpS"),pnps_out)
}
openxlsx::saveWorkbook(wb,out_xlsx,overwrite=TRUE)
cat("Excel saved:", basename(out_xlsx),"\n")

cat("\n✓ Done!\n")
cat(sprintf("  Fig1 IRR (two sets) : %s\n",basename(fig1_out)))
cat(sprintf("  Fig4 pN/pS (two sets): %s\n",basename(fig4_out)))
cat(sprintf("  Excel               : %s\n",basename(out_xlsx)))