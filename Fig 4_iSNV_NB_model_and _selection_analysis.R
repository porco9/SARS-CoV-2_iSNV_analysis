################################################################################
# Two-Set dNdS Analysis
#
# 두 세트(A, B)를 각각 분석하여 동일한 figure에 나란히 표기
#
# Figure 1 (combined): IRR forest (Total/NS/S) + IRR ratio bar
#   - facet by Gene, 각 세트별 IRR + ratio
#
# Figure 4 (combined): pN/pS 평균 ± SD — Set A (a) / Set B (b) 동일 디자인
#   - S iSNV = 0인 샘플: genome-wide pS (ORF8/ORF10 제외) fallback 적용
#   - fallback 적용 샘플은 열린 기호(open symbol)로 구분
#
# 분석 유전자  : Spike, ORF1a, ORF1b, ORF3a, ORF8, ORF10 (6개)
# GW pS 계산  : Spike, ORF1a, ORF1b, ORF3a, ORF6 (ORF8/ORF10 제외)
#
# Input: 각 세트당 2개 파일
#   [1] iSNV summary  — Set A: 7__ (filtered),  Set B: 6__ (unfiltered)
#   [2] callable sites — Set A: 10__ (filtered), Set B: 9__ (unfiltered)
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

# 분석 대상 유전자 (Figure에 표시)
genes_of_interest <- c("Spike","ORF1a","ORF1b","ORF3a","ORF8","ORF10")
gene_ref_lengths  <- c(Spike=3822,ORF1a=13218,ORF1b=8085,
                       ORF3a=828,ORF8=366,ORF10=117)

# genome-wide pS 계산용 추가 유전자 (ORF6 포함, ORF8/ORF10 제외)
GW_GENES_EXTRA   <- c("ORF6")          # genes_of_interest에 없지만 GW pS에 포함
GW_GENE_LENGTHS  <- c(ORF6=186)
GW_EXCLUDE       <- c("ORF8","ORF10")  # GW pS 계산에서 제외

# ── 전역 색상/모양/레이블 정의 ─────────────────────────────────────────────────
SET_COLORS    <- c("A"="#c0392b", "B"="#2980b9")
SET_SHAPES    <- c("A"=16,        "B"=17)
SET_LABELS    <- c("A"="Filtered iSNVs",
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

# genes_to_read: 기본값은 genes_of_interest, GW용 추가 유전자도 읽을 때 확장
parse_sheet <- function(infile, sheet_name, meta,
                        genes_to_read=genes_of_interest,
                        lengths_to_use=gene_ref_lengths) {
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
  genes_avail <- intersect(genes_to_read,
                           names(lengths_to_use)[sapply(names(lengths_to_use), function(g)
                             paste0(g,"_iSNVs") %in% names(rd) & paste0(g,"_callable") %in% names(rd))])
  purrr::map_dfr(genes_avail, function(g) {
    rl <- lengths_to_use[[g]]
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
  
  # ── possible sites (분석 유전자 기준) ───────────────────────────────────────
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
  
  # ── NB regression용 count_data ───────────────────────────────────────────
  # log_offset: log(callable sites) — NB regression의 exposure 보정용
  # Possible_sites는 pN/pS 계산에만 사용하며 offset에 포함하지 않음
  count_data <- purrr::map_dfr(c("NS","S","Total"),
                               ~ parse_sheet(callable_file,.x,meta)) %>%
    dplyr::left_join(
      gene_possible %>% dplyr::mutate(Total_poss_all=Total_poss_N+Total_poss_S),
      by="Gene") %>%
    dplyr::mutate(
      log_offset=log(Callable_norm)
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
  
  # ── pN/pS: genome-wide pS fallback (S=0 샘플 처리) ───────────────────────
  # 분석 유전자 S/NS counts
  count_all <- purrr::map_dfr(c("NS","S"),
                              ~ parse_sheet(callable_file,.x,meta)) %>%
    dplyr::left_join(gene_possible,by="Gene")
  
  wide <- count_all %>%
    dplyr::filter(MutType %in% c("NS","S")) %>%
    dplyr::select(Sample,Dose,Gene,MutType,Count,Total_poss_N,Total_poss_S) %>%
    tidyr::pivot_wider(names_from=MutType,
                       values_from=c(Count,Total_poss_N,Total_poss_S),
                       names_sep="_")
  
  # ORF6 등 GW_GENES_EXTRA도 callable file에서 읽기
  all_gw_lengths <- c(gene_ref_lengths, GW_GENE_LENGTHS)
  all_gw_genes   <- c(genes_of_interest, GW_GENES_EXTRA)
  
  count_gw_extra <- purrr::map_dfr(c("NS","S"), function(mt) {
    tryCatch(
      parse_sheet(callable_file, mt, meta,
                  genes_to_read  = GW_GENES_EXTRA,
                  lengths_to_use = GW_GENE_LENGTHS),
      error=function(e) NULL)
  })
  
  # GW pS 계산용: 분석 유전자 + ORF6, GW_EXCLUDE 제외, S sheet만 사용
  count_S_gw <- dplyr::bind_rows(
    count_all %>% dplyr::filter(MutType == "S"),
    count_gw_extra %>% dplyr::filter(MutType == "S")
  ) %>%
    dplyr::filter(!Gene %in% GW_EXCLUDE)
  
  # 샘플별 genome-wide S iSNV 합계
  gw_S_counts <- count_S_gw %>%
    dplyr::group_by(Sample, Dose) %>%
    dplyr::summarise(
      gw_Count_S = sum(Count, na.rm=TRUE),
      .groups    = "drop"
    )
  
  # genome-wide possible S sites (GW_EXCLUDE 제외) — 고정값
  poss_gw <- dplyr::bind_rows(
    gene_possible %>% dplyr::filter(!Gene %in% GW_EXCLUDE),
    tryCatch({
      purrr::map_dfr(GW_GENES_EXTRA, function(g) {
        raw <- readxl::read_excel(isnv_file, sheet=g, col_names=FALSE)
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
      }) %>%
        dplyr::distinct(Gene,Position,Codon) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(sites=list(calc_possible_sites(Codon)),
                      poss_S=sites$S, poss_N=sites$N) %>%
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(poss_S)) %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(Total_poss_S=sum(poss_S),
                         Total_poss_N=sum(poss_N), .groups="drop")
    }, error=function(e) NULL)
  )
  
  gw_total_poss_S <- sum(poss_gw$Total_poss_S, na.rm=TRUE)
  cat(sprintf("  Genome-wide possible S sites (excl. ORF8/ORF10): %.1f\n",
              gw_total_poss_S))
  
  # 샘플별 genome-wide pS rate
  gw_pS <- gw_S_counts %>%
    dplyr::mutate(
      gw_pS_rate = ifelse(gw_total_poss_S > 0,
                          gw_Count_S / gw_total_poss_S, NA_real_)
    )
  
  # gene별 pN/pS: 모든 샘플에 genome-wide pS 통일 적용
  indiv_pnps <- wide %>%
    dplyr::left_join(gw_pS, by=c("Sample","Dose")) %>%
    dplyr::mutate(
      pN         = Count_NS / Total_poss_N_NS,
      pS         = gw_pS_rate,
      pS_source  = "genome_wide",
      pNpS       = ifelse(!is.na(pS) & pS > 0 & !is.na(pN), pN/pS, NA_real_),
      Dose_label = ifelse(Dose==0,"Pos",as.character(Dose))
    ) %>%
    dplyr::filter(!is.na(pNpS), is.finite(pNpS))
  
  cat(sprintf("  pN/pS computed (genome-wide pS): %d sample-gene observations\n",
              nrow(indiv_pnps)))
  
  # Dose factor levels
  has_pos  <- any(indiv_pnps$Dose==0)
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
cat("=== Set A (Filtered) ===\n")
cat("[A-1] Select iSNV summary  : 7__iSNV_summary...Syn분류.xlsx\n")
isnv_A <- file.choose()
cat("[A-2] Select callable sites: 10__최종정리본...xlsx\n")
callable_A <- file.choose()

cat("\n=== Set B (Unfiltered) ===\n")
cat("[B-1] Select iSNV summary  : 6__iSNV_summary...Syn분류.xlsx\n")
isnv_B <- file.choose()
cat("[B-2] Select callable sites: 9__최종정리본...xlsx\n")
callable_B <- file.choose()

out_dir <- dirname(callable_A)
out_pfx <- file.path(out_dir, "TwoSet_Analysis")

# ── 분석 실행 ─────────────────────────────────────────────────────────────────
res_A <- analyze_set(isnv_A, callable_A, "A")
res_B <- analyze_set(isnv_B, callable_B, "B")

# ════════════════════════════════════════════════════════════════
# FIGURE 1: IRR Forest + IRR Ratio
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
        plot.caption      = element_text(size=8, color="gray40", hjust=0),
        plot.title=element_blank(), plot.subtitle=element_blank(),
        plot.background=element_rect(fill="white",color=NA),
        panel.spacing=unit(0.8,"lines"), plot.margin=margin(8,12,8,12))

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
  
  all_ci <- c(df$CI_lo,df$CI_hi)
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
    geom_text(data=ratio_df,
              aes(x=xmx^0.95,y=y_pos,label=label,color=Set),
              hjust=1,size=3.5,fontface="bold",
              inherit.aes=FALSE,show.legend=FALSE) +
    scale_color_manual(values=SET_COLORS,labels=SET_LABELS,name=NULL) +
    scale_shape_manual(values=SET_SHAPES,labels=SET_LABELS,name=NULL) +
    scale_x_log10(limits=c(xmn,xmx),breaks=brk,
                  labels=scales::label_number(accuracy=0.1),
                  name="Incidence Rate Ratio (IRR, log scale)") +
    scale_y_continuous(breaks=c(1,2,3),labels=c("S","NS","Total"),
                       limits=c(0.5,3.8)) +
    facet_wrap(~Gene,ncol=3,scales="free_x") +
    theme_bw(base_size=SET_BASE_SIZE) +
    theme(
      axis.title.y       = element_blank(),
      axis.text.y        = element_text(size=SET_BASE_SIZE,face="bold"),
      axis.text.x        = element_text(size=SET_BASE_SIZE-2),
      strip.text         = element_text(face="bold",size=SET_BASE_SIZE),
      strip.background   = element_rect(fill="gray90"),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(size=SET_BASE_SIZE),
      legend.key.size    = unit(0.8,"cm"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color="gray90",linewidth=0.3),
      plot.background    = element_rect(fill="white",color=NA),
      plot.title=element_blank(),plot.subtitle=element_blank(),
      panel.spacing      = unit(1,"lines"),
      plot.margin        = margin(10,15,10,15)
    ) +
    guides(color=guide_legend(override.aes=list(size=5)),
           shape=guide_legend(override.aes=list(size=5)))
}

p_fig1 <- make_irr_combined(res_A, res_B)
fig1_out <- paste0(out_pfx,"_Fig1_IRR_TwoSets.png")
ggsave(fig1_out, p_fig1, width=16, height=9, dpi=300, bg="white")
cat("Figure 1 saved:", basename(fig1_out),"\n")

# ════════════════════════════════════════════════════════════════
# FIGURE 4: pN/pS 평균 ± SD
# ════════════════════════════════════════════════════════════════
cat("\nGenerating Figure 4 (pN/pS mean±SD, two sets)...\n")

pnps_to_bg <- function(val, alpha=0.35) {
  if (is.na(val)||!is.finite(val)) return(grDevices::adjustcolor("gray85",alpha))
  if (val<1) {
    intensity <- max(0,min(1,1-val))
    col <- grDevices::colorRamp(c("gray95","steelblue"))(intensity)
  } else {
    intensity <- max(0,min(1,(val-1)/1.5))
    col <- grDevices::colorRamp(c("gray95","#c0392b"))(intensity)
  }
  grDevices::adjustcolor(grDevices::rgb(col[1]/255,col[2]/255,col[3]/255),alpha)
}

pnps_to_bg_A <- function(val, alpha=0.45) {
  if (is.na(val)||!is.finite(val))
    return(grDevices::adjustcolor("gray90", alpha))
  if (val < 1) {
    intensity <- max(0, min(1, 1 - val))
    col <- grDevices::colorRamp(c("gray95","#2471a3"))(intensity)
  } else {
    intensity <- max(0, min(1, (val-1)/1.5))
    col <- grDevices::colorRamp(c("gray95","#c0392b"))(intensity)
  }
  grDevices::adjustcolor(
    grDevices::rgb(col[1]/255, col[2]/255, col[3]/255), alpha)
}

pnps_to_bg_B <- function(val, alpha=0.25) {
  if (is.na(val)||!is.finite(val))
    return(grDevices::adjustcolor("gray90", alpha))
  if (val < 1) {
    intensity <- max(0, min(1, 1 - val))
    col <- grDevices::colorRamp(c("gray95","#2471a3"))(intensity)
  } else {
    intensity <- max(0, min(1, (val-1)/1.5))
    col <- grDevices::colorRamp(c("gray95","#c0392b"))(intensity)
  }
  grDevices::adjustcolor(
    grDevices::rgb(col[1]/255, col[2]/255, col[3]/255), alpha)
}

# ── Fig4 본체 ─────────────────────────────────────────────────────────────────
make_pnps_combined <- function(res_a, res_b) {
  
  all_lvs_ordered <- c("Pos","0.01","0.1","1")
  dlvs <- all_lvs_ordered[all_lvs_ordered %in%
                            unique(c(res_a$dose_lvs, res_b$dose_lvs))]
  
  # 요약 통계
  prep_ms <- function(res, set_nm) {
    res$indiv_pnps %>%
      dplyr::group_by(Gene, Dose_f) %>%
      dplyr::summarise(
        mean_pNpS = mean(pNpS, na.rm=TRUE),
        sd_pNpS   = sd(pNpS,   na.rm=TRUE),
        .groups   = "drop"
      ) %>%
      dplyr::mutate(
        Gene      = factor(Gene, levels=genes_of_interest),
        x_pos_new = match(as.character(Dose_f), dlvs),
        x_jit     = x_pos_new + ifelse(set_nm=="A", -0.15, 0.15),
        Set       = set_nm,
        pt_color  = ifelse(set_nm=="A", "#c0392b", "#2980b9")
      ) %>%
      dplyr::filter(!is.na(x_pos_new))
  }
  
  ms_A   <- prep_ms(res_a, "A")
  ms_B   <- prep_ms(res_b, "B")
  ms_all <- dplyr::bind_rows(ms_A, ms_B)
  
  # 분할 배경색
  bg_A <- ms_A %>% dplyr::mutate(
    xmin=x_pos_new-0.5, xmax=x_pos_new,
    bg_col=mapply(pnps_to_bg_A, mean_pNpS))
  bg_B <- ms_B %>% dplyr::mutate(
    xmin=x_pos_new, xmax=x_pos_new+0.5,
    bg_col=mapply(pnps_to_bg_B, mean_pNpS))
  bg_df <- dplyr::bind_rows(bg_A, bg_B)
  
  # 개별 점
  set.seed(42)
  indiv_all <- dplyr::bind_rows(
    res_a$indiv_pnps %>% dplyr::mutate(Set="A", pt_color="#c0392b"),
    res_b$indiv_pnps %>% dplyr::mutate(Set="B", pt_color="#2980b9")
  ) %>%
    dplyr::mutate(
      Gene      = factor(Gene, levels=genes_of_interest),
      x_pos_new = match(as.character(Dose_f), dlvs),
      x_jit     = x_pos_new + ifelse(Set=="A",-0.15,0.15) +
        runif(dplyr::n(),-0.06,0.06),
      pt_shape  = ifelse(Set=="A", 16L, 17L),
      pt_fill   = pt_color
    ) %>%
    dplyr::filter(!is.na(x_pos_new))
  
  # KW (각 Set별 dose effect) — 빨강(Filtered)/파랑(Unfiltered)
  # ★ 수정: x_mid를 0.6으로 고정 (왼쪽 정렬)
  kw_df <- purrr::map_dfr(c("A","B"), function(set_nm) {
    res <- if(set_nm=="A") res_a else res_b
    purrr::map_dfr(genes_of_interest, function(g) {
      dat <- res$indiv_pnps %>% dplyr::filter(Gene==g, is.finite(pNpS))
      if(nrow(dat)<4||length(unique(dat$Dose_f))<2) return(NULL)
      kw <- tryCatch(kruskal.test(pNpS~Dose_f,data=dat),error=function(e)NULL)
      if(is.null(kw)) return(NULL)
      kp <- kw$p.value
      data.frame(Gene=g, Set=set_nm,
                 kp=kp,
                 sig=ifelse(kp<0.001,"***",ifelse(kp<0.01,"**",ifelse(kp<0.05,"*","ns"))),
                 stringsAsFactors=FALSE)
    })
  }) %>%
    dplyr::mutate(
      Gene    = factor(Gene, levels=genes_of_interest),
      sym     = "\u25c7",
      label   = paste0(sym, "p=", formatC(kp, format="f", digits=3), "(", sig, ")"),
      kw_col  = ifelse(Set=="A", "#c0392b", "#2980b9"),
      vjust_v = ifelse(Set=="A", 1.6, 3.4),
      x_mid   = 0.6   # ★ 수정: 패널 왼쪽 고정 (왼쪽 정렬)
    )
  
  # Wilcoxon signed-rank test (paired): Filtered vs Unfiltered
  wx_df <- purrr::map_dfr(genes_of_interest, function(g) {
    purrr::map_dfr(dlvs, function(dv) {
      da <- res_a$indiv_pnps %>%
        dplyr::filter(Gene==g, as.character(Dose_f)==dv, is.finite(pNpS)) %>%
        dplyr::select(Sample, pNpS_A=pNpS)
      db <- res_b$indiv_pnps %>%
        dplyr::filter(Gene==g, as.character(Dose_f)==dv, is.finite(pNpS)) %>%
        dplyr::select(Sample, pNpS_B=pNpS)
      paired <- dplyr::inner_join(da, db, by="Sample")
      if(nrow(paired) < 3) return(NULL)
      diff <- paired$pNpS_A - paired$pNpS_B
      if(all(diff == 0)) return(NULL)
      wx <- tryCatch(
        wilcox.test(paired$pNpS_A, paired$pNpS_B,
                    paired=TRUE, exact=FALSE),
        error=function(e) NULL)
      if(is.null(wx)) return(NULL)
      p <- wx$p.value
      ms_A_top <- ms_A %>%
        dplyr::filter(Gene==g, as.character(Dose_f)==dv) %>%
        dplyr::mutate(top=mean_pNpS+sd_pNpS) %>%
        dplyr::pull(top)
      ms_B_top <- ms_B %>%
        dplyr::filter(Gene==g, as.character(Dose_f)==dv) %>%
        dplyr::mutate(top=mean_pNpS+sd_pNpS) %>%
        dplyr::pull(top)
      y_top <- max(c(ms_A_top, ms_B_top, 0), na.rm=TRUE)
      data.frame(
        Gene      = g,
        Dose_lbl  = dv,
        x_pos_new = match(dv, dlvs),
        p_val     = p,
        sig       = ifelse(p<0.001,"***",ifelse(p<0.01,"**",
                                                ifelse(p<0.05,"*",""))),
        n_pairs   = nrow(paired),
        y_top     = y_top,
        stringsAsFactors=FALSE
      )
    })
  }) %>%
    dplyr::filter(sig != "") %>%
    dplyr::mutate(
      Gene      = factor(Gene, levels=genes_of_interest),
      # ★ 수정: bracket 위치를 더 위로
      y_bracket = y_top * 1.55,
      y_tick    = y_top * 1.43,
      y_label   = y_top * 1.68
    )
  
  ggplot() +
    geom_rect(data=bg_df,
              aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf),
              fill=bg_df$bg_col, inherit.aes=FALSE) +
    
    geom_hline(yintercept=1, linetype="dashed",
               color="gray40", linewidth=0.6) +
    
    # 개별 점
    geom_point(data=indiv_all, aes(x=x_jit, y=pNpS),
               color=indiv_all$pt_color,
               shape=indiv_all$pt_shape,
               fill =indiv_all$pt_fill,
               size=1.8, alpha=0.70, stroke=0.8,
               inherit.aes=FALSE) +
    
    # 에러바
    geom_errorbar(data=ms_all,
                  aes(x=x_jit,
                      ymin=pmax(mean_pNpS-sd_pNpS,0),
                      ymax=mean_pNpS+sd_pNpS),
                  color=ms_all$pt_color,
                  width=0.18, linewidth=0.9,
                  inherit.aes=FALSE) +
    
    # 평균 점 (열린 다이아몬드, 흰속)
    geom_point(data=ms_all, aes(x=x_jit, y=mean_pNpS),
               color=ms_all$pt_color,
               shape=23, fill="white",
               size=4.5, stroke=1.4,
               inherit.aes=FALSE) +
    
    # ★ 수정: KW p-value 왼쪽 정렬 (hjust=0 추가)
    geom_text(data=kw_df,
              aes(x=x_mid, y=Inf, label=label),
              vjust     = kw_df$vjust_v,
              hjust     = 0,
              color     = kw_df$kw_col,
              size=4.0, fontface="bold",
              inherit.aes=FALSE, show.legend=FALSE) +
    
    # Wilcoxon bracket with asterisks
    {if(nrow(wx_df)>0) {
      list(
        geom_segment(data=wx_df,
                     aes(x=x_pos_new-0.22, xend=x_pos_new+0.22,
                         y=y_bracket, yend=y_bracket),
                     linewidth=0.55, color="black",
                     inherit.aes=FALSE),
        geom_segment(data=wx_df,
                     aes(x=x_pos_new-0.22, xend=x_pos_new-0.22,
                         y=y_tick, yend=y_bracket),
                     linewidth=0.55, color="black",
                     inherit.aes=FALSE),
        geom_segment(data=wx_df,
                     aes(x=x_pos_new+0.22, xend=x_pos_new+0.22,
                         y=y_tick, yend=y_bracket),
                     linewidth=0.55, color="black",
                     inherit.aes=FALSE),
        geom_text(data=wx_df,
                  aes(x=x_pos_new, y=y_label, label=sig),
                  size=4.2, fontface="bold",
                  color="black", inherit.aes=FALSE,
                  show.legend=FALSE)
      )
    }} +
    
    scale_x_continuous(breaks=seq_along(dlvs), labels=dlvs,
                       expand=expansion(add=0.7)) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.45))) +
    facet_wrap(~Gene, ncol=3, scales="free_y") +
    labs(x="Dose group", y="pN/pS (mean \u00b1 SD)") +
    theme_bw(base_size=SET_BASE_SIZE) +
    theme(
      axis.text.x        = element_text(size=SET_BASE_SIZE,   face="bold"),
      axis.text.y        = element_text(size=SET_BASE_SIZE-2),
      axis.title         = element_text(size=SET_BASE_SIZE,   face="bold"),
      strip.text         = element_text(face="bold",        size=SET_BASE_SIZE),
      strip.background   = element_rect(fill="gray90"),
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color="gray90", linewidth=0.3),
      panel.background   = element_rect(fill="white",   color=NA),
      plot.background    = element_rect(fill="white",   color=NA),
      plot.title=element_blank(), plot.subtitle=element_blank(),
      plot.margin        = margin(10,15,10,15),
      panel.spacing      = unit(0.8,"lines")
    )
}

# ── legend 별도 패널 ─────────────────────────────────────────────────────────
# ★ 수정: plot 영역(정사각형 점) 완전히 숨기고 legend만 표시
make_legend_panel <- function() {
  df_leg <- data.frame(
    label = factor(c("Filtered iSNVs","Unfiltered iSNVs"),
                   levels=c("Filtered iSNVs","Unfiltered iSNVs")),
    x=c(1,2), y=c(1,1)
  )
  ggplot(df_leg, aes(x=x, y=y, color=label, fill=label, shape=label)) +
    geom_point(size=5) +
    scale_color_manual(values=c("Filtered iSNVs"="#c0392b",
                                "Unfiltered iSNVs"="#2980b9"),
                       name=NULL) +
    scale_fill_manual(values=c("Filtered iSNVs"="#c0392b",
                               "Unfiltered iSNVs"="#2980b9"),
                      name=NULL) +
    scale_shape_manual(values=c("Filtered iSNVs"=15,
                                "Unfiltered iSNVs"=15),
                       name=NULL) +
    # ★ plot 영역을 완전히 숨겨 정사각형 점이 보이지 않도록
    coord_cartesian(xlim=c(0,0), ylim=c(0,0), clip="off") +
    theme_void() +
    theme(
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.title       = element_blank(),
      legend.text        = element_text(size=SET_BASE_SIZE),
      legend.key.size    = unit(0.8,"cm"),
      plot.background    = element_rect(fill="white", color=NA),
      plot.margin        = margin(0,0,0,0)
    ) +
    guides(color=guide_legend(override.aes=list(size=5, shape=15)),
           shape=guide_legend(override.aes=list(size=5)),
           fill ="none")
}

# ── colorbar 별도 패널 ────────────────────────────────────────────────────────
# ★ 수정: barwidth, barheight 절반으로 축소
make_colorbar_panel <- function() {
  df_cb <- data.frame(x=1, y=seq(0, 2, length.out=200))
  ggplot(df_cb, aes(x=x, y=y, fill=y)) +
    geom_tile() +
    scale_fill_gradient2(
      low="steelblue", mid="white", high="#c0392b",
      midpoint=1, limits=c(0,2),
      name="pN/pS",
      breaks=c(0, 0.5, 1.0, 1.5, 2.0),
      guide=guide_colorbar(
        title.position="top",
        barwidth =unit(0.175,"cm"),   # ★ 수정: 0.35 → 0.175 (절반)
        barheight=unit(1.5,"cm"),     # ★ 수정: 3.0  → 1.5   (절반)
        ticks.colour="gray40"
      )
    ) +
    theme_void() +
    theme(
      legend.position="right",
      legend.title=element_text(size=SET_BASE_SIZE-2, face="bold"),
      legend.text =element_text(size=SET_BASE_SIZE-3),
      plot.background=element_rect(fill="white", color=NA)
    )
}

# ── 조합 출력 ─────────────────────────────────────────────────────────────────
p_main <- make_pnps_combined(res_A, res_B)
p_cb   <- make_colorbar_panel()
p_leg  <- make_legend_panel()

# ★ 수정: legend 패널 높이 0.06 → 0.04 (더 위로)
p_fig4 <- patchwork::wrap_plots(
  patchwork::wrap_plots(p_main, p_cb, ncol=2, widths=c(1, 0.06)),
  patchwork::wrap_plots(p_leg),
  ncol=1, heights=c(1, 0.04)
)

fig4_out <- paste0(out_pfx,"_Fig4_pNpS_TwoSets.png")
ggsave(fig4_out, p_fig4, width=16, height=10, dpi=300, bg="white")
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
    dplyr::select(Sample,Dose,Gene,Count_NS,Count_S,pN,pS,pNpS,pS_source)
  openxlsx::writeData(wb,paste0("Set",set_nm,"_pNpS"),pnps_out)
}
openxlsx::saveWorkbook(wb,out_xlsx,overwrite=TRUE)
cat("Excel saved:", basename(out_xlsx),"\n")

cat("\n\u2713 Done!\n")
cat(sprintf("  Fig1 IRR (two sets)  : %s\n",basename(fig1_out)))
cat(sprintf("  Fig4 pN/pS (two sets): %s\n",basename(fig4_out)))
cat(sprintf("  Excel                : %s\n",basename(out_xlsx)))
