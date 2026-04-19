###############################################################################
## CARTWHEEL 모션캡쳐 분석 파이프라인 - (B) 국소 3D (torso-local XYZ)
##
##  과적합 대응 리팩토링 및 초고속 최적화:
##   1) 단일 파일 분석 (49_08.csv)
##   2) 16개 핵심 관절 마커만 사용 (Collinearity 제거)
##   3) 프레임 Downsample (연산량 반토막)
##   4) BFGS 알고리즘 및 trace=1 로 실시간 확인
##   5) 모형 적합도(Fit Indices) 콘솔 요약 및 자동 저장
##   6) print() 함수 partial matching 에러 리팩토링 완료
###############################################################################

pkgs <- c("MARSS", "lavaan", "dsem", "dplyr", "readr", "signal")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
    library(MARSS); library(lavaan); library(dsem)
    library(dplyr); library(readr)
})

## =========================================================================
## 0. 경로/옵션
## =========================================================================
FILES <- c(
    trial_49_08 = "C:/Users/yjl59/Downloads/49_08.csv"
)
OUT_DIR <- "C:/Users/yjl59/Downloads/results_xyz_local"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SMOOTH          <- TRUE
SG_N            <- 11; SG_P <- 3
DETREND         <- TRUE
GAP_FRAMES      <- 5
R_FORM          <- "diagonal and equal"
M_GRID          <- 3:4
USE_CORE_SUBSET <- TRUE
DOWNSAMPLE      <- 2

MARSS_CONTROL   <- list(maxit = 400, trace = 1, reltol = 1e-4)

save_csv <- function(df, path) {
    readr::write_csv(as.data.frame(df), path); cat("  saved:", path, "\n")
}


## =========================================================================
## 1) 국소 3D 좌표 변환
## =========================================================================
ANCHORS <- list(
    origin_candidates = list(c("PEL0"), c("LFWT","RFWT","LBWT","RBWT")),
    up_candidates     = list(c("C7"), c("CLAV")),
    front_pair        = c("STRN","T10")
)

.get_xyz <- function(dat, marker) {
    cols <- paste0(marker, c("_X","_Y","_Z"))
    if (!all(cols %in% names(dat))) return(NULL)
    as.matrix(dat[, cols])
}
.resolve_anchor <- function(dat, candidates) {
    for (cand in candidates) {
        mats <- lapply(cand, function(m) .get_xyz(dat, m))
        if (all(!sapply(mats, is.null))) {
            arr <- simplify2array(mats); return(apply(arr, c(1,2), mean))
        }
    }
    stop("앵커 마커를 찾지 못했습니다.")
}
.orthonormal_frame <- function(up, front) {
    y <- up / (sqrt(sum(up^2)) + 1e-12)
    x <- front - sum(front*y)*y
    x <- x / (sqrt(sum(x^2)) + 1e-12)
    z <- c(y[2]*x[3]-y[3]*x[2],
           y[3]*x[1]-y[1]*x[3],
           y[1]*x[2]-y[2]*x[1])
    cbind(x, y, z)
}
to_local_frame <- function(dat) {
    stopifnot("Frame" %in% names(dat))
    origin <- .resolve_anchor(dat, ANCHORS$origin_candidates)
    up_pt  <- .resolve_anchor(dat, ANCHORS$up_candidates)
    s <- .get_xyz(dat, ANCHORS$front_pair[1])
    e <- .get_xyz(dat, ANCHORS$front_pair[2])
    if (is.null(s) || is.null(e)) stop("front 앵커(STRN/T10) 필요.")
    Tn <- nrow(dat)
    marker_cols <- grep("_[XYZ]$", names(dat), value = TRUE)
    markers <- unique(sub("_[XYZ]$","", marker_cols))
    mk_mats <- setNames(lapply(markers, function(m) .get_xyz(dat, m)), markers)
    loc <- array(NA_real_, dim = c(Tn, length(markers), 3))
    for (t in seq_len(Tn)) {
        up    <- up_pt[t, ] - origin[t, ]; front <- s[t, ] - e[t, ]
        Rt <- tryCatch(.orthonormal_frame(up, front), error = function(e) NULL)
        if (is.null(Rt)) next
        for (k in seq_along(markers)) {
            mm <- mk_mats[[k]]; if (is.null(mm)) next
            loc[t, k, ] <- crossprod(Rt, mm[t, ] - origin[t, ])
        }
    }
    out <- dat["Frame"]
    for (k in seq_along(markers)) {
        m <- markers[k]
        out[[paste0(m,"_X")]] <- loc[,k,1]
        out[[paste0(m,"_Y")]] <- loc[,k,2]
        out[[paste0(m,"_Z")]] <- loc[,k,3]
    }
    out
}


## =========================================================================
## 2) 마커 선택 + 전처리
## =========================================================================
SUBSET_FULL_3D <- c(
    "RFIN","LFIN","RWRB","LWRB","RWRA","LWRA",
    "RFRM","LFRM","RELB","LELB","RUPA","LUPA",
    "STRN","CLAV","C7","T10",
    "LFWT","RFWT","LBWT","RBWT",
    "LTHI","RTHI","LKNE","RKNE","LSHN","RSHN","LANK","RANK"
)
SUBSET_CORE_3D <- c(
    "RFIN", "LFIN", "RWRB", "LWRB",
    "RELB", "LELB",
    "CLAV", "C7", "STRN", "T10",
    "LFWT", "RFWT",
    "RKNE", "LKNE",
    "RANK", "LANK"
)

MARKERS_3D <- if (USE_CORE_SUBSET) SUBSET_CORE_3D else SUBSET_FULL_3D
cat("\n[Markers] n =", length(MARKERS_3D),
    "→ channels =", 3 * length(MARKERS_3D), "(X/Y/Z local)\n")

preprocess_matrix <- function(mat) {
    mat <- apply(mat, 2, function(x) {
        if (all(is.na(x))) return(x)
        stats::approx(seq_along(x), x, xout = seq_along(x), rule = 2)$y
    })
    if (DETREND) {
        tt <- seq_len(nrow(mat))
        mat <- apply(mat, 2, function(x) residuals(lm(x ~ tt)))
    }
    if (SMOOTH && requireNamespace("signal", quietly = TRUE))
        mat <- apply(mat, 2, function(x) signal::sgolayfilt(x, p = SG_P, n = SG_N))
    scale(mat)
}

load_trial_xyz_local <- function(path, markers = MARKERS_3D) {
    raw <- readr::read_csv(path, show_col_types = FALSE)
    loc <- to_local_frame(raw)
    wanted <- c(outer(markers, c("X","Y","Z"), paste, sep = "_"))
    wanted <- wanted[wanted %in% names(loc)]
    if (!length(wanted)) stop("지정한 마커가 없습니다.")
    mat <- as.matrix(loc[, wanted, drop = FALSE])

    if (DOWNSAMPLE > 1) {
        idx <- seq(1, nrow(mat), by = DOWNSAMPLE)
        mat <- mat[idx, , drop = FALSE]
    }

    preprocess_matrix(mat)
}

trials_3d <- lapply(FILES, load_trial_xyz_local)
names(trials_3d) <- names(FILES)

cat("\n[Preprocessed]\n")
for (nm in names(trials_3d))
    cat(sprintf("  %-12s frames=%4d  channels=%3d\n",
                nm, nrow(trials_3d[[nm]]), ncol(trials_3d[[nm]])))


## =========================================================================
## 3) Pooled matrix
## =========================================================================
build_pooled_matrix <- function(trials, gap = GAP_FRAMES) {
    ch <- ncol(trials[[1]]); gap_mat <- matrix(NA_real_, nrow = gap, ncol = ch)
    pieces <- list(); ranges <- list(); pos <- 0L
    for (i in seq_along(trials)) {
        tm <- trials[[i]]
        ranges[[names(trials)[i]]] <- (pos + 1):(pos + nrow(tm))
        pieces[[length(pieces)+1]] <- tm; pos <- pos + nrow(tm)
        if (i < length(trials)) { pieces[[length(pieces)+1]] <- gap_mat; pos <- pos + gap }
    }
    pooled <- do.call(rbind, pieces); colnames(pooled) <- colnames(trials[[1]])
    list(mat = pooled, ranges = ranges)
}

pool_3d       <- build_pooled_matrix(trials_3d)
pooled_mat    <- pool_3d$mat
pooled_ranges <- pool_3d$ranges

n_ch <- ncol(pooled_mat); n_fr <- nrow(pooled_mat)
cat(sprintf("\n[Pooled]\n  pooled frames=%d  channels=%d  frame/channel=%.2f\n",
            n_fr, n_ch, n_fr / n_ch))


## =========================================================================
## 4) DFA 적합 + m 자동 선택
## =========================================================================
fit_dfa <- function(mat_time_by_ch, m, R_form = R_FORM) {
    Y <- t(mat_time_by_ch)
    MARSS::MARSS(Y,
                 model   = list(m = m, R = R_form),
                 form    = "dfa",
                 z.score = FALSE, demean = FALSE,
                 method  = "BFGS",
                 control = MARSS_CONTROL,
                 silent  = FALSE)
}

cat("\n[DFA] m 자동 선택 (AICc grid)\n")
aicc_rows <- list()
for (mm in M_GRID) {
    fit_mm <- try(fit_dfa(pooled_mat, m = mm), silent = TRUE)
    if (inherits(fit_mm, "try-error")) { cat(sprintf("  m=%d: fit failed\n", mm)); next }
    aicc_rows[[as.character(mm)]] <- data.frame(
        m = mm, logLik = fit_mm$logLik, AICc = fit_mm$AICc,
        converged = isTRUE(fit_mm$convergence == 0)
    )
    cat(sprintf("  m=%d  AICc=%.2f  logLik=%.2f  converged=%s\n",
                mm, fit_mm$AICc, fit_mm$logLik,
                isTRUE(fit_mm$convergence == 0)))
}
aic_tbl <- dplyr::bind_rows(aicc_rows)
save_csv(aic_tbl, file.path(OUT_DIR, "dfa_m_selection_AICc.csv"))

best_m <- aic_tbl$m[which.min(aic_tbl$AICc)]
cat("  → 선택된 m =", best_m, "\n")
fit_best <- fit_dfa(pooled_mat, m = best_m)

extract_loadings <- function(fit) {
    Z <- coef(fit, type = "matrix")$Z
    if (ncol(Z) >= 2) { rot <- varimax(Z); Z <- Z %*% rot$rotmat }
    rn <- rownames(fit$model$data); if (is.null(rn)) rn <- attr(fit$model$data, "Y.names")
    rownames(Z) <- rn; colnames(Z) <- paste0("F", seq_len(ncol(Z)))
    Z
}
extract_states <- function(fit) {
    s <- fit$states; rownames(s) <- paste0("F", seq_len(nrow(s))); s
}

Zpooled <- extract_loadings(fit_best)
Xpooled <- extract_states(fit_best)

prop_var <- apply(Xpooled, 1, var); prop_var <- prop_var / sum(prop_var)
cat("\n[Factor variance share]\n")
for (k in seq_along(prop_var)) cat(sprintf("  F%d: %.3f\n", k, prop_var[k]))


## =========================================================================
## 5) 저장: pooled loading + trial별 states
## =========================================================================
Z_df <- data.frame(
    channel = rownames(Zpooled),
    marker  = sub("_[XYZ]$", "", rownames(Zpooled)),
    axis    = sub(".*_", "", rownames(Zpooled)),
    as.data.frame(Zpooled), check.names = FALSE
)
save_csv(Z_df, file.path(OUT_DIR, "pooled_loadings.csv"))

for (nm in names(pooled_ranges)) {
    idx <- pooled_ranges[[nm]]
    Sk  <- t(Xpooled[, idx, drop = FALSE])
    save_csv(data.frame(frame = seq_along(idx), as.data.frame(Sk), check.names = FALSE),
             file.path(OUT_DIR, paste0(nm, "_states.csv")))
}


## =========================================================================
## 6) 손과의 cosine similarity
## =========================================================================
region_of <- function(marker) {
    dplyr::case_when(
        marker %in% c("RFIN","LFIN","RWRB","LWRB","RWRA","LWRA") ~ "Hand",
        marker %in% c("RFRM","LFRM","RELB","LELB","RUPA","LUPA") ~ "Upper",
        marker %in% c("STRN","CLAV","C7","T10",
                      "LFWT","RFWT","LBWT","RBWT")               ~ "Torso",
        marker %in% c("LTHI","RTHI","LKNE","RKNE",
                      "LSHN","RSHN","LANK","RANK")               ~ "Lower",
        TRUE ~ NA_character_
    )
}
hand_markers <- c("RFIN","LFIN","RWRB","LWRB","RWRA","LWRA")

similarity_to_hand <- function(Z, hand_markers) {
    mk <- sub("_[XYZ]$","", rownames(Z))
    Zm <- t(sapply(split(as.data.frame(Z), mk),
                   function(df) colMeans(as.matrix(df))))
    hp  <- colMeans(Zm[rownames(Zm) %in% hand_markers, , drop = FALSE])
    cos <- apply(Zm, 1, function(v) sum(v*hp) /
                     (sqrt(sum(v^2))*sqrt(sum(hp^2)) + 1e-12))
    data.frame(marker = rownames(Zm),
               region = region_of(rownames(Zm)),
               cos_sim = cos,
               is_hand = rownames(Zm) %in% hand_markers) |>
        dplyr::arrange(dplyr::desc(cos_sim))
}

sim_3d <- similarity_to_hand(Zpooled, hand_markers)
save_csv(sim_3d, file.path(OUT_DIR, "pooled_similarity.csv"))

# na.print 에러 방지: data.frame에 n = Inf 전달하면 print.default의 na.print로
# 부분매칭되어 "유효하지 않은 'na.print' 지정입니다" 에러 발생. 따라서 제거.
cat("\n[Local XYZ] 손과 loading 유사도 (pooled DFA):\n"); print(sim_3d)


## =========================================================================
## 7) Region scores + lavaan + dsem 패키지
## =========================================================================
region_scores_from_pool <- function(Z, X, idx) {
    Yhat <- Z %*% X[, idx, drop = FALSE]
    rgn  <- region_of(sub("_[XYZ]$","", rownames(Yhat)))
    df   <- as.data.frame(sapply(split(as.data.frame(Yhat), rgn),
                                 function(d) colMeans(as.matrix(d))))
    as.data.frame(scale(df))
}
regions_3d <- lapply(pooled_ranges,
                     function(idx) region_scores_from_pool(Zpooled, Xpooled, idx))

cat("\n[Save] region scores per trial\n")
for (nm in names(regions_3d)) {
    save_csv(data.frame(frame = seq_len(nrow(regions_3d[[nm]])), regions_3d[[nm]]),
             file.path(OUT_DIR, paste0(nm, "_region_scores.csv")))
}

build_cl_data <- function(df, lag = 1) {
    Tn <- nrow(df); cur <- df[1:(Tn-lag),,drop=FALSE]; nxt <- df[(1+lag):Tn,,drop=FALSE]
    names(cur) <- paste0(names(cur), "_t"); names(nxt) <- paste0(names(nxt), "_t1")
    cbind(cur, nxt)
}
cl_model <- '
  Hand_t1  ~ a1*Hand_t
  Upper_t1 ~ a2*Upper_t
  Lower_t1 ~ a3*Lower_t
  Torso_t1 ~ a4*Torso_t

  Hand_t1  ~ b1*Upper_t + b2*Lower_t + b3*Torso_t
  Upper_t1 ~ c1*Hand_t
  Lower_t1 ~ c2*Hand_t
  Torso_t1 ~ c3*Hand_t
'
fit_lavaan_to_tbl <- function(fit) {
    pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
    dplyr::select(pe, lhs, op, rhs, label, est, se, z, pvalue,
                  ci.lower, ci.upper, std.all)
}

cat("\n[Save] lavaan SEM per trial\n")
for (nm in names(regions_3d)) {
    d <- build_cl_data(regions_3d[[nm]], lag = 1)
    f <- lavaan::sem(cl_model, data = d, missing = "fiml")
    save_csv(fit_lavaan_to_tbl(f),
             file.path(OUT_DIR, paste0(nm, "_dsem_lavaan.csv")))
}
cl_data_pooled <- dplyr::bind_rows(lapply(regions_3d, build_cl_data), .id = "trial")
fit_cl_pooled  <- lavaan::sem(cl_model, data = cl_data_pooled, missing = "fiml")
save_csv(fit_lavaan_to_tbl(fit_cl_pooled),
         file.path(OUT_DIR, "all_dsem_lavaan_pooled.csv"))


## =========================================================================
## ★ 모형 적합도(Fit Indices) 콘솔 요약 및 저장
## =========================================================================
cat("\n[Model Fit] 핵심 적합도 지수 요약 (lavaan):\n")
fit_measures <- lavaan::fitMeasures(fit_cl_pooled,
                                    c("chisq", "df", "pvalue",
                                      "cfi", "tli", "rmsea", "srmr",
                                      "aic", "bic"))
print(round(fit_measures, 3))

fit_df <- as.data.frame(t(fit_measures))
save_csv(fit_df, file.path(OUT_DIR, "all_dsem_lavaan_fit_indices.csv"))


cat("\n[Local XYZ] Cross-lagged SEM (3 trial pooled):\n")
summary(fit_cl_pooled, standardized = TRUE, fit.measures = TRUE)


## ---- dsem 패키지 DSEM (trial별) ----
dsem_spec_xyz <- "
  Hand  -> Hand,  1, arH
  Upper -> Upper, 1, arU
  Lower -> Lower, 1, arL
  Torso -> Torso, 1, arT

  Upper -> Hand,  1, UtoH
  Lower -> Hand,  1, LtoH
  Torso -> Hand,  1, TtoH

  Hand  -> Upper, 1, HtoU
  Hand  -> Lower, 1, HtoL
  Hand  -> Torso, 1, HtoT

  Hand  <-> Hand,  0, v_H
  Upper <-> Upper, 0, v_U
  Lower <-> Lower, 0, v_L
  Torso <-> Torso, 0, v_T
"
fit_dsem_trial <- function(df_regions, sem_text = dsem_spec_xyz) {
    cols <- c("Hand","Upper","Torso","Lower")
    cols <- cols[cols %in% names(df_regions)]
    tsd  <- ts(as.matrix(df_regions[, cols, drop = FALSE]))
    dsem::dsem(
        sem = sem_text, tsdata = tsd,
        family = rep("fixed", ncol(tsd)),
        control = dsem::dsem_control(getsd = TRUE, quiet = TRUE)
    )
}
dsem_to_tbl <- function(fit_d) {
    s <- summary(fit_d); if (is.list(s) && !is.data.frame(s)) s <- as.data.frame(s); s
}

cat("\n[Save] dsem package DSEM per trial\n")
for (nm in names(regions_3d)) {
    res <- try(fit_dsem_trial(regions_3d[[nm]]), silent = TRUE)
    if (inherits(res, "try-error")) {
        cat("  ! dsem failed for", nm, "-", attr(res, "condition")$message, "\n"); next
    }
    save_csv(dsem_to_tbl(res), file.path(OUT_DIR, paste0(nm, "_dsem_package.csv")))
}

cat("\nDONE. Output directory:\n  ", OUT_DIR, "\n")
###############################################################################
