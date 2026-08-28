#Step 1: load the required data:
setwd("C:\\Users\\seabd4\\OneDrive - Université Laval\\Documents\\ULAVAL\\tcga\\GSE62944\\rhythmate")

rhythmate_cancer_data <- readRDS("rhythmate_cancer_data.rds")
rhythmate_normal_data <- readRDS("rhythmate_normal_data.rds")
rhythmate_cancer_clinical <- readRDS("rhythmate_cancer_clinical.rds")
rhythmate_normal_clinical <- readRDS("rhythmate_normal_clinical.rds")
rhythmate_survival_data <- read.xlsx("rhythmate_survival_data.xlsx")

#the quantile normalized data
rhythmate_cancer_qn <- readRDS("rhythmate_cancer_qn.rds")
rhythmate_normal_qn <- readRDS("rhythmate_normal_qn.rds")


group1 <- toupper(c("Clock","Arntl","Npas2","Nfil3","Cry1","Cry2", "Per1","Per2","Per3",
                    "Nr1d1","Nr1d2" ,"Bhlhe40","Dbp" ,"Hlf" ,"Tef"))
group1 <- append(group1, "C1orf51", after = which(group1 == "CRY2"))

group2 <- toupper(c("Rb1","Cdk4","E2f1","Pcna","Cdk1", "Ccna2","Cdc25c","Cenpe",
                    "Aurka", "Ccnb1", "Cdc20", "Rpa1", "Plk1", "Mcm2","Cdc6","Ccnd1" ))

group3 <- c("MYC", "MYCN", "MKI67")



#Step2: run run_rhythmate_normal_analysis in two objects named "out_norm" and "out_norm_cluster"
#define your desired gene names
#this function should be run twice: once with apply_log2_transform equal FALSE and once equal TRUE (base on 
#whether you have applied log2 to your expression data). 

run_rhythmate_normal_analysis <- function(
    rhythmate_normal_data,
    rhythmate_normal_clinical,
    group1,
    group2,
    group3 = NULL,
    iterations                       = 20,
    min_number_of_samples_per_subtype = 30,
    n_samples_per_iteration           = 24,
    seed                              = 130,
    correlation_method                = c("spearman", "pearson"),
    apply_log2_transform              = TRUE,
    exclude_self_tissue               = TRUE
) {
  
  #dependencies
  required_pkgs <- c("dplyr", "purrr", "tidyr", "tibble", "philentropy", "magrittr")
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(required_pkgs, library, character.only = TRUE))
  correlation_method <- match.arg(correlation_method)
  
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  # coerce and validate
  if (is.data.frame(rhythmate_normal_data)) {
    rhythmate_normal_data <- as.matrix(rhythmate_normal_data)
  }
  if (!is.matrix(rhythmate_normal_data)) stop("`rhythmate_normal_data` must be a matrix or coercible to matrix.")
  if (is.null(rownames(rhythmate_normal_data))) stop("`rhythmate_normal_data` must have gene symbols as rownames.")
  if (is.null(colnames(rhythmate_normal_data))) stop("`rhythmate_normal_data` must have barcodes as colnames.")
  
  # clinical table checks
  req_cols <- c("Cancer_type", "Barcode")
  if (!all(req_cols %in% colnames(rhythmate_normal_clinical))) {
    stop("`rhythmate_normal_clinical` must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  # gene lists
  if (!is.character(group1) || !is.character(group2)) stop("`group1` and `group2` must be character vectors of gene symbols.")
  if (length(group1) < 2) warning("`group1` has < 2 genes; within-group correlations will be trivial.")
  if (length(group2) < 2) warning("`group2` has < 2 genes; within-group correlations will be trivial.")
  
  # Optional group3
  if (!is.null(group3)) {
    if (!is.character(group3)) stop("`group3` must be a character vector of gene symbols.")
    g3 <- intersect(group3, rownames(rhythmate_normal_data))
    if (length(g3) == 0) {
      warning("None of the `group3` genes are present in the data; group3 will be ignored.")
      g3 <- character(0)
    }
  } else {
    g3 <- character(0)
  }
  
  # iterations / sampling
  iterations <- as.integer(iterations)
  n_samples_per_iteration <- as.integer(n_samples_per_iteration)
  min_number_of_samples_per_subtype <- as.integer(min_number_of_samples_per_subtype)
  if (iterations < 1L) stop("`iterations` must be >= 1.")
  if (n_samples_per_iteration < 2L) stop("`n_samples_per_iteration` must be >= 2.")
  if (min_number_of_samples_per_subtype < n_samples_per_iteration) {
    warning("`min_number_of_samples_per_subtype` < `n_samples_per_iteration`; some subtypes may be skipped if they cannot provide enough valid barcodes.")
  }
  
  #  helper functions 
  mean_expression_per_gene <- function(list_of_expr_mats) {
    all_cbind <- do.call(cbind, list_of_expr_mats)
    if (is.null(all_cbind) || ncol(all_cbind) == 0) return(setNames(numeric(0), character(0)))
    rowMeans(all_cbind, na.rm = TRUE)
  }
  
  reduceApplyListOfArrays <- function(list_of_mats) {
    # element-wise mean across a list of identically shaped matrices
    if (length(list_of_mats) == 0) return(NULL)
    ref_dim <- dim(list_of_mats[[1]])
    if (!all(vapply(list_of_mats, function(m) identical(dim(m), ref_dim), logical(1)))) {
      stop("All matrices in `list_of_mats` must have identical dimensions.")
    }
    Reduce("+", list_of_mats) / length(list_of_mats)
  }
  
  lower_tri_vector <- function(mat) {
    if (is.null(mat) || any(dim(mat) < 2)) return(numeric(0))
    mat[lower.tri(mat)]
  }
  pairwise_distance <- function(vec_a, vec_b) {
    if (!length(vec_a) || !length(vec_b) || length(vec_a) != length(vec_b)) return(NA_real_)
    as.numeric(philentropy::distance(rbind(vec_a, vec_b), method = "euclidean",
                                     use.row.names = FALSE, mute.message = TRUE))
  }
  
  # prepare subtype barcode list 
  dplyr <- asNamespace("dplyr")
  purrr <- asNamespace("purrr")
  tidyr <- asNamespace("tidyr")
  tibble <- asNamespace("tibble")
  
  barcode_list <-
    rhythmate_normal_clinical %>%
    dplyr::group_by(.data$Cancer_type) %>%
    dplyr::filter(dplyr::n() >= min_number_of_samples_per_subtype) %>%
    dplyr::group_split() %>%
    stats::setNames(
      rhythmate_normal_clinical %>%
        dplyr::group_by(.data$Cancer_type) %>%
        dplyr::filter(dplyr::n() >= min_number_of_samples_per_subtype) %>%
        dplyr::distinct(.data$Cancer_type) %>%
        dplyr::pull(.data$Cancer_type)) %>%
    purrr::map(~ .x$Barcode)
  
  set.seed(seed)
  
  # fix gene sets to what's available in the data 
  all_available <- intersect(rownames(rhythmate_normal_data), union(group1, group2))
  g1 <- intersect(group1, all_available)
  g2 <- intersect(group2, all_available)
  if (!length(g1)) stop("None of the `group1` genes are present in `rhythmate_normal_data`.")
  if (!length(g2)) stop("None of the `group2` genes are present in `rhythmate_normal_data`.")
  union_genes <- union(union(g1, g2), g3)
  
  # pre-compute row labels to make outputs readable
  cross_row_labels <- as.vector(outer(g1, g2, function(a, b) paste0(a, "::", b)))
  cir_row_labels   <- if (length(g1) >= 2) apply(utils::combn(g1, 2), 2, paste, collapse = "::") else character(0)
  cell_row_labels  <- if (length(g2) >= 2) apply(utils::combn(g2, 2), 2, paste, collapse = "::") else character(0)
  
  # containers
  exp.normal <- list()
  exp.normal.g3 <- list()
  normal     <- list()
  normal.cir <- list()
  normal.cell <- list()
  
  # per subtype iterations
  for (tyno in names(barcode_list)) {
    barcodes <- barcode_list[[tyno]]
    valid_barcodes <- intersect(barcodes, colnames(rhythmate_normal_data))
    
    if (length(valid_barcodes) >= n_samples_per_iteration) {
      exp.norm <- vector("list", iterations)
      exp.norm.g3 <- vector("list", iterations)
      normal.it <- vector("list", iterations)
      normal.it.cir <- vector("list", iterations)
      normal.it.cell <- vector("list", iterations)
      
      for (i in seq_len(iterations)) {
        sampled_barcodes <- sample(valid_barcodes, n_samples_per_iteration, replace = FALSE)
        n.samples <- rhythmate_normal_data[, sampled_barcodes, drop = FALSE]
        
        if (isTRUE(apply_log2_transform)) {
          exp.data <- log2(1 + n.samples[union_genes, , drop = FALSE])
        } else {
          exp.data <- n.samples[union_genes, , drop = FALSE]
        }
        
        g1_data <- exp.data[g1, , drop = FALSE]
        g2_data <- exp.data[g2, , drop = FALSE]
        
        if (length(g3) > 0) {
          g3_data <- exp.data[g3, , drop = FALSE]
        }
        
        tcganorm      <- stats::cor(t(g1_data), t(g2_data),
                                    method = correlation_method,
                                    use = "pairwise.complete.obs")
        tcganorm.cir  <- stats::cor(t(g1_data),
                                    method = correlation_method,
                                    use = "pairwise.complete.obs")
        tcganorm.cell <- stats::cor(t(g2_data),
                                    method = correlation_method,
                                    use = "pairwise.complete.obs")
        
        rownames(tcganorm)     <- rownames(g1_data)  
        colnames(tcganorm)     <- rownames(g2_data)  
        rownames(tcganorm.cir) <- rownames(g1_data)
        colnames(tcganorm.cir) <- rownames(g1_data)
        rownames(tcganorm.cell)<- rownames(g2_data)
        colnames(tcganorm.cell)<- rownames(g2_data)
        
        colnames(exp.data) <- paste0(colnames(exp.data), "_it", i)
        
        if (length(g3) > 0) {
          exp.norm.g3[[i]] <- g3_data
        }
        
        exp.norm[[i]]      <- exp.data
        normal.it[[i]]     <- tcganorm
        normal.it.cir[[i]] <- tcganorm.cir
        normal.it.cell[[i]]<- tcganorm.cell
      }
      
      if (length(g3) > 0) {
        exp.normal.g3[[tyno]] <- exp.norm.g3
      }
      
      if (length(normal.it) > 0) {
        exp.normal[[tyno]] <- exp.norm
        normal[[tyno]]     <- normal.it
        normal.cir[[tyno]] <- normal.it.cir
        normal.cell[[tyno]]<- normal.it.cell
      }
    } else {
      warning("Not enough valid samples for cancer type: ", tyno,
              " (have ", length(valid_barcodes), ", need >= ", n_samples_per_iteration, ")")
    }
  }
  
  # aggregate per subtype across iterations
  expr.normal  <- list()
  normm        <- list()
  normm.cir    <- list()
  normm.cell   <- list()
  normccd.cir  <- list()
  normccd.cell <- list()
  normccd.both <- list()
  disn.cir     <- list()
  disn.cell    <- list()
  disn.both    <- list()
  
  # cross-correlation matrix flattened per subtype
  normal_matrix <- NULL
  subtype_names <- character(0)
  
  for (tissue in names(normal)) {
    enormt    <- exp.normal[[tissue]]
    normt     <- normal[[tissue]]
    normt.cir <- normal.cir[[tissue]]
    normt.cell<- normal.cell[[tissue]]
    
    exp.normal.mean <- mean_expression_per_gene(enormt)
    norm            <- reduceApplyListOfArrays(normt)
    norm.cir        <- reduceApplyListOfArrays(normt.cir)
    norm.cell       <- reduceApplyListOfArrays(normt.cell)
    
    reshaped_vector <- as.vector(norm)
    normal_matrix   <- cbind(normal_matrix, reshaped_vector)
    subtype_names   <- c(subtype_names, tissue)
    
    expr.normal[[tissue]] <- exp.normal.mean
    normm[[tissue]]       <- norm
    normm.cir[[tissue]]   <- norm.cir
    normm.cell[[tissue]]  <- norm.cell
    
    normccd.cir[[tissue]]  <- lower_tri_vector(normm.cir[[tissue]])
    normccd.cell[[tissue]] <- lower_tri_vector(normm.cell[[tissue]])
    normccd.both[[tissue]] <- as.vector(normm[[tissue]])
  }
  
  # build the reference profile FROM THE NORMAL SAMPLES themselves 
  ref.mat.mean      <- reduceApplyListOfArrays(normm)       # Group1 x Group2, averaged over ALL normal subtypes
  ref.mat.mean.cir  <- reduceApplyListOfArrays(normm.cir)   # Group1 x Group1, averaged over ALL normal subtypes
  ref.mat.mean.cell <- reduceApplyListOfArrays(normm.cell)  # Group2 x Group2, averaged over ALL normal subtypes
  
  ref.ccd.cir.mean  <- if (!is.null(ref.mat.mean.cir))  lower_tri_vector(ref.mat.mean.cir)  else numeric(0)
  ref.ccd.cell.mean <- if (!is.null(ref.mat.mean.cell)) lower_tri_vector(ref.mat.mean.cell) else numeric(0)
  ref.ccd.both.mean <- if (!is.null(ref.mat.mean))      as.vector(ref.mat.mean)             else numeric(0)
  
  # distance of every normal subtype to each normal reference subtype
  for (tissue in names(normal)) {
    other_tissues <- if (isTRUE(exclude_self_tissue)) {
      setdiff(names(normccd.cir), tissue)
    } else {
      names(normccd.cir)
    }
    
    disn.cir[[tissue]]  <- vapply(other_tissues,
                                  function(ot) pairwise_distance(normccd.cir[[tissue]],  normccd.cir[[ot]]),  numeric(1))
    disn.cell[[tissue]] <- vapply(other_tissues,
                                  function(ot) pairwise_distance(normccd.cell[[tissue]], normccd.cell[[ot]]), numeric(1))
    disn.both[[tissue]] <- vapply(other_tissues,
                                  function(ot) pairwise_distance(normccd.both[[tissue]], normccd.both[[ot]]), numeric(1))
    
    
  }
  
  # build final data frames
  normal_dataframe <- if (!is.null(normal_matrix)) {
    normal_df <- as.data.frame(normal_matrix, stringsAsFactors = FALSE)
    colnames(normal_df) <- subtype_names
    rownames(normal_df) <- cross_row_labels
    normal_df
  } else {
    data.frame()
  }
  
  as_df_from_equal_length_vectors <- function(lst, row_labels) {
    if (!length(lst)) return(data.frame())
    lens <- vapply(lst, length, integer(1))
    if (!all(lens == lens[1])) stop("Within-group vectors differ in length across subtypes; check gene availability.")
    out <- as.data.frame(do.call(cbind, lst), stringsAsFactors = FALSE)
    colnames(out) <- names(lst)
    if (length(row_labels) == nrow(out)) rownames(out) <- row_labels
    out
  }
  
  normal_dataframe.cir  <- as_df_from_equal_length_vectors(normccd.cir,  cir_row_labels)
  normal_dataframe.cell <- as_df_from_equal_length_vectors(normccd.cell, cell_row_labels)
  
  # expression-by-sample long frame across all subtypes/iterations
  all_tissue_expr <- list()
  all_tissue_expr_g3 <- list()
  for (tissue in names(exp.normal)) {
    combined <- do.call(cbind, exp.normal[[tissue]])
    if (!is.null(combined) && ncol(combined) > 0) {
      df <- combined %>%
        tibble::as_tibble(rownames = "Gene", .name_repair = "unique") %>%
        tidyr::pivot_longer(-"Gene", names_to = "Sample", values_to = "Expression") %>%
        dplyr::filter(.data$Gene %in% union_genes) %>%
        dplyr::mutate(
          Type = dplyr::if_else(.data$Gene %in% g1, "Circadian", "CellCycle"),
          Subtype = tissue
        )
      all_tissue_expr[[tissue]] <- df
    }
    
    if (length(g3) > 0) {
      combined_g3 <- do.call(cbind, exp.normal.g3[[tissue]])
      if (!is.null(combined_g3) && ncol(combined_g3) > 0) {
        df3 <- combined_g3 %>%
          tibble::as_tibble(rownames = "Gene", .name_repair = "unique") %>%
          tidyr::pivot_longer(-"Gene", names_to = "Sample", values_to = "Expression") %>%
          dplyr::mutate(
            Type = "Group3",
            Subtype = tissue
          )
        all_tissue_expr_g3[[tissue]] <- df3
      }
    }
  }
  exp.normal.samples <- if (length(all_tissue_expr)) dplyr::bind_rows(all_tissue_expr) else tibble::tibble()
  
  exp.normal.samples.group3 <- 
    if (length(all_tissue_expr_g3)) dplyr::bind_rows(all_tissue_expr_g3)
  else tibble::tibble()
  
  # package distances as data frames (one row per subtype) 
  to_df_distance <- function(lst) {
    if (!length(lst)) return(data.frame())
    out <- data.frame(
      Subtype = names(lst),
      Distance = vapply(lst, function(v) mean(v, na.rm = TRUE), numeric(1)),
      stringsAsFactors = FALSE)
    out
  }
  disn.cir.df  <- to_df_distance(disn.cir)
  disn.cell.df <- to_df_distance(disn.cell)
  disn.both.df <- to_df_distance(disn.both)
  
  # return
  invisible(list(
    cor_matrix = normm,
    cor_matrix_cir = normm.cir,
    cor_matrix_cell = normm.cell,
    normal_dataframe      = normal_dataframe,
    normal_dataframe.cir  = normal_dataframe.cir,
    normal_dataframe.cell = normal_dataframe.cell,
    exp.normal.samples    = exp.normal.samples,
    exp.normal.samples.group3 = exp.normal.samples.group3,
    disn.cir              = disn.cir.df,
    disn.cell             = disn.cell.df,
    disn.both             = disn.both.df,
    ref.mat.mean          = ref.mat.mean,
    ref.mat.mean.cir       = ref.mat.mean.cir,
    ref.mat.mean.cell      = ref.mat.mean.cell,
    ref.ccd.cir.mean       = ref.ccd.cir.mean,
    ref.ccd.cell.mean      = ref.ccd.cell.mean,
    ref.ccd.both.mean      = ref.ccd.both.mean
  ))
}

#here we input quantile normalized data where log2 already applied
out_norm <- run_rhythmate_normal_analysis(
  rhythmate_normal_qn, rhythmate_normal_clinical, group1, group2, group3,
  iterations = 20, apply_log2_transform = FALSE)

#here we input NON-normalized data where log2 has not been applied
out_norm_cluster <- run_rhythmate_normal_analysis(
  rhythmate_normal_data, rhythmate_normal_clinical, group1, group2, group3,
  iterations = 20, apply_log2_transform = TRUE)



#Step3: run run_rhythmate_analysis in an object named "out_can" and "out_can_cluster"
#this function should be run twice: once with apply_log2_transform equal FALSE and once equal TRUE (base on 
#whether you have applied log2 to your expression data). 

run_rhythmate_analysis <- function(
    rhythmate_cancer_data,      
    rhythmate_cancer_clinical,  
    mat.mean,                  
    ccd.cir.mean,               
    ccd.cell.mean,               
    group1,                     
    group2,                     
    group3,                     
    iterations = 20,
    min_number_of_samples_per_subtype = 30,
    samples_per_subtype = 24,
    seed = 130,
    apply_log2_transform = TRUE,
    normal_reference_cross = NULL,
    normal_reference_cir   = NULL,
    normal_reference_cell  = NULL
) {
  # safety checks
  required_cols <- c("cancer.type","Subtype_Selected","com")
  optional_cols <- c("ID","sample_type_id")
  
  miss_required <- setdiff(required_cols, colnames(rhythmate_cancer_clinical))
  if (length(miss_required) > 0) {
    stop("Missing required clinical columns: ", paste(miss_required, collapse=", "))
  }
  for (col in optional_cols) {
    if (!(col %in% colnames(rhythmate_cancer_clinical))) {
      rhythmate_cancer_clinical[[col]] <- NA_character_
    }
  }
  if (is.null(rownames(rhythmate_cancer_data))) stop("rhythmate_cancer_data must have gene symbols as rownames.")
  if (is.null(colnames(rhythmate_cancer_data))) stop("rhythmate_cancer_data must have sample IDs as colnames.")
  
  
  required_pkgs <- c("dplyr", "reshape2", "tidyr", "tibble", "philentropy")
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(required_pkgs, library, character.only = TRUE))
  
  set.seed(seed)
  
  # gene handling 
  g1 <- toupper(as.character(group1))
  g2 <- toupper(as.character(group2))
  g3 <- toupper(as.character(group3))  
  symbols <- unique(c(g1, g2, g3))     # include g3 so we can export it later; excluded from correlations
  
  rn <- rownames(rhythmate_cancer_data)
  rn_up <- toupper(rn)
  s_types <- unique(rhythmate_cancer_clinical$cancer.type)
  
  # helpers
  count_samples <- function(type, subtype, clin) {
    if (is.null(subtype) || is.na(subtype)) return(0L)
    tab <- table(clin$cancer.type, clin$Subtype_Selected)
    if (!(type %in% rownames(tab)) || !(subtype %in% colnames(tab))) return(0L)
    as.integer(tab[type, subtype])
  }
  ch_sample <- function(mat, n) {
    mat[, sample(seq_len(ncol(mat)), n, replace = FALSE), drop = FALSE]
  }
  reshape_group_expr <- function(mat, subtype_name, group1, group2) {
    df <- data.frame(Gene = rownames(mat), mat, check.names = FALSE)
    long <- reshape2::melt(df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
    long <- long %>%
      dplyr::filter(.data$Gene %in% c(group1, group2)) %>%
      dplyr::mutate(
        Type = ifelse(.data$Gene %in% group1, "Circadian", "CellCycle"),
        Subtype = rep_len(subtype_name, nrow(long))
      )
    long
  }
  reshape_group3_expr <- function(mat, subtype_name, group3) {
    if (length(group3) == 0) return(NULL)
    df <- data.frame(Gene = rownames(mat), mat, check.names = FALSE)
    long <- reshape2::melt(df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
    long <- long %>%
      dplyr::filter(.data$Gene %in% group3) %>%
      dplyr::mutate(Type = "Group3", Subtype = rep_len(subtype_name, nrow(long)))
    long
  }
  sanitize_subtype <- function(type, subtype) {
    if (startsWith(subtype, paste0(type, "."))) return(subtype)
    sub_clean <- gsub("[^A-Za-z0-9]+", ".", subtype)
    paste0(type, ".", sub_clean)
  }
  ccd_to_reference <- function(obs_vec, ref_table = NULL, ref_mean_vec = NULL) {
    if (!is.null(ref_table)) {
      ref_mat <- as.matrix(ref_table)
      if (nrow(ref_mat) == length(obs_vec) && ncol(ref_mat) > 0) {
        dvals <- apply(ref_mat, 2, function(ref_vec) {
          as.numeric(philentropy::distance(rbind(as.numeric(ref_vec), as.numeric(obs_vec)),
                                           method = "euclidean"))
        })
        return(as.numeric(dvals))
      }
    }
    if (!is.null(ref_mean_vec) && length(ref_mean_vec) == length(obs_vec)) {
      return(as.numeric(philentropy::distance(rbind(as.numeric(ref_mean_vec), as.numeric(obs_vec)),
                                              method = "euclidean")))
    }
    numeric(0)
  }
  
  
  match_sample_ids <- function(clin_subset, expr_colnames) {
    common_id <- character(0)
    for (k in seq_len(nrow(clin_subset))) {
      base   <- clin_subset$com[k]
      id_val <- clin_subset$ID[k]
      type_id<- clin_subset$sample_type_id[k]
      if (is.na(base)) next
      
      candidate_ids <- character(0)
      if (!is.na(type_id)) {
        candidate_ids <- c(candidate_ids,
                           paste0(base, "-", sprintf("%02d", suppressWarnings(as.numeric(type_id)))),
                           paste0(base, "-", type_id))
      }
      if (!is.na(id_val)) candidate_ids <- c(candidate_ids, paste0(base, id_val))
      candidate_ids <- c(candidate_ids, base)
      
      matched <- intersect(candidate_ids, expr_colnames)
      if (length(matched) == 0) {
        matched <- expr_colnames[startsWith(expr_colnames, base)]
      }
      if (length(matched) > 0) common_id <- c(common_id, matched)
    }
    unique(common_id)
  }
  
  # containers 
  ls.exp.pr      <- list()
  ls.it.mat      <- list()
  ls.it.mat.cir  <- list()
  ls.it.mat.cell <- list()
  disss.cir      <- list()
  disss.cell     <- list()
  disss.both     <- list()
  
  ccd_scalar_cir_by_type  <- list()  
  ccd_scalar_cell_by_type <- list()
  ccd_scalar_both_by_type <- list()
  
  
  # main loop over cancer types 
  for (type in s_types) {
    subtypes <- rhythmate_cancer_clinical %>%
      dplyr::filter(.data$cancer.type == !!type, !is.na(.data$Subtype_Selected)) %>%
      dplyr::pull(.data$Subtype_Selected) %>% unique()
    
    exp.iter.list     <- vector("list", length = iterations)
    it.mat.iter       <- vector("list", length = iterations)    # Group1 x Group2
    it.mat.cir.iter   <- vector("list", length = iterations)    # Group1 x Group1
    it.mat.cell.iter  <- vector("list", length = iterations)    # Group2 x Group2
    
    for (i in seq_len(iterations)) {
      exp.pro     <- list()
      ls.cor      <- list()
      ls.cor.cir  <- list()
      ls.cor.cell <- list()
      ccd.cir     <- list()
      ccd.cell    <- list()
      
      for (sub in subtypes) {
        if (count_samples(type, sub, rhythmate_cancer_clinical) < min_number_of_samples_per_subtype) next
        
        clin_sub <- rhythmate_cancer_clinical %>%
          dplyr::filter(.data$cancer.type == !!type, .data$Subtype_Selected == !!sub)
        
        common_ids <- match_sample_ids(clin_sub, colnames(rhythmate_cancer_data))
        if (length(common_ids) < samples_per_subtype) next
        
        expr_sub <- rhythmate_cancer_data[, common_ids, drop = FALSE]
        rownames(expr_sub) <- rn_up
        expr_sub <- expr_sub[rownames(expr_sub) %in% symbols, , drop = FALSE]
        if (ncol(expr_sub) < samples_per_subtype) next
        
        expr_draw <- ch_sample(expr_sub, samples_per_subtype)
        
        g1_present <- intersect(g1, rownames(expr_draw))
        g2_present <- intersect(g2, rownames(expr_draw))
        if (length(g1_present) < 2 || length(g2_present) < 2) next
        
        g1_mat <- expr_draw[g1_present, , drop = FALSE]
        g2_mat <- expr_draw[g2_present, , drop = FALSE]
        
        cor_cross <- suppressWarnings(cor(t(g1_mat), t(g2_mat),
                                          method = "spearman", use = "pairwise.complete.obs"))
        cor_g1    <- suppressWarnings(cor(t(g1_mat),
                                          method = "spearman", use = "pairwise.complete.obs"))
        cor_g2    <- suppressWarnings(cor(t(g2_mat),
                                          method = "spearman", use = "pairwise.complete.obs"))
        
        key <- sanitize_subtype(type, sub)
        exp.pro[[key]]     <- expr_draw                # stores g1+g2+g3 rows
        ls.cor[[key]]      <- cor_cross
        ls.cor.cir[[key]]  <- cor_g1
        ls.cor.cell[[key]] <- cor_g2
        ccd.cir[[key]]     <- cor_g1[lower.tri(cor_g1)]
        ccd.cell[[key]]    <- cor_g2[lower.tri(cor_g2)]
        
        # BOTH OSCILLATORS (cross matrix): not symmetric, so no lower.tri
        # every (circadian gene, cell cycle gene) cell is a distinct value,
        # so flatten the whole matrix.
        both_vec_obs <- as.vector(cor_cross)
        both_vec_ref <- if (!is.null(mat.mean) && identical(dim(cor_cross), dim(mat.mean))) as.vector(mat.mean) else NULL
        s_both <- ccd_to_reference(both_vec_obs, ref_table = normal_reference_cross, ref_mean_vec = both_vec_ref)
        if (length(s_both) > 0) {
          if (is.null(ccd_scalar_both_by_type[[type]])) ccd_scalar_both_by_type[[type]] <- list()
          ccd_scalar_both_by_type[[type]][[key]] <- c(ccd_scalar_both_by_type[[type]][[key]], as.numeric(s_both))
        }
        
        # CIRCADIAN
        if (length(ccd.cir[[key]]) > 0) {
          s_cir <- ccd_to_reference(ccd.cir[[key]], ref_table = normal_reference_cir, ref_mean_vec = ccd.cir.mean)
          if (length(s_cir) > 0) {
            if (is.null(ccd_scalar_cir_by_type[[type]])) ccd_scalar_cir_by_type[[type]] <- list()
            ccd_scalar_cir_by_type[[type]][[key]] <- c(ccd_scalar_cir_by_type[[type]][[key]], as.numeric(s_cir))
          }
        }
        
        # CELL
        if (length(ccd.cell[[key]]) > 0) {
          s_cell <- ccd_to_reference(ccd.cell[[key]], ref_table = normal_reference_cell, ref_mean_vec = ccd.cell.mean)
          if (length(s_cell) > 0) {
            if (is.null(ccd_scalar_cell_by_type[[type]])) ccd_scalar_cell_by_type[[type]] <- list()
            ccd_scalar_cell_by_type[[type]][[key]] <- c(ccd_scalar_cell_by_type[[type]][[key]], as.numeric(s_cell))
          }
        }
        
      } # end subtypes loop
      
      exp.iter.list[[i]]     <- exp.pro
      it.mat.iter[[i]]       <- ls.cor
      it.mat.cir.iter[[i]]   <- ls.cor.cir
      it.mat.cell.iter[[i]]  <- ls.cor.cell
    } # end iterations
    
    # summarize CCD values to one mean per subtype for this cancer type
    final_cir <- NULL
    if (!is.null(ccd_scalar_cir_by_type[[type]])) {
      sub_names <- names(ccd_scalar_cir_by_type[[type]])
      if (length(sub_names) > 0) {
        vals <- vapply(sub_names, function(s) {
          mean(ccd_scalar_cir_by_type[[type]][[s]], na.rm = TRUE)
        }, numeric(1))
        final_cir <- vals
        names(final_cir) <- sub_names
      }
    }
    
    final_cell <- NULL
    if (!is.null(ccd_scalar_cell_by_type[[type]])) {
      sub_names <- names(ccd_scalar_cell_by_type[[type]])
      if (length(sub_names) > 0) {
        vals <- vapply(sub_names, function(s) {
          mean(ccd_scalar_cell_by_type[[type]][[s]], na.rm = TRUE)
        }, numeric(1))
        final_cell <- vals
        names(final_cell) <- sub_names
      }
    }
    
    final_both <- NULL
    if (!is.null(ccd_scalar_both_by_type[[type]])) {
      sub_names <- names(ccd_scalar_both_by_type[[type]])
      if (length(sub_names) > 0) {
        vals <- vapply(sub_names, function(s) {
          mean(ccd_scalar_both_by_type[[type]][[s]], na.rm = TRUE)
        }, numeric(1))
        final_both <- vals
        names(final_both) <- sub_names
      }
    }
    
    # store 
    disss.cir[[type]]  <- final_cir
    disss.cell[[type]] <- final_cell
    disss.both[[type]] <- final_both
    
    ls.exp.pr[[type]]      <- exp.iter.list
    ls.it.mat[[type]]      <- it.mat.iter
    ls.it.mat.cir[[type]]  <- it.mat.cir.iter
    ls.it.mat.cell[[type]] <- it.mat.cell.iter
  } # end types loop
  
  # build exp.tcga.samples 
  subtype_keys <- unique(unlist(lapply(ls.exp.pr, function(iter_list) {
    unique(unlist(lapply(iter_list, names), use.names = FALSE))
  }), use.names = FALSE))
  
  exp_long_list <- list()
  exp_long_list_g3 <- list()  
  
  for (k in subtype_keys) {
    mats_k <- list()
    for (type in names(ls.exp.pr)) {
      for (i in seq_along(ls.exp.pr[[type]])) {
        if (!is.null(ls.exp.pr[[type]][[i]][[k]])) {
          mats_k[[length(mats_k) + 1]] <- ls.exp.pr[[type]][[i]][[k]]
        }
      }
    }
    if (length(mats_k) == 0) next
    big_mat <- do.call(cbind, mats_k)
    
    keep_genes <- intersect(rownames(big_mat), c(g1, g2))
    if (length(keep_genes) > 0) {
      big_mat_g12 <- big_mat[keep_genes, , drop = FALSE]
      exp_long_list[[k]] <- reshape_group_expr(big_mat_g12, k, g1, g2)
    }
    
    keep_genes_g3 <- intersect(rownames(big_mat), g3)
    if (length(keep_genes_g3) > 0) {
      big_mat_g3 <- big_mat[keep_genes_g3, , drop = FALSE]
      exp_long_list_g3[[k]] <- reshape_group3_expr(big_mat_g3, k, g3)
    }
  }
  
  exp.tcga.samples <- dplyr::bind_rows(exp_long_list)
  if (nrow(exp.tcga.samples) > 0 && isTRUE(apply_log2_transform)) {
    exp.tcga.samples$Expression <- log2(exp.tcga.samples$Expression + 1)
  }
  exp.tcga.samples$Subtype <- as.character(exp.tcga.samples$Subtype)
  
  exp.tcga.samples.group3 <- dplyr::bind_rows(exp_long_list_g3)
  if (nrow(exp.tcga.samples.group3) > 0 && isTRUE(apply_log2_transform)) {
    exp.tcga.samples.group3$Expression <- log2(exp.tcga.samples.group3$Expression + 1)
  }
  exp.tcga.samples.group3$Subtype <- as.character(exp.tcga.samples.group3$Subtype)
  
  # build all.cancer / all.cancer.cir / all.cancer.cell
  mean_over_iterations <- function(ls_by_type_and_iter) {
    out <- list()
    for (type in names(ls_by_type_and_iter)) {
      iter_list <- ls_by_type_and_iter[[type]]
      if (length(iter_list) == 0 || !is.list(iter_list)) next
      subs <- unique(unlist(lapply(iter_list, function(x) {
        if (is.null(x) || length(x) == 0) return(character(0))
        names(x)
      }), use.names = FALSE))
      for (s in subs) {
        mats <- lapply(iter_list, function(x) {
          if (is.null(x) || length(x) == 0) return(NULL)
          x[[s]]
        })
        mats <- mats[!vapply(mats, is.null, logical(1))]
        if (length(mats) == 0) next
        ok <- vapply(mats, function(m) is.matrix(m) && nrow(m) > 0 && ncol(m) > 0, logical(1))
        mats <- mats[ok]
        if (length(mats) == 0) next
        dims <- unique(vapply(mats, function(m) paste(dim(m), collapse = "x"), character(1)))
        if (length(dims) > 1) next
        if (is.null(rownames(mats[[1]]))) rownames(mats[[1]]) <- paste0("r", seq_len(nrow(mats[[1]])))
        if (is.null(colnames(mats[[1]]))) colnames(mats[[1]]) <- paste0("c", seq_len(ncol(mats[[1]])))
        target_rn <- rownames(mats[[1]]); target_cn <- colnames(mats[[1]])
        mats <- lapply(mats, function(m) {
          if (is.null(rownames(m))) rownames(m) <- paste0("r", seq_len(nrow(m)))
          if (is.null(colnames(m))) colnames(m) <- paste0("c", seq_len(ncol(m)))
          m[target_rn, target_cn, drop = FALSE]
        })
        nr <- nrow(mats[[1]]); nc <- ncol(mats[[1]])
        arr <- array(unlist(mats, use.names = FALSE), dim = c(nr, nc, length(mats)))
        Mmean <- apply(arr, c(1, 2), mean, na.rm = TRUE)
        rownames(Mmean) <- target_rn; colnames(Mmean) <- target_cn
        out[[s]] <- Mmean
      }
    }
    out
  }
  
  mean_cross <- mean_over_iterations(ls.it.mat)       # Group1 x Group2
  mean_cir   <- mean_over_iterations(ls.it.mat.cir)   # Group1 x Group1
  mean_cell  <- mean_over_iterations(ls.it.mat.cell)  # Group2 x Group2
  
  .get_first_matrix <- function(L) {
    if (length(L) == 0) return(NULL)
    for (nm in names(L)) {
      M <- L[[nm]]
      if (is.matrix(M) && nrow(M) > 0 && ncol(M) > 0) return(M)
    }
    NULL
  }
  
  {
    first_cross <- .get_first_matrix(mean_cross)
    if (is.null(first_cross)) {
      all.cancer <- data.frame()
    } else {
      rn <- rownames(first_cross); cn <- colnames(first_cross)
      if (is.null(rn)) rn <- paste0("G1_", seq_len(nrow(first_cross)))
      if (is.null(cn)) cn <- paste0("G2_", seq_len(ncol(first_cross)))
      row_pairs <- as.vector(outer(rn, cn, FUN = function(a,b) paste(a, "-", b)))
      all.cancer <- as.data.frame(lapply(mean_cross, function(M) as.vector(M)))
      rownames(all.cancer) <- row_pairs
      colnames(all.cancer) <- names(mean_cross)
    }
  }
  {
    first_cir <- .get_first_matrix(mean_cir)
    if (is.null(first_cir)) {
      all.cancer.cir <- data.frame()
    } else {
      lt_idx <- lower.tri(first_cir)
      rn <- rownames(first_cir); cn <- colnames(first_cir)
      if (is.null(rn)) rn <- paste0("G1_", seq_len(nrow(first_cir)))
      if (is.null(cn)) cn <- paste0("G1_", seq_len(ncol(first_cir)))
      lt_names <- paste(rep(rn, each = nrow(first_cir))[lt_idx],
                        rep(cn, times = nrow(first_cir))[lt_idx], sep = " - ")
      all.cancer.cir <- as.data.frame(lapply(mean_cir, function(M) M[lower.tri(M)]))
      rownames(all.cancer.cir) <- lt_names
      colnames(all.cancer.cir) <- names(mean_cir)
    }
  }
  {
    first_cell <- .get_first_matrix(mean_cell)
    if (is.null(first_cell)) {
      all.cancer.cell <- data.frame()
    } else {
      lt_idx2 <- lower.tri(first_cell)
      rn <- rownames(first_cell); cn <- colnames(first_cell)
      if (is.null(rn)) rn <- paste0("G2_", seq_len(nrow(first_cell)))
      if (is.null(cn)) cn <- paste0("G2_", seq_len(ncol(first_cell)))
      lt_names2 <- paste(rep(rn, each = nrow(first_cell))[lt_idx2],
                         rep(cn, times = nrow(first_cell))[lt_idx2], sep = " - ")
      all.cancer.cell <- as.data.frame(lapply(mean_cell, function(M) M[lower.tri(M)]))
      rownames(all.cancer.cell) <- lt_names2
      colnames(all.cancer.cell) <- names(mean_cell)
    }
  }
  
  if (ncol(all.cancer) > 0)       colnames(all.cancer)      <- names(mean_cross)
  if (ncol(all.cancer.cir) > 0)   colnames(all.cancer.cir)  <- names(mean_cir)
  if (ncol(all.cancer.cell) > 0)  colnames(all.cancer.cell) <- names(mean_cell)
  
  
  # return 
  return(invisible(list(
    disss.cir        = disss.cir,
    disss.cell       = disss.cell,
    disss.both       = disss.both,
    cor_matrix = mean_cross,
    cor_matrix_cir = mean_cir,
    cor_matrix_cell = mean_cell,
    exp.tcga.samples = exp.tcga.samples,
    exp.tcga.samples.group3 = exp.tcga.samples.group3,
    all.cancer       = all.cancer,
    all.cancer.cir   = all.cancer.cir,
    all.cancer.cell  = all.cancer.cell
  )))
}

#here we input quantile normalized data where log2 already applied
out_can <- run_rhythmate_analysis(
  rhythmate_cancer_qn, rhythmate_cancer_clinical,
  out_norm$ref.mat.mean, out_norm$ref.ccd.cir.mean, out_norm$ref.ccd.cell.mean,
  group1, group2, group3,
  apply_log2_transform = FALSE,
  normal_reference_cross = out_norm$normal_dataframe,
  normal_reference_cir = out_norm$normal_dataframe.cir,
  normal_reference_cell = out_norm$normal_dataframe.cell)

#here we input NON-normalized data where log2 has not been applied
out_can_cluster <- run_rhythmate_analysis(
  rhythmate_cancer_data, rhythmate_cancer_clinical,
  out_norm_cluster$ref.mat.mean, out_norm_cluster$ref.ccd.cir.mean, out_norm_cluster$ref.ccd.cell.mean,
  group1, group2, group3,
  apply_log2_transform = TRUE,
  normal_reference_cross = out_norm_cluster$normal_dataframe,
  normal_reference_cir = out_norm_cluster$normal_dataframe.cir,
  normal_reference_cell = out_norm_cluster$normal_dataframe.cell)



# Step 4: run the following helper functions
rhythmate_survival_analysis <- function(rhythmate_survival_data, out_can,
                                        clusters_df, cluster_order) {
  message("Running survival analysis using clinical file: ", rhythmate_survival_data)
  
  first_non_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA else x[1]
  }
  
  #### TCGA CELL-CYCLE / CIRCADIAN RATIO PER SAMPLE
  mean_expr_by_sample_tcga <- out_can$exp.tcga.samples %>%
    group_by(Sample, Subtype, Type) %>%
    summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
  
  wide_expr_sample_tcga <- mean_expr_by_sample_tcga %>%
    pivot_wider(names_from = Type, values_from = MeanExpression)
  
  ratio_per_sample_tcga <- wide_expr_sample_tcga %>%
    filter(!is.na(CellCycle), !is.na(Circadian), Circadian != 0) %>%
    mutate(Ratio = CellCycle / Circadian, Log2Ratio = log2(Ratio))
  
  #### CLINICAL DATA
  rhythmate_survival_data <- read.xlsx(rhythmate_survival_data)
  
  #### GENE-LEVEL MARKERS (MYC / MYCN / MKI67) FOR TCGA SAMPLES
  ratio_per_Group3 <- out_can$exp.tcga.samples.group3 %>%
    pivot_wider(names_from = Gene, values_from = Expression)
  can.mark.3 <- ratio_per_Group3 %>%
    dplyr::mutate(com = str_sub(Sample, 1, 12), ID = str_sub(Sample, 14, 15))
  can.mark <- ratio_per_sample_tcga %>%
    dplyr::mutate(com = str_sub(Sample, 1, 12), ID = str_sub(Sample, 14, 15))
  can.mark <- can.mark %>%
    group_by(com) %>%
    summarise(Subtype = dplyr::first(Subtype),
              across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop")
  can.mark.3 <- can.mark.3 %>%
    group_by(com) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  
  can.mark.mer <- merge(can.mark, can.mark.3, by = "com", all = TRUE)
  
  #### SUBTYPE / CLINICAL MERGE
  subtypes <- PanCancerAtlas_subtypes()
  s.pan <- subtypes %>%
    dplyr::select(pan.samplesID, cancer.type, Subtype_mRNA, Subtype_Selected) %>%
    dplyr::mutate(com = str_sub(pan.samplesID, 1, 12), ID = str_sub(pan.samplesID, 14, 15))
  
  merged2 <- s.pan %>%
    left_join(rhythmate_survival_data, by = "com")
  
  can.mark.merged <- can.mark.mer %>%
    group_by(com) %>%
    summarise(Subtype = dplyr::first(Subtype),
              across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
  
  merged2 <- merged2 %>%
    left_join(can.mark.merged, by = "com") %>%
    dplyr::rename(circadian_mean = Circadian, cell_cycle_mean = CellCycle,
                  ratio = Ratio, log_ratio = Log2Ratio)
  
  #### SURVIVAL ENDPOINTS
  merged2 <- merged2 %>%
    dplyr::mutate(
      dsfs = case_when(DSS_STATUS == "1:DEAD WITH TUMOR" ~ 1, DSS_STATUS == "0:ALIVE OR DEAD TUMOR FREE" ~ 0),
      dsft = as.integer(merged2$DSS_MONTHS),
      dfs  = case_when(DFS_STATUS == "1:Recurred/Progressed" ~ 1, DFS_STATUS == "0:DiseaseFree" ~ 0),
      dft  = as.integer(merged2$DFS_MONTHS),
      oss  = case_when(OS_STATUS == "1:DECEASED" ~ 1, OS_STATUS == "0:LIVING" ~ 0),
      ost  = as.integer(merged2$OS_MONTHS),
      pfs  = case_when(PFS_STATUS == "1:PROGRESSION" ~ 1, PFS_STATUS == "0:CENSORED" ~ 0),
      pft  = as.integer(merged2$PFS_MONTHS))
  
  clinical.merge <- split(merged2, merged2$cancer.type)
  
  for (pro in names(clinical.merge)) {
    message(pro, "   ", paste(unique(clinical.merge[[pro]]$Subtype_Selected), collapse = ", "))
  }
  
  sanitize_subtype <- function(type, subtype) {
    mapply(function(ty, sub) {
      if (is.na(ty) || is.na(sub)) return(NA_character_)
      if (startsWith(sub, paste0(ty, "."))) return(sub)
      sub_clean <- gsub("[^A-Za-z0-9]+", ".", sub)
      paste0(ty, ".", sub_clean)
    }, type, subtype, USE.NAMES = FALSE)
  }
  
  mer.tcga1 <- list()
  for (pro in names(clinical.merge)) {
    tbl.tcga <- clinical.merge[[pro]] %>%
      dplyr::mutate(Subtype_Selected = sanitize_subtype(pro, Subtype_Selected)) %>%
      dplyr::select(pan.samplesID, CANCER_TYPE_ACRONYM, Subtype_Selected, Subtype,
                    dsfs, dsft, dfs, dft, oss, ost, pfs, pft, ratio, log_ratio,
                    circadian_mean, cell_cycle_mean, MYC, MYCN, MKI67)
    mer.tcga1[[pro]] <- tbl.tcga
  }
  
  mer.tcga2 <- data.table::rbindlist(mer.tcga1, fill = TRUE)
  mer.tcga2 <- as.data.frame(mer.tcga2)
  
  #### CLUSTER GROUPS: labeled C1 to Ck in the SAME left to right heatmap
  samples_by_cluster <- split(as.character(clusters_df$Sample), clusters_df$Cluster)
  k_eff <- length(cluster_order)
  
  group_list <- vector("list", k_eff)
  for (i in seq_len(k_eff)) {
    subs_i <- samples_by_cluster[[cluster_order[i]]]
    df_i <- mer.tcga2 %>% dplyr::filter(Subtype_Selected %in% subs_i)
    df_i$group   <- paste0("Group ", i)
    df_i$Cluster <- paste0("C", i)
    group_list[[i]] <- df_i
  }
  names(group_list) <- paste0("C", seq_len(k_eff))
  
  all.tcga      <- dplyr::bind_rows(group_list)
  combined_data <- all.tcga
  
  all.tcga$Cluster     <- factor(all.tcga$Cluster, levels = paste0("C", seq_len(k_eff)))
  combined_data$group  <- factor(combined_data$group, levels = paste0("Group ", seq_len(k_eff)))
  all.tcga$tissue_good <- as.character(lapply(strsplit(all.tcga$Subtype_Selected, "\\."), function(x) x[[1]]))
  
  #### COX MODELS
  fit_cox <- function(time_var, status_var) {
    cols <- c(time_var, status_var, "Cluster", "ratio", "tissue_good",
              "circadian_mean", "cell_cycle_mean", "MYC", "MYCN", "MKI67")
    dd <- all.tcga %>% dplyr::select(dplyr::all_of(cols)) %>% na.omit()
    f <- as.formula(paste0(
      "Surv(", time_var, ", ", status_var, ") ~ Cluster + ratio + tissue_good + ",
      "circadian_mean + cell_cycle_mean + MYC + MYCN + MKI67"))
    coxph(f, data = dd)
  }
  
  cox.os  <- fit_cox("ost",  "oss")
  cox.pf  <- fit_cox("pft",  "pfs")
  cox.df  <- fit_cox("dft",  "dfs")
  cox.dsf <- fit_cox("dsft", "dsfs")
  
  print(anova(cox.os));  print(summary(cox.os))
  print(anova(cox.pf));  print(summary(cox.pf))
  print(anova(cox.df));  print(summary(cox.df))
  print(anova(cox.dsf)); print(summary(cox.dsf))
  
  cox_models <- list(OS = cox.os, PFS = cox.pf, DFS = cox.df, DSS = cox.dsf)
  
  #### LOG-RANK TESTS
  ologrank_test  <- survdiff(Surv(ost, oss) ~ group, data = combined_data)
  plogrank_test  <- survdiff(Surv(pft, pfs) ~ group, data = combined_data)
  dlogrank_test  <- survdiff(Surv(dft, dfs) ~ group, data = combined_data)
  dslogrank_test <- survdiff(Surv(dsft, dsfs) ~ group, data = combined_data)
  
  oglobal_pval  <- ologrank_test$pvalue
  pglobal_pval  <- plogrank_test$pvalue
  dglobal_pval  <- dlogrank_test$pvalue
  dsglobal_pval <- dslogrank_test$pvalue
  
  logrank_pvals <- data.frame(
    Endpoint = c("Overall survival", "Progression-free survival",
                 "Disease-free survival", "Disease-specific survival"),
    GlobalLogRank_p = c(oglobal_pval, pglobal_pval, dglobal_pval, dsglobal_pval))
  
  #### KAPLAN-MEIER CURVES, ONE PANEL PER CLUSTER, PER ENDPOINT
  make_ssfit <- function(time_var, status_var) {
    f <- as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ 1"))
    fits <- lapply(group_list, function(df) survfit(f, data = df))
    names(fits) <- paste0("CC", seq_along(fits))
    fits
  }
  
  make_surv_plot <- function(time_var, status_var, ylab, main, global_pval) {
    ssfit <- make_ssfit(time_var, status_var)
    p <- ggsurvplot(ssfit, data = mer.tcga2, combine = TRUE, risk.table = TRUE, conf.int = FALSE,
                    conf.int.style = "step", censor = FALSE, tables.theme = theme_cleantable(),
                    palette = "jco", xlim = c(0, 60), break.time.by = 12,
                    xlab = "Time (month)", ylab = ylab, main = main)
    p$plot <- p$plot +
      annotate("text", x = 10, y = 0.2,
               label = paste("Global log-rank p =", format(global_pval, scientific = TRUE, digits = 2)),
               size = 5, color = "black")
    p
  }
  
  p_os  <- make_surv_plot("ost",  "oss",  "Overall-survival probability",
                          "Overall-survival curve", oglobal_pval)
  p_pfs <- make_surv_plot("pft",  "pfs",  "Progression-free survival probability",
                          "Progression-free survival curve", pglobal_pval)
  p_dfs <- make_surv_plot("dft",  "dfs",  "Disease-free survival probability",
                          "Disease-free survival curve", dglobal_pval)
  p_dss <- make_surv_plot("dsft", "dsfs", "Disease-specific survival probability",
                          "Disease-specific survival curve", dsglobal_pval)
  
  list(
    mer.tcga2     = mer.tcga2,
    all.tcga      = all.tcga,
    cox_models    = cox_models,
    logrank_pvals = logrank_pvals,
    plots         = list(os = p_os, pfs = p_pfs, dfs = p_dfs, dss = p_dss)
  )
}

rhythmate_export_cox_excel <- function(cox_models, logrank_pvals, xlsx_path) {
  wb <- openxlsx::createWorkbook()
  
  for (nm in names(cox_models)) {
    model <- cox_models[[nm]]
    
    tab <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
    names(tab)[names(tab) == "estimate"]  <- "HR"
    names(tab)[names(tab) == "conf.low"]  <- "HR_lower95"
    names(tab)[names(tab) == "conf.high"] <- "HR_upper95"
    
    smry <- summary(model)
    tab$n          <- smry$n
    tab$n_events   <- smry$nevent
    tab$concordance <- unname(smry$concordance["C"])
    
    anova_tab <- as.data.frame(anova(model))
    anova_tab <- cbind(Term = rownames(anova_tab), anova_tab)
    rownames(anova_tab) <- NULL
    
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, "Cox model (HR, 95% CI, p-values)", startRow = 1)
    openxlsx::writeData(wb, nm, tab, startRow = 2)
    
    anova_start_row <- 2 + nrow(tab) + 2
    openxlsx::writeData(wb, nm, "ANOVA (sequential likelihood ratio tests)", startRow = anova_start_row)
    openxlsx::writeData(wb, nm, anova_tab, startRow = anova_start_row + 1)
  }
  
  openxlsx::addWorksheet(wb, "LogRank_pvalues")
  openxlsx::writeData(wb, "LogRank_pvalues", logrank_pvals)
  
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
}


# Step 5: run run_rhythmate_clustering
# there are 3 output files with "rhythmate" as prefix:
# 1 ––> A multipage pdf file including:
# clustered heatmap, expression box plot, CCD plots and survival plots
# 2 --> An excel file including COX models 
# 3 --> An csv file inlcuding the cluster assignments

run_rhythmate_clustering <- function(out_norm, out_can, k = 5, out_prefix = "Rhythmate",
                                     rhythmate_survival_data = NULL,
                                     out_norm_cluster = NULL, out_can_cluster = NULL) {
  required_pkgs <- c("dplyr", "purrr", "tidyr", "tibble", "ComplexHeatmap", "tidyverse",
                     "circlize", "reshape2", "ggrepel", "survival", "survminer", "data.table",
                     "openxlsx", "broom", "TCGAbiolinks")
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(required_pkgs, library, character.only = TRUE))
  
  message("Starting Rhythmate clustering with k = ", k)
  
  #### 1) PREPARE MATRIX and CLUSTERING
  cluster_norm <- if (is.null(out_norm_cluster)) out_norm else out_norm_cluster
  cluster_can  <- if (is.null(out_can_cluster))  out_can  else out_can_cluster
  new.both <- cbind(cluster_norm$normal_dataframe,      cluster_can$all.cancer)
  new.cir  <- cbind(cluster_norm$normal_dataframe.cir,  cluster_can$all.cancer.cir)
  new.cell <- cbind(cluster_norm$normal_dataframe.cell, cluster_can$all.cancer.cell)
  all.three <- rbind(new.both, new.cir, new.cell)
  t.all.three <- t(all.three)
  
  d  <- dist(t.all.three)
  hc <- hclust(d, method = "ward.D2")
  clusters <- cutree(hc, k = k)
  
  clusters_df <- data.frame(
    Sample  = rownames(t.all.three),
    Cluster = as.factor(clusters))
  
  out_csv <- paste0(out_prefix, ".cluster_assignments_k", k, ".csv")
  write.csv(clusters_df, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
  
  #### 2) SAMPLE TYPE ANNOTATION
  n_normal <- ncol(cluster_norm$normal_dataframe)
  sample_types <- c(rep("Normal", n_normal),
                    rep("Cancer", ncol(cluster_can$all.cancer)))
  
  sample_anno <- HeatmapAnnotation(
    SampleType = sample_types,
    col = list(SampleType = c(Normal="#006666", Cancer="red")),
    annotation_legend_param = list(title = "Sample Type"),
    height = unit(3, "mm"))
  
  #### 3) EXTRACT HEATMAP ORDER
  n_both <- nrow(new.both); n_cir <- nrow(new.cir); n_cell <- nrow(new.cell)
  rowsplit <- c(
    rep("1.Both Oscillators", n_both),    rep("2.Circadian Clock",  n_cir),
    rep("3.Cell Cycle",       n_cell))
  
  ht_temp <- Heatmap(
    as.matrix(all.three),
    clustering_distance_columns = "euclidean",
    clustering_method_columns   = "ward.D2",
    column_split = factor(clusters),
    row_split    = rowsplit,
    show_row_names   = FALSE,
    show_column_dend = FALSE,
    show_row_dend    = FALSE)
  
  ht_drawn <- draw(ht_temp)
  ordered_cluster_indices <- column_order(ht_drawn)
  
  ordered_cluster_subtypes <- lapply(ordered_cluster_indices, function(idx) {
    colnames(all.three)[idx]
  })
  message("Extracted correct heatmap cluster order.")
  
  #left to right order of the ACTUAL cutree cluster ids, as they appear
  #on the heatmap. C1 would be the leftmost cluster.
  cluster_order <- names(ordered_cluster_indices)
  
  #### 4) COMPUTE CLUSTER-LEVEL MEAN LOG2 RATIOS (ORDERED)
  mean_expr_by_sample_tcga <- out_can$exp.tcga.samples %>%
    group_by(Sample, Subtype, Type) %>%
    summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
  
  wide_expr_sample_tcga <- mean_expr_by_sample_tcga %>%
    pivot_wider(names_from = Type, values_from = MeanExpression)
  
  ratio_per_sample_tcga <- wide_expr_sample_tcga %>%
    dplyr::filter(!is.na(CellCycle), !is.na(Circadian), Circadian != 0) %>%
    mutate(Ratio = CellCycle / Circadian, Log2Ratio = log2(Ratio))
  
  mean_expr_by_sample_normal <- out_norm$exp.normal.samples %>%
    group_by(Sample, Subtype, Type) %>%
    summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")
  
  wide_expr_sample_normal <- mean_expr_by_sample_normal %>%
    pivot_wider(names_from = Type, values_from = MeanExpression)
  
  ratio_per_sample_normal <- wide_expr_sample_normal %>%
    dplyr::filter(!is.na(CellCycle), !is.na(Circadian), Circadian != 0) %>%
    mutate(Ratio = CellCycle / Circadian, Log2Ratio = log2(Ratio))
  
  ratio_exp_per_sample_both_normal <- rbind(ratio_per_sample_tcga, ratio_per_sample_normal)
  
  overall_means <- sapply(ordered_cluster_subtypes, function(subs) {
    ratio_exp_per_sample_both_normal %>%
      dplyr::filter(Subtype %in% subs) %>%
      summarise(mean_log2_ratio = mean(Log2Ratio, na.rm = TRUE)) %>%
      pull(mean_log2_ratio)
  })
  
  block_colors <- ifelse(overall_means >= 0, "#FF6600", "#00CC99")
  log_labels   <- sprintf("%.3f", overall_means)
  
  log_anno <- HeatmapAnnotation(
    logratio = anno_block(
      gp = gpar(fill = block_colors),
      labels = log_labels,
      labels_gp = gpar(col = "white", fontsize = 12)))
  
  log_legend <- Legend(
    labels = c("Negative", "Positive"),
    title = "log2(cell cycle / circadian)",
    legend_gp = gpar(fill = c("#00CC99", "#FF6600")))
  
  if (!exists("col_fun")) {
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
  }
  
  ht_final <- Heatmap(
    as.matrix(all.three),
    name = "Spearman's rho",
    clustering_distance_columns = "euclidean",
    clustering_method_columns   = "ward.D2",
    column_split      = factor(clusters),
    top_annotation    = log_anno,
    bottom_annotation = sample_anno,
    row_split         = rowsplit,
    show_row_names    = FALSE,
    show_column_dend  = FALSE,
    show_row_dend     = FALSE,
    row_title_gp      = gpar(border = "black"),
    col               = col_fun)
  
  #### 6) PAGE 2: EXPRESSION BOXPLOT IN ORDER
  ordered_subtypes <- unlist(ordered_cluster_subtypes)
  
  melted_expr <- ratio_exp_per_sample_both_normal %>%
    pivot_longer(cols = c(Circadian, CellCycle),
                 names_to = "Type", values_to = "Expression")
  
  melted_expr$Subtype <- factor(melted_expr$Subtype, levels = ordered_subtypes)
  
  p_expr <- ggplot(melted_expr, aes(x = Subtype, y = Expression, color = Type)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
    scale_color_manual(values = c("#DC3220", "#005AB5"))+
    labs(title="Expression per Subtype", y="Normalized Expression")
  
  #### 7) PAGES 3-5: CIRCADIAN / CELL-CYCLE / BOTH-OSCILLATOR CCD
  normalize_named_numeric <- function(x) {
    if (is.null(x)) return(setNames(numeric(0), character(0)))
    if (is.data.frame(x) && all(c("Subtype", "Distance") %in% colnames(x))) {
      return(stats::setNames(as.numeric(x$Distance), as.character(x$Subtype)))
    }
    v <- if (is.list(x)) unlist(x, use.names = TRUE) else x
    v_num <- suppressWarnings(as.numeric(v))
    nm <- names(v)
    if (is.null(nm)) nm <- as.character(seq_along(v_num))
    names(v_num) <- nm
    v_num
  }
  normal_mean <- function(fin_df) {
    mean(fin_df$CCD[fin_df$Type == "Normal"], na.rm = TRUE)
  }
  normal_p95 <- function(fin_df) {
    stats::quantile(fin_df$CCD[fin_df$Type == "Normal"], probs = 0.95, na.rm = TRUE, names = FALSE)
  }
  normal_p05 <- function(fin_df) {
    stats::quantile(fin_df$CCD[fin_df$Type == "Normal"], probs = 0.05, na.rm = TRUE, names = FALSE)
  }
  
  flatten_disss_to_df <- function(disss) {
    if (is.null(disss) || length(disss) == 0) {
      return(data.frame(Type=character(0), Subtype=character(0), CCD=numeric(0),
                        stringsAsFactors = FALSE))
    }
    type_col <- character(0)
    subtype_col <- character(0)
    ccd_col <- numeric(0)
    
    for (ty in names(disss)) {
      v <- normalize_named_numeric(disss[[ty]])
      type_col    <- c(type_col, rep(ty, length(v)))
      subtype_col <- c(subtype_col, ifelse(is.na(names(v)) | names(v) == "", paste0(ty, "_", seq_along(v)), names(v)))
      ccd_col     <- c(ccd_col, as.numeric(v))
    }
    
    data.frame(Type = type_col, Subtype = subtype_col, CCD = ccd_col, stringsAsFactors = FALSE)
  }
  
  build_ccd_df <- function(disn, disss) {
    disn_vec <- normalize_named_numeric(disn)
    df_norm <- data.frame(
      Type = "Normal",
      Subtype = ifelse(is.na(names(disn_vec)) | names(disn_vec) == "",
                       paste0("Normal_", seq_along(disn_vec)),
                       names(disn_vec)),
      CCD = as.numeric(disn_vec),
      stringsAsFactors = FALSE)
    
    df_cancer <- flatten_disss_to_df(disss)
    
    fin <- rbind(df_norm, df_cancer)
    fin <- fin[!is.na(fin$CCD), , drop = FALSE]
    fin
  }
  
  order_types_by_median <- function(fin_df) {
    ord <- fin_df %>%
      dplyr::group_by(Type) %>%
      dplyr::summarise(med = median(CCD, na.rm = TRUE)) %>%
      dplyr::arrange(med) %>%
      dplyr::pull(Type)
    fin_df$Type <- factor(fin_df$Type, levels = ord)
    fin_df
  }
  
  build_cancer_colors <- function(all_types) {
    cancers <- setdiff(all_types, "Normal")
    cols <- if (length(cancers) > 0) {
      setNames(grDevices::hcl.colors(length(cancers), palette = "Dark 3"), cancers)
    } else {
      character(0)
    }
    c(Normal = "#006666", cols)
  }
  
  make_ccd_plot <- function(fin_df, title, col_map) {
    fin_df <- order_types_by_median(fin_df)
    norm_mean <- normal_mean(fin_df)
    norm_p95_val <- normal_p95(fin_df)
    norm_p05_val <- normal_p05(fin_df)
    
    ggplot(fin_df, aes(x = Type, y = CCD)) +
      geom_hline(yintercept = norm_mean, linetype = "dotted",
                 color = "black", linewidth = 0.6) +
      geom_hline(yintercept = c(norm_p95_val, norm_p05_val), linetype = "dotted",
                 color = "#555555", linewidth = 0.6) +
      geom_point(data = fin_df,
                 aes(color = Type),
                 position = position_jitter(width = 0.18, height = 0),
                 size = 2) +
      ggrepel::geom_text_repel(data = fin_df,
                               aes(label = Subtype, color = Type),
                               size = 2,
                               max.overlaps = 50,
                               position = position_jitter(width = 0.18, height = 0)) +
      scale_color_manual(values = col_map, drop = FALSE) +
      theme_minimal(base_size = 16) +
      labs(x = "", y = "CCD", title = title) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            legend.position = "none")
  }
  
  add_cluster_info <- function(fin_df, clusters_df, cluster_order) {
    merged <- dplyr::left_join(fin_df, clusters_df, by = c("Subtype" = "Sample"))
    n_dropped <- sum(is.na(merged$Cluster))
    if (n_dropped > 0) {
      warning(n_dropped, " rows had no matching Sample in clusters_df and will be dropped")
    }
    merged <- merged[!is.na(merged$Cluster), , drop = FALSE]
    merged$Cluster <- factor(merged$Cluster, levels = cluster_order)
    merged
  }
  
  make_ccd_plot_cluster <- function(fin_df, title, col_map, norm_mean, norm_p95_val, norm_p05_val) {
    ggplot(fin_df, aes(x = Cluster, y = CCD)) +
      geom_hline(yintercept = norm_mean, linetype = "dotted",
                 color = "black", linewidth = 0.6) +
      geom_hline(yintercept = c(norm_p95_val, norm_p05_val), linetype = "dotted",
                 color = "#555555", linewidth = 0.6) +
      geom_point(data = fin_df,
                 aes(color = Type),
                 position = position_jitter(width = 0.18, height = 0, seed = 42),
                 size = 2) +
      ggrepel::geom_text_repel(data = fin_df,
                               aes(label = Subtype, color = Type),
                               size = 2, max.overlaps = 50,
                               position = position_jitter(width = 0.18, height = 0, seed = 42)) +
      scale_color_manual(values = col_map, drop = FALSE) +
      theme_minimal(base_size = 16) +
      labs(x = " ", y = "CCD", title = title) +
      theme(axis.text.x = element_blank(),
            legend.position = "none")
  }
  
  # reference comes from out_norm_cluster (non-normalized samples)
  fin_cir  <- build_ccd_df(out_norm_cluster$disn.cir,  out_can_cluster$disss.cir)
  fin_cell <- build_ccd_df(out_norm_cluster$disn.cell, out_can_cluster$disss.cell)
  fin_both <- build_ccd_df(out_norm_cluster$disn.both, out_can_cluster$disss.both)
  
  all_types_union <- Reduce(union, list(unique(fin_cir$Type), unique(fin_cell$Type), unique(fin_both$Type)))
  col_map <- build_cancer_colors(all_types_union)
  
  p_cir  <- make_ccd_plot(fin_cir,  "Circadian Correlation Distance (CCD)", col_map)
  p_cell <- make_ccd_plot(fin_cell, "Cell-Cycle Correlation Distance (CCD)", col_map)
  p_both <- make_ccd_plot(fin_both, "Both-Oscillator Correlation Distance (CCD)", col_map)
  
  fin_cir_clust  <- add_cluster_info(fin_cir,  clusters_df, cluster_order)
  fin_cell_clust <- add_cluster_info(fin_cell, clusters_df, cluster_order)
  fin_both_clust <- add_cluster_info(fin_both, clusters_df, cluster_order)
  
  p_cir_cluster  <- make_ccd_plot_cluster(fin_cir_clust,  "Circadian CCD by Cluster", col_map, normal_mean(fin_cir), normal_p95(fin_cir), normal_p05(fin_cir))
  p_cell_cluster <- make_ccd_plot_cluster(fin_cell_clust, "Cell-Cycle CCD by Cluster", col_map, normal_mean(fin_cell), normal_p95(fin_cell), normal_p05(fin_cell))
  p_both_cluster <- make_ccd_plot_cluster(fin_both_clust, "Both-Oscillator CCD by Cluster", col_map, normal_mean(fin_both), normal_p95(fin_both), normal_p05(fin_both))
  
  #### 9) SURVIVAL ANALYSIS (only runs if a clinical csv path was given)
  survival_out <- NULL
  if (!is.null(rhythmate_survival_data)) {
    survival_out <- rhythmate_survival_analysis(
      rhythmate_survival_data     = rhythmate_survival_data,
      out_can                = out_can,
      clusters_df             = clusters_df,
      cluster_order            = cluster_order)
  }
  
  #### 10) OUTPUT MULTIPAGE PDF
  pdf_name <- paste0(out_prefix, "_FULL_RESULTS_k", k, ".pdf")
  pdf(pdf_name, width = 18, height = 10)
  draw(ht_final, annotation_legend_list = list(log_legend))# page 1
  print(p_expr) #page 2
  print(p_cir) #page 3
  print(p_cell) #page 4
  print(p_both) #page 5
  print(p_cir_cluster) #page 6
  print(p_cell_cluster) #page 7
  print(p_both_cluster) #page 8
  if (!is.null(survival_out)) {
    print(survival_out$plots$os)  #page 9 
    print(survival_out$plots$pfs) #page 10 
    print(survival_out$plots$dfs) #page 11 
    print(survival_out$plots$dss) #page 12
  }
  
  dev.off()
  message("Full multipage PDF saved: ", pdf_name)
  
  #### 11) EXPORT COX MODELS TO EXCEL
  if (!is.null(survival_out)) {
    xlsx_name <- paste0(out_prefix, "_cox_models_k", k, ".xlsx")
    rhythmate_export_cox_excel(survival_out$cox_models, survival_out$logrank_pvals, xlsx_name)
    message("Cox model results saved: ", xlsx_name)
  }
  
  message("COMPLETED.")
  
  invisible(list(
    cluster_assignments      = clusters_df,
    ordered_cluster_subtypes = ordered_cluster_subtypes,
    cluster_means            = overall_means,
    ccd_circadian             = fin_cir,
    ccd_cellcycle              = fin_cell,
    ccd_both                   = fin_both,
    survival                   = survival_out
  ))
}

run_rhythmate_clustering(
  out_norm, out_can,
  rhythmate_survival_data = "rhythmate_survival_data.xlsx",
  out_norm_cluster = out_norm_cluster,
  out_can_cluster = out_can_cluster)
