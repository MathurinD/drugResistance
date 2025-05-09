#' @import synergyfinder

#' @title Get a valid input for ReshapeData
#'
#' @description Get a synergy table that can be used by the function ReshapeData synergyfinder.
#' @param pdata A tibble as outputed from process_growth_curves. Must have columns 'Treatment' with format drugRow_rowConcentration+drugCol_colConcentation and 'Viability' 
#' @param restrict Whether the data should be strictly restricted between 0 and 100. Apply a cutoff, not a scaling.
#' @param col_drug A drug name that must be the column drug. Otherwise choosen as the first drug to appear.
#' @param control The name of the control.
#' @return A tibble with column names as required by ReshapeData
#' @export
get_synergy_table <- function(pdata, restrict=FALSE, col_drug="", control="DMSO") {
    #pigfrazd %>% filter(!grepl("1|8", Well)) %>%
    pdata = pdata %>%
        mutate(Viability=case_when(is.infinite(Viability)&Viability<0~0, is.infinite(Viability)&Viability>0~1, TRUE~Viability)) %>% # Remove infinite
        mutate(Response = 100*Viability) # Convert Viability in ~[0,1] to a reponse in ~[0,100]
    if (restrict) { # Scale the Response strictly between 0 and 100
        pdata = pdata %>%
            mutate(Response = case_when(Response<0~0, Response>100~100, TRUE~Response))
    }

    # Convert the treatments to 2 columns with the names expected by synergyfinder
    treatments = str_split_fixed(pdata$Treatment, "\\+|_and_", 2) %>% as_tibble %>% rename(Drug1=V1, Drug2=V2) %>% separate(Drug1, c("DrugRow", "ConcRow"), "_") %>% separate(Drug2, c("DrugCol", "ConcCol"), "_")
    drugs_col = treatments %>% distinct(DrugCol) %>% pull(DrugCol)
    drugs_row = treatments %>% distinct(DrugRow) %>% pull(DrugRow)
    if (col_drug == "") {
        col_drug = find_first_drug(drugs_row, control)
    }
    treatments %>% apply(1, function(x){
                            if (x["DrugRow"] == col_drug) { # Swap
                                tmp = x[c("DrugRow", "ConcRow")]
                                x[c("DrugRow", "ConcRow")] = x[c("DrugCol", "ConcCol")]
                                x[c("DrugCol", "ConcCol")] = tmp
                            }
                            if (x["DrugRow"] %in% control) {
                                x[c("DrugRow", "ConcRow")]=c(NA, 0)
                            }
                            return(x) }) %>% t %>% as.tibble %>%
         mutate(ConcRow=case_when(is.na(ConcRow)~"0nM", TRUE~ConcRow), ConcCol=case_when(is.na(ConcCol)~"0nM", TRUE~ConcCol)) %>%
         mutate(ConcRow = concentration_values(ConcRow)*1e9, ConcRowUnit = "nM", ConcCol = concentration_values(ConcCol)*1e9, ConcColUnit = "nM") -> treatments
     row_drug = treatments %>% pull(DrugRow) %>% find_first_drug(c(col_drug, control))
     treatments = treatments %>% mutate(DrugRow=row_drug, DrugCol=col_drug)

    synergy_data = pdata %>%
        bind_cols(treatments) %>% # Add the treatment columns
        group_by(DrugCol, ConcCol, DrugRow, ConcRow) %>%
        summarize(Response=mean(Response)) %>% mutate(BlockId=1, ConcColUnit="nM", ConcRowUnit="nM") %>% # Make sure each concentration combination corresponds to 1 data point (not done earlier because controls can have different concentrations thus different names but all correspond to the same drug combination of 0)
        mutate(BlockID=1) %>% # For synergyFinder
        ungroup() %>%
        select(BlockID, DrugCol, ConcCol, ConcColUnit, DrugRow, ConcRow, ConcRowUnit, Response) # Sanity check that all column expected by synergyfinder are present

    return(synergy_data)
}

find_first_drug <- function(drugs_list, exclude=c("DMSO","CTRL")) {
    for (cc in drugs_list) {
        if (! cc %in% c("", exclude) & !is.na(cc)) {
            return(cc)
        }
    }
    return(NA)
}

#' Helper to get bliss score matrix
#'
#' Helper to get a nice bliss score matrix from processed data
#' @param pdata Proccessed incucyte data as outputed by process_growth_curves
#' @param restrict Whether the data should be restricted between 0 and 100
#' @export
bliss_score <- function(pdata, restrict=TRUE, col_drug="", control="DMSO") {
    pdata %>% get_synergy_table(restrict=restrict, col_drug=col_drug, control=control) %>% mutate(ConcUnit=ConcColUnit) %>% ReshapeData(noise=FALSE) -> drm
    if (!is.null(drm$response_statistics)) { # Make it flexible to synergyfinder version
        drm$response_statistics %>% mutate(response=response_mean) -> shaped_data
        drm$response_matrix = drm$response_statistics %>% dplyr::select(conc1, conc2, response_mean) %>% pivot_wider(values_from=response_mean, names_from=conc1) %>% column_to_rownames('conc2')
        Bliss(shaped_data) %>% dplyr::select(conc1, conc2, Bliss_synergy) %>% pivot_wider(values_from=Bliss_synergy, names_from=conc1) %>% column_to_rownames('conc2') -> drm$synergy # From synergyFinder
    } else if (!is.null(drm$response)) {
        drm$response -> shaped_data
        drm$response_matrix = drm$response %>% dplyr::select(conc1, conc2, response) %>% pivot_wider(values_from=response, names_from=conc1) %>% column_to_rownames('conc2')
        Bliss(shaped_data) %>% dplyr::select(conc1, conc2, Bliss_synergy) %>% pivot_wider(values_from=Bliss_synergy, names_from=conc1) %>% column_to_rownames('conc2') -> drm$synergy # From synergyFinder
    } else {
        drm %>% .$dose.response.mats %>% .[[1]] -> shaped_data
        Bliss(shaped_data) -> drm$synergy # From synergyFinder
        drm$response_matrix = drm$dose.response.mats[[1]]
    }
   # azd = shaped_data[1,]/100
   # aew = shaped_data[,1]/100
   # matrix(rep(azd, length(aew)), nrow=length(aew), byrow=TRUE) + matrix(rep(aew, length(azd)), ncol=length(azd)) - t(t(aew)) %*% t(azd) -> ybliss
    #return(synergy)
    return(drm)
}

#' Plot starHeatmaps for a group of synergy data
#' @param pdata_list Processed data or named list thereof. Should contain columns Well, Treatment, Ref_T, Elapsed and Analysis_Job
#' @param restrict Whether the data should be restricted between 0 and 100%.
#' @param col_drug Force a drug name to be in the columns.
#' @param control Name of the control condition.
#' @param ... Extra arguments to pass to starHeatmap2
#' @export
plot_bliss_scores <- function(pdata_list, restrict=TRUE, col_drug="", control="DMSO", name="Drug effect", ...) {
    if (!'list' %in% class(pdata_list)) {
        pdata_list = list("all"=pdata_list)
    }
    template_matrix = do.call(bind_rows, pdata_list) %>% bliss_score(restrict=restrict, col_drug=col_drug, control=control) %>% .$response_matrix  # Compute the maximum matrix size to scale all others
    hl = c()
    annotation = c()
    for (pp in names(pdata_list)) {
        bss = bliss_score(pdata_list[[pp]], restrict=restrict, col_drug=col_drug, control=control)
        if (suppressWarnings(!is.null(bss$drug.pairs$drug.col))) { drug_col = bss$drug.pairs$drug.col } else if (suppressWarnings(!is.null(bss$drug_pairs$drug1))) { drug_col = bss$drug_pairs$drug1} else { drug_col = bss$drug.pairs$drug_col } # To work with both versions of SynergyFinder
        if (suppressWarnings(!is.null(bss$drug.pairs$drug.row))) { drug_row = bss$drug.pairs$drug.row } else if (suppressWarnings(!is.null(bss$drug_pairs$drug2))) { drug_row = bss$drug_pairs$drug2} else { drug_row = bss$drug.pairs$drug_row } # To work with both versions of SynergyFinder
        bliss_heatmap = starHeatmap2(match_matrices(bss$response_matrix,template_matrix), match_matrices(bss$synergy,template_matrix), name=pp, column_title=drug_col, row_title=drug_row, return_list=TRUE, heatmap_legend_param = list(title = name), ...)
        hl = hl + bliss_heatmap$object
        annotation = bliss_heatmap$annotation_legend_list
    }
    draw(hl, padding = unit(c(2, 2, 10, 2), "mm"), annotation_legend_list=annotation)
    lapply(names(pdata_list), function(nn) { decorate_heatmap_body(nn, { grid.text(nn, y = unit(1, "npc") + unit(2, "mm"), just = "bottom") }) }) -> null
}

# Match the data matrix to the dimension of the template matrix and fill the blanks with NAs
match_matrices <- function(data_matrix, template_matrix) {
    template_matrix[] = NA
    if (!all(rownames(data_matrix) %in% rownames(template_matrix)) | !all(colnames(data_matrix) %in% colnames(template_matrix))) { print(data_matrix); print(template_matrix); stop("data_matrix not a subset of template in 'match_matrices'") }
    template_matrix[rownames(data_matrix), colnames(data_matrix)] = data_matrix
    return(template_matrix)
}

#' Heatmap with stars for extra data
#'
#' Plot heatmap with stars in the tile to encode for another data (significance, effect size, etc)
#' @param hm_values Values for the colors of the heatmap tiles
#' @param data Values for the text (or number of stars) of the heatmap tiles
#' @param encoding One of c('size', 'count', 'multi_count') for whether the value in data should be encoded by the size of one star or a number of stars on one or multiple lines
#' @param return_list If FALSE the heatmap is plotted with a legend, if TRUE returns a list with fields 'object' for the Heatmap object and 'annotation_legend_list' for the annotation corresponding to the encoding.
#' @export
#' @rdname starHeatmap
starHeatmap2 <- function(hm_values, data, encoding="size", return_list=FALSE, ...) {
    data = data %>% .[nrow(.):1,,drop=FALSE]
    if (encoding == 'size') {
        out = list(
            object=hm_values %>% .[nrow(.):1,,drop=FALSE] %>%
                Heatmap(cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(data[i,j]) && data[i,j] > 10) { grid.text("*", x, y-unit(0.2, "lines"), gp=gpar(fontsize=10*floor(data[i,j]/10))) }
                }, cluster_rows=FALSE, cluster_columns=FALSE, row_names_side="left", col=circlize::colorRamp2(c(0, 100), c("white", "red")), border=TRUE, column_title_side="bottom", ...),
            annotation_legend_list = list(Legend(title="Bliss score", at=10*1:5, labels=10*1:5, pch="*", legend_gp=gpar(fontsize=10*1:5), type="points", background="white", grid_height=unit(2, "lines")))
        )
    } else if (encoding == 'multi_count') {
        out = list(
            object = hm_values %>% .[nrow(.):1,,drop=FALSE] %>%
                Heatmap(cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(data[i,j]) && data[i,j] > 10) { grid.text(paste0(rep("*", floor(data[i,j]/10)), collapse=""), x, y-unit(0.2, "lines"), gp=gpar(fontsize=20)) }
                }, cluster_rows=FALSE, cluster_columns=FALSE, row_names_side="left", col=circlize::colorRamp2(c(0, 100), c("white", "red")), border=TRUE, column_title_side="bottom", ...),
                annotation_legend_list = list(Legend(title="Bliss score", at=10*1:5, labels=10*1:5, pch=1:5 %>% sapply(function(nn){paste0(rep("*", nn), collapse="")}), type="points"))
        )
    } else if (encoding == 'count') {
        out = list(
            object=hm_values %>% .[nrow(.):1,,drop=FALSE] %>%
                Heatmap(cell_fun=function(j, i, x, y, width, height, fill) {
                    if(!is.na(data[i,j]) && data[i,j] > 10) { grid.text(gsub("(.{3})", "\\1\n", paste0(rep("*", floor(data[i,j]/10)), collapse="")), x, y-unit(0.2, "lines"), gp=gpar(fontsize=20)) }
                }, cluster_rows=FALSE, cluster_columns=FALSE, row_names_side="left", col=circlize::colorRamp2(c(0, 100), c("white", "red")), border=TRUE, column_title_side="bottom", ...),
                annotation_legend_list = list(Legend(title="Bliss score", at=10*1:5, labels=10*1:5, pch=1:5 %>% sapply(function(nn){paste0(rep("*", nn), collapse="")}), type="points"))
            )
    }
    if (return_list) {
        return(out)
    } else {
        return(do.call(draw, out))
    }
}

#' Plot heatmap with values as stars
#'
#' Helper for starHeatmap2 with hm_values = data
#' @export
#' @rdname starHeatmap
starHeatmap <- function(data, ...) {
    starHeatmap2(data, data, ...)
}

#' Heatmap with values in the tiles
#'
#' @param hm_values Values for the heatmap colors
#' @param data Values for the heatmap text
#' @export
#' @rdname starHeatmap
valueHeatmap2 <- function(hm_values, data, ...) {
    data = data %>% .[nrow(.):-1:1,,drop=FALSE]
    hm_values %>% .[nrow(.):-1:1,,drop=FALSE] %>% Heatmap(cell_fun=function(j, i, x, y, width, height, fill) { if(!is.na(data[i,j]) && data[i,j] > 10) { grid.text(signif(data[i,j], 2), x, y, gp=gpar(fontsize=10)) } }, cluster_rows=FALSE, cluster_columns=FALSE, row_names_side="left", col=circlize::colorRamp2(c(0, 100), c("white", "red")), border=TRUE, column_title_side="bottom", ...)
}

