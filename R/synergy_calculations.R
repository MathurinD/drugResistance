#' @import synergyfinder

#' @title Get a valid input for ReshapeData
#'
#' Get a synergy table that can be used by the function ReshapeData synergyfinder.
#' @param pdata A tibble as outputed from process_growth_curves. Must have columns 'Treatment' with format drugRow_rowConcentration+drugCol_colConcentation and 'Viability' 
#' @param scale Whether the data should be strictly restricted between 0 and 100. Apply a cutoff, not a scaling.
#' @return A tibble with column names as required by ReshapeData
#' @export
get_synergy_table <- function(pdata, restrict=FALSE) {
    #pigfrazd %>% filter(!grepl("1|8", Well)) %>%
    pdata = pdata %>%
        mutate(Viability=case_when(is.infinite(Viability)&Viability<0~0, is.infinite(Viability)&Viability>0~1, TRUE~Viability)) %>% # Remove infinite
        mutate(Response = 100*Viability) # Convert Viability in ~[0,1] to a reponse in ~[0,100]
    if (restrict) { # Scale the Response strictly between 0 and 100
        pdata = pdata %>%
            mutate(Response = case_when(Response<0~0, Response>100~100, TRUE~Response))
    }

    # Convert the treatments to 2 columns with the names expected by synergyfinder
    treatments = str_split_fixed(pdata$Treatment, "\\+", 2) %>% as_tibble %>% rename(Drug1=V1, Drug2=V2) %>% separate(Drug1, c("DrugRow", "ConcRow"), "_") %>% separate(Drug2, c("DrugCol", "ConcCol"), "_")
    drugs_col = treatments %>% distinct(DrugCol) %>% pull(DrugCol)
    drugs_row = treatments %>% distinct(DrugRow) %>% pull(DrugRow)
    row_drug = find_first_drug(drugs_row, "DMSO")
    col_drug = find_first_drug(drugs_col, c(row_drug, "DMSO"))
    treatments %>% apply(1, function(x){
                            if(x["DrugRow"] == "DMSO") {
                                 x[c("DrugRow", "ConcRow")]=c("", NA)
                            } else if(x["DrugRow"] %in% c(col_drug)) {
                                x[c("DrugCol", "ConcCol")]=x[c("DrugRow", "ConcRow")]
                                x[c("DrugRow", "ConcRow")]=c("", NA)
                            }
                            return(x) }) %>% t %>% as.tibble %>% mutate(DrugRow=row_drug, DrugCol=col_drug) %>%
         mutate(ConcRow=case_when(is.na(ConcRow)~"0nM", TRUE~ConcRow), ConcCol=case_when(is.na(ConcCol)~"0nM", TRUE~ConcCol)) %>%
         mutate(ConcRow = concentration_values(ConcRow)*1e9, ConcRowUnit = "nM", ConcCol = concentration_values(ConcCol)*1e9, ConcColUnit = "nM") -> treatments

    synergy_data = pdata %>%
        bind_cols(treatments) %>% # Add the treatment columns
        group_by(DrugCol, ConcCol, DrugRow, ConcRow) %>%
        summarize(Response=mean(Response)) %>% mutate(BlockId=1, ConcColUnit="nM", ConcRowUnit="nM") %>% # Make sure each concentration combination corresponds to 1 data point (not done earlier because controls can have different DMSO concentration thus different names but all correspond to the same drug combination of 0)
        mutate(BlockID=1) %>%
        ungroup() %>%
        select(BlockID, DrugCol, ConcCol, ConcColUnit, DrugRow, ConcRow, ConcRowUnit, Response) # Sanity check that all column expected by synergyfinder are present

    return(synergy_data)
}

find_first_drug <- function(drugs_list, exclude=c("DMSO")) {
    for (cc in drugs_list) {
        if (! cc %in% c("", exclude)) {
            return(cc)
        }
    }
    return(NA)
}
