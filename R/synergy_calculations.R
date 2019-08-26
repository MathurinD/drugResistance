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
        mutate(response = 100*Viability) # Convert Viability in ~[0,1] to a reponse in ~[0,100]
    if (restrict) { # Scale the response strictly between 0 and 100
        pdata = pdata %>%
            mutate(response = case_when(response<0~0, response>100~100, TRUE~response))
    }

    # Convert the treatments to 2 columns with the names expected by synergyfinder
    treatments = str_split_fixed(pdata$Treatment, "\\+", 2) %>% as_tibble %>% rename(Drug1=V1, Drug2=V2) %>% separate(Drug1, c("drug_row", "conc_r"), "_") %>% separate(Drug2, c("drug_col", "conc_c"), "_")
    drugs_col = treatments %>% distinct(drug_col) %>% pull(drug_col)
    drugs_row = treatments %>% distinct(drug_row) %>% pull(drug_row)
    row_drug = find_first_drug(drugs_row, "DMSO")
    col_drug = find_first_drug(drugs_col, c(row_drug, "DMSO"))
    treatments %>% apply(1, function(x){
                            if(x["drug_row"] == "DMSO") {
                                 x[c("drug_row", "conc_r")]=c("", NA)
                            } else if(x["drug_row"] %in% c(col_drug)) {
                                x[c("drug_col", "conc_c")]=x[c("drug_row", "conc_r")]
                                x[c("drug_row", "conc_r")]=c("", NA)
                            }
                            return(x) }) %>% t %>% as.tibble %>% mutate(drug_row=row_drug, drug_col=col_drug) %>%
         mutate(conc_r=case_when(is.na(conc_r)~"0nM", TRUE~conc_r), conc_c=case_when(is.na(conc_c)~"0nM", TRUE~conc_c)) %>%
         mutate(conc_r = concentration_values(conc_r)*1e9, conc_r_unit = "nM", conc_c = concentration_values(conc_c)*1e9, conc_c_unit = "nM") -> treatments

    synergy_data = pdata %>%
        bind_cols(treatments) %>% # Add the treatment columns
        group_by(drug_col, conc_c, drug_row, conc_r) %>%
        summarize(response=mean(response)) %>% mutate(block_id=1, conc_c_unit="nM", conc_r_unit="nM") %>% # Make sure each concentration combination corresponds to 1 data point (not done earlier because controls can have different DMSO concentration thus different names but all correspond to the same drug combination of 0)
        mutate(block_id=1) %>%
        ungroup() %>%
        select(block_id, drug_col, conc_c, conc_c_unit, drug_row, conc_r, conc_r_unit, response) # Sanity check that all column expected by synergyfinder are present

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
