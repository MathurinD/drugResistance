#' Process growth curves to viability values
#'
#' Extract growth rate, normalise and averages growth curves from the experiment
#' @param experiment Dataset as return by read_incu from incucyter
#' @param keep_controls Whether the controls should be present in the output. Control values are DMSO, medium and no_cells.
#' @return A tibble with new columns Inhibitor, Concentration and Viability
#' @export
# Note: this function is still buggy with keep_controls==TRUE because of concentration_values function applied to no conformant values extracted as "Inhibitor_Concentration"
process_growth_curve <- function(experiment, keep_controls=FALSE) {
    well_agg = experiment %>% distinct(Analysis_Job, Well, img) %>%
        apply(1, function(x){ experiment %>%
                            filter(Analysis_Job==x["Analysis_Job"], Well==x["Well"], img==x["img"]) %>%
                            process_single_growth_curve()
                        }) %>%
        bind_rows %>%
        group_by(Analysis_Job, Treatment, Reference, Well, Feature, Unit, Metric, Ref_T) %>%
        summarise(Mean_well=mean(Value), Sd_well=sd(Value, na.rm=T)) %>% 
        ungroup()

    controls = well_agg %>%
            distinct(Ref_T) %>%
            as_vector %>% map(function(xx){ well_agg %>%
                              filter(Treatment==xx) %>%
                              group_by(Treatment) %>%
                              summarise(Mean=mean(Mean_well), Sd=sd(Mean_well))
                    }) %>%
            bind_rows()

    norm_values = controls %>% apply(1, function(xx) {well_agg %>% filter(Ref_T==xx["Treatment"]) %>% mutate(Viability=Mean_well/as.numeric(xx["Mean"]), Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) }) %>% bind_rows() %>% mutate(Concentration_value=concentration_values(Concentration)) %>% return()
}

#' Process one growth curve to a viability value
#'
#' Process one growth curve to a viability value
#' @param trace A tibble with columns Elapsed for the time and Value for the confluency at those times
#' @return A one line tibble with the same columns where value is now the growth rate
#' @export
process_single_growth_curve <- function(trace) {
    if (max(trace$Value) < 95) {
        growth_rate = (trace %>% filter(Elapsed==max(Elapsed)) %>% .$Value - trace %>% filter(Elapsed==min(Elapsed)) %>% .$Value) / (max(trace$Elapsed)-min(trace$Elapsed))
    } else {
        stop = 10 # Arbitrary value, should depend on the sampling interval (estimation of the time it takes the cell to gain 5% confluence
        while (mean(trace$Value[(stop-9):stop]) < 95 & stop < length(trace$Value)) { stop = stop + 1 }
        growth_rate = (trace$Value[stop]-trace$Value[1]) / (trace$Elapsed[stop]-trace$Elapsed[1])
    }
    return( trace[1,] %>% mutate(Value=growth_rate, Elapsed=NULL, Date_Time=NULL, Feature="Confluence growth rate", Unit="%/h", Description=paste0(Description, " confluence growth rate"), Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) )
}
