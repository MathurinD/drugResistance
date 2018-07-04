############################## growth_curves.R #############################################
#' @import propagate
#' @import tidyverse
#' @import cowplot

#' @title Process growth curves to viability values
#'
#' Extract growth rate, normalise and averages growth curves from the experiment
#' @param experiment Dataset as return by read_incu from incucyter
#' @param keep_controls Whether the controls should be present in the output. Control values are DMSO, medium and no_cells.
#' @return A tibble with new columns Inhibitor, Concentration and Viability
#' @export
# Note: this function is still buggy with keep_controls==TRUE because of concentration_values function applied to no conformant values extracted as "Inhibitor_Concentration"
process_growth_curve <- function(experiment, keep_controls=TRUE) {
    well_agg = experiment %>% filter(Value > 1) %>% # Remove drops due to lack of focus or other technical artefacts
        distinct(Analysis_Job, Well, img) %>%
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

    norm_values = controls %>% apply(1, function(xx) {well_agg %>% filter(Ref_T==xx["Treatment"]) %>% mutate(Viability=Mean_well/as.numeric(xx["Mean"]), Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) }) %>% bind_rows() %>% mutate(Concentration_value=concentration_values(Concentration))
    
    if (keep_controls) {
        return(norm_values)
    } else {
        norm_values  %>% filter(!grepl("DMSO|medium|cells", Treatment)) %>% return()
    }
}

#' Process one growth curve to a viability value
#'
#' Process one growth curve to a viability value
#' @param trace A tibble with columns Elapsed for the time and Value for the confluency at those times
#' @param max_confluency Value above which the confluency is not growing exponentially anymore
#' @return A one line tibble with the same columns where value is now the growth rate
#' @export
process_single_growth_curve <- function(trace, max_confluency=75) {
    if (max(trace$Value) < max_confluency) {
        growth_rate = log(trace %>% filter(Elapsed==max(Elapsed)) %>% .$Value / trace %>% filter(Elapsed==min(Elapsed)) %>% .$Value) * log(2) / (max(trace$Elapsed)-min(trace$Elapsed))
    } else {
        stop = 10 # Arbitrary value, should depend on the sampling interval (estimation of the time it takes the cell to gain 5% confluence
        while (mean(trace$Value[(stop-9):stop]) < max_confluency & stop < length(trace$Value)) { stop = stop + 1 }
        growth_rate = log(trace$Value[stop]/trace$Value[1]) * log(2) / (trace$Elapsed[stop]-trace$Elapsed[1])
    }
    return( trace[1,] %>% mutate(Value=growth_rate, Elapsed=NULL, Date_Time=NULL, Feature="Confluence growth rate", Unit="%/h", Description=paste0(Description, " confluence growth rate"), Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) )
}

#' Fit sigmoid to drug sensitivities
#'
#' @export
fit_drug_sensitivity <- function(pexp, controls=c("DMSO", "control", "medium", "cell")) {
    pcontrol = pexp %>% filter(grepl( paste0(controls, collapse="|"), Treatment ))
    pexp = pexp %>% filter(!grepl( paste0(controls, collapse="|"), Treatment ))

    treatment_plots = list()
    fits_treatment = list()
    log_breaks = 10^-seq(1:10)
    log_breaks = c(log_breaks, 3*log_breaks)
    for (tt in distinct(pexp, Inhibitor)$Inhibitor) {
        t_data = pexp %>% filter(Inhibitor == tt)
        t_control = pcontrol %>% filter(Analysis_Job %in% pexp$Analysis_Job, Treatment %in% pexp$Ref_T)
        t_control = t_control %>% mutate( Concentration_value = apply(t_control, 1, function(XX){ median(pexp %>% filter(Ref_T == XX["Treatment"]) %>% .$Concentration_value) }) )
        # Equidistant points in log space
        crange=log(range(t_data$Concentration_value))
        xx = exp(seq(crange[1], crange[2], length.out=100))
        dumx=tibble(xx=xx)

        fits_treatment[[tt]] = fit_sigmoid(t_data)
        t_data = tryCatch({
            suppressMessages({ ipred = predictNLS(fits_treatment[[tt]]$model, t_data[,"Concentration_value"])$summary })
            t_data %>% mutate(Vmin=ipred$"Prop.2.5%", Vmax=ipred$"Prop.97.5%")
        }, error = function(e){
            t_data %>% mutate(Vmin=t_data$Viability, Vmax=t_data$Viability)
        })
        t_data = t_data %>% mutate(Vmin = sapply(Vmin, max, -1, na.rm=TRUE), Vmax = sapply(Vmax, min, 1.5, na.rm=TRUE)) # Prevent the confidence interval from blowing the limits of the y-axis
        new_plot = t_data %>% ggplot() + scale_x_log10(breaks=log_breaks) +
                geom_point(aes( Concentration_value, Viability, col=substr(Well, 1, 1), shape=substr(Well, 2, 2) )) +
                geom_ribbon(aes(Concentration_value, ymin=Vmin, ymax=Vmax), fill="grey", alpha=0.5) +
                geom_point(aes(Concentration_value, Viability, color="controls", alpha=0.5), show.legend=FALSE, data=t_control) +
                geom_line(aes(xx, sigmoid(fits_treatment[[tt]]$coef, log10(xx))), data=dumx) +
                ggtitle(tt)
        print(new_plot)
        treatment_plots[[tt]] = new_plot
    }
    control_plot = pcontrol %>% filter(grepl(paste0("DMSO"), Treatment)) %>% mutate(Concentration_value=1/Concentration_value) %>% ggplot(aes(Concentration_value, Mean_well, color=Analysis_Job)) + geom_point() + scale_x_log10(breaks=log_breaks) + xlab("Carrier dilution") + ggtitle("Controls growth rates")
    print(control_plot)
    return(list(fits=fits_treatment, plots=treatment_plots, controls=control_plot))
}
