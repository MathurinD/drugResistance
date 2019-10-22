############################## growth_curves.R #############################################
#' @import propagate
#' @import minpack.lm
#' @import tidyverse
#' @import cowplot
#' @import ggrepel

#' @title Process growth curves to viability values
#'
#' Extract growth rate, normalise and averages growth curves from the experiment
#' @param experiment Dataset as return by read_incu from incucyter
#' @param keep_controls Whether the controls should be present in the output. Control values are DMSO, medium and no_cells.
#' @return A tibble with new columns Inhibitor, Concentration and Viability
#' @export
# Note: this function is still buggy with keep_controls==TRUE because of concentration_values function applied to no conformant values extracted as "Inhibitor_Concentration"
process_growth_curve <- function(experiment, keep_controls=TRUE, max_confluency=75) {
    well_agg = experiment %>% filter(Value > 0) %>% # Remove drops due to lack of focus or other technical artefacts
        filter(Value < max_confluency) %>%
        mutate(Value=log(Value)) %>% # Linear in log space
        group_by(Analysis_Job, Well, img, Description, Treatment, Ref_T, Reference, Metric) %>%
        group_map(function(.x, .y) { .y %>% mutate(Value=lm(Value~1+Elapsed, .x)$coefficients["Elapsed"], Feature="Confluence growth rate", Unit="%/h", Description=paste0(Description, " confluence growth rate"), Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) }) %>%
       # distinct(Analysis_Job, Well, img) %>%
       # apply(1, function(x){ experiment %>%
       #                     filter(Analysis_Job==x["Analysis_Job"], Well==x["Well"], img==x["img"]) %>%
       #                     process_single_growth_curve(max_confluency)
       #                 }) %>%
        bind_rows %>%
        group_by(Analysis_Job, Treatment, Reference, Well, Feature, Unit, Metric, Ref_T) %>%
        summarise(Mean_well=mean(Value), Sd_well=sd(Value, na.rm=T)) %>% 
        ungroup()

    control_names = well_agg %>% distinct(Ref_T) %>% as_vector
    controls = well_agg %>%
            distinct(Ref_T) %>%
            as_vector %>% map(function(xx){ well_agg %>%
                              filter(Treatment==xx) %>%
                              group_by(Treatment) %>%
                              summarise(Mean=mean(Mean_well), Sd=sd(Mean_well))
                    }) %>%
            bind_rows()
    if (nrow(controls) < length(control_names)) {
        not_found = control_names[!control_names %in% (controls %>% pull(Treatment))]
        warning(paste("The controls ", paste(not_found, collapse=","), " have not been found, can't normalize all conditions. Check the Treatment name for the controls and that all relevant wells are present"))
    }

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
        growth_rate = log(trace$Value[stop]/trace$Value[1], 2) / (trace$Elapsed[stop]-trace$Elapsed[1])
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
        message(crange)
        if (is.na(crange[1])) {crange[1]=-17}
        if (is.na(crange[2])) {crange[2]=-9}
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
                #geom_line(aes(xx, fit$curve[[1]](xx)), data=dumx) + # If fit=drc::drm(Viability~Concentration_value, data=t_data, fct=drc::LL2.2(names=c("Slope", "log(IC50)")))
                #geom_line(aes(xx, sigmoid(c(-fit$coefficients["Slope:(Intercept)"], fit$coefficients["log(IC50):(Intercept)"] %>% exp %>% log10), log10(xx))), data=dumx) + # If fit=drc::drm(Viability~Concentration_value, data=t_data, fct=drc::LL2.2(names=c("Slope", "log(IC50)")))
                #geom_line(aes(xx, sigmoid(c(-confint(fit)["Slope:(Intercept)", "2.5 %"], confint(fit)["log(IC50):(Intercept)", "2.5 %"] %>% exp %>% log10), log10(xx))), data=dumx, col="grey") + geom_line(aes(xx, sigmoid(c(-confint(fit)["Slope:(Intercept)", "97.5 %"], confint(fit)["log(IC50):(Intercept)", "97.5 %"] %>% exp %>% log10), log10(xx))), data=dumx, col="grey") # Confidence intervals
                coord_cartesian(ylim=c(-0.5, 1.5)) +
                ggtitle(tt)
        if (nrow(t_data %>% filter(Viability < -0.5)) > 0) {
            new_plot = new_plot + geom_text_repel(aes(Concentration_value, -0.5, label=signif(Viability, 2)), data=.%>%filter(Viability < -0.5))
        }
        if (nrow(t_data %>% filter(Viability > 1.5)) > 0) {
            new_plot = new_plot + geom_text_repel(aes(Concentration_value, 1.5, label=signif(Viability, 2)), data=.%>%filter(Viability > 1.5))
        }
        print(new_plot)
        treatment_plots[[tt]] = new_plot
    }
    control_plot = pcontrol %>% filter(grepl(paste0("DMSO"), Treatment)) %>% mutate(Concentration_value=1/Concentration_value) %>% ggplot(aes(Concentration_value, Mean_well, color=Analysis_Job)) + geom_point() + scale_x_log10(breaks=log_breaks) + xlab("Carrier dilution") + ggtitle("Controls growth rates")
    print(control_plot)
    return(list(fits=fits_treatment, plots=treatment_plots, controls=control_plot))
}

#' Compute the slope of the confluency in a window
#'
#' Compute the slope of the confluency in a window using the max-min approximation. Does the direct sum, considering the values are log-transformed
#' @param experiment An experiment tibble as returned by process_growth_curve
#' @param winsize The half-size of the window in time units
#' @param do_plot Boolean whether the slopes microplate should be plotted
#' @return Invisibly the experiment tibble without the first 'winsize' data points and an extra column 'Slope'.
#' @export
compute_window_slope <- function(experiment, winsize=20, do_plot=TRUE) {
    experiment %>%
        filter(Elapsed>winsize) %>%
        group_by(Well, img, Elapsed) %>%
        mutate(Slope=NA) %>%
        group_map(function(.x, .y){
            time = as.numeric(.y$Elapsed)
            experiment %>% filter(Elapsed > time-winsize & Elapsed < time+winsize) -> tmp
            max_time = max(tmp$Elapsed)
            min_time = min(tmp$Elapsed)
            .x$Slope = ((experiment %>% filter(img==.y$img, Well==.y$Well, Elapsed==max_time) %>% pull(Value))-(experiment %>% filter(img==.y$img, Well==.y$Well, Elapsed==min_time) %>% pull(Value)))/(max_time-min_time)
            return(bind_cols(.y, .x)) }) %>%
        bind_rows -> with_slope
    if (do_plot) {
        plot(with_slope %>% mutate(Value=Slope) %>% plot_microplate(plot_name=paste("Windowed slope (", winsize, ")", sep=" ")))
    }

    invisible(with_slope)
}

