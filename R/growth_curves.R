############################## growth_curves.R #############################################
#' @import propagate
#' @import minpack.lm
#' @import tidyverse
#' @import cowplot
#' @import ggrepel

#' @title Process growth curves to viability values
#'
#' @description Extract growth rate, normalise and averages growth curves from the experiment
#' @param experiment Dataset as return by read_incu from incucyter
#' @param keep_controls Whether the controls should be present in the output. Control values are DMSO, medium and no_cells.
#' @param max_confluency The maximum value to consider valid. Used to remove data points when the growth tappers out because of contact inhibition. Adjust according to your data type.
#' @param space A function to modify the fitting space (log (default), identity, etc)
#' @return A tibble with new columns Inhibitor, Concentration and Viability
#' @export
# Note: this function is still buggy with keep_controls==TRUE because of concentration_values function applied to no conformant values extracted as "Inhibitor_Concentration"
process_growth_curve <- function(experiment, keep_controls=TRUE, max_confluency=75, space=log) {
    data_size = nrow(experiment)
    throwout = experiment %>% dplyr::filter(Value > max_confluency) %>% nrow
    if (throwout/data_size > 0.5) {
        warning("Throwing out >50% of the data because of the 'max_confluency' threshold. Make sure the threshold is compatible with your data type.")
    }
    if (length(experiment %>% ungroup %>% distinct(Elapsed)) < 2) { stop('Trying to fit a single time point in process_growth_curve, check your filters.') }
    well_agg = experiment %>% dplyr::filter(Value > 0) %>% # Remove drops due to lack of focus or other technical artefacts
        dplyr::select(Analysis_Job, Well, Treatment, Ref_T, Value, Elapsed) %>%
        dplyr::filter(Value < max_confluency) %>%
        mutate(Value=space(Value)) %>% # Linear in log space
        group_by(Analysis_Job, Well, Treatment, Ref_T) %>%
        group_map(function(.x, .y) { .y %>% mutate(Value=lm(Value~1+Elapsed, .x)$coefficients["Elapsed"], Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) }) %>%
        bind_rows %>%
        group_by(Analysis_Job, Treatment, Well, Ref_T) %>%
        summarise(Mean_well=mean(Value), Sd_well=sd(Value, na.rm=T)) %>% 
        ungroup()

    control_names = well_agg %>% distinct(Ref_T) %>% as_vector
    controls = well_agg %>%
            distinct(Ref_T) %>%
            as_vector %>% map(function(xx){ well_agg %>%
                              dplyr::filter(Treatment==xx) %>%
                              group_by(Treatment) %>%
                              summarise(Mean=mean(Mean_well), Sd=sd(Mean_well))
                    }) %>%
            bind_rows()
    if (nrow(controls) < length(control_names)) {
        not_found = control_names[!control_names %in% (controls %>% pull(Treatment))]
        warning(paste("The controls ", paste(not_found, collapse=","), " have not been found, can't normalize all conditions. Check the Treatment name for the controls and that all relevant wells are present"))
    }

    norm_values = controls %>% apply(1, function(xx) {well_agg %>% dplyr::filter(Ref_T==xx["Treatment"]) %>% mutate(Viability=Mean_well/as.numeric(xx["Mean"]), Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) }) %>% bind_rows() %>% mutate(Concentration_value=concentration_values(Concentration))
    
    if (keep_controls) {
        return(norm_values)
    } else {
        norm_values  %>% dplyr::filter(!grepl("DMSO|medium|cells", Treatment)) %>% return()
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
        growth_rate = log(trace %>% dplyr::filter(Elapsed==max(Elapsed)) %>% .$Value / trace %>% dplyr::filter(Elapsed==min(Elapsed)) %>% .$Value) * log(2) / (max(trace$Elapsed)-min(trace$Elapsed))
    } else {
        stop = 10 # Arbitrary value, should depend on the sampling interval (estimation of the time it takes the cell to gain 5% confluence
        while (mean(trace$Value[(stop-9):stop]) < max_confluency & stop < length(trace$Value)) { stop = stop + 1 }
        growth_rate = log(trace$Value[stop]/trace$Value[1], 2) / (trace$Elapsed[stop]-trace$Elapsed[1])
    }
    return( trace[1,] %>% mutate(Value=growth_rate, Elapsed=NULL, Date_Time=NULL, Feature="Growth rate", Unit="/h", Inhibitor=gsub("_.*", "", Treatment), Concentration=gsub(".*_", "", Treatment)) )
}

#' Fit sigmoid to drug sensitivities
#'
#' @export
fit_drug_sensitivity <- function(pexp, controls=c("DMSO", "control", "CTRL", "medium", "cell"), verbose=FALSE, color='WellX', shape='WellY') {
    pcontrol = pexp %>% dplyr::filter(grepl( paste0(controls, collapse="|"), Treatment ))
    pexp = pexp %>% dplyr::filter(!grepl( paste0(controls, collapse="|"), Treatment ))

    treatment_plots = list()
    fits_treatment = list()
    log_breaks = 10^-seq(1:10)
    log_breaks = c(log_breaks, 3*log_breaks)
    for (tt in distinct(pexp, Inhibitor)$Inhibitor) {
        t_data = pexp %>% dplyr::filter(Inhibitor == tt)
        t_control = pcontrol %>% dplyr::filter(Analysis_Job %in% pexp$Analysis_Job, Treatment %in% pexp$Ref_T)
        t_control = t_control %>% mutate( Concentration_value = apply(t_control, 1, function(XX){ median(pexp %>% dplyr::filter(Ref_T == XX["Treatment"]) %>% .$Concentration_value) }) )
        # Equidistant points in log space
        crange=log10(range(t_data$Concentration_value, na.rm=TRUE))
        if (verbose){ message(paste0(crange, collapse=",")) }
        if (is.na(crange[1])) {crange[1]=-7}
        if (is.na(crange[2])) {crange[2]=-4}
        xx = 10^(seq(crange[1], crange[2], length.out=100))
        dumx=tibble(xx=xx)

        fits_treatment[[tt]] = fit_sigmoid(t_data)
        t_data = tryCatch({
            suppressMessages({ ipred = predictNLS(fits_treatment[[tt]]$model, t_data[,"Concentration_value"])$summary })
            t_data %>% mutate(Vmin=ipred$"Prop.2.5%", Vmax=ipred$"Prop.97.5%")
        }, error = function(e){
            t_data %>% mutate(Vmin=t_data$Viability, Vmax=t_data$Viability)
        })
        t_data = t_data %>% mutate(Vmin = sapply(Vmin, max, -1, na.rm=TRUE), Vmax = sapply(Vmax, min, 1.5, na.rm=TRUE)) # Prevent the confidence interval from blowing the limits of the y-axis
        # Select the features to use for annotation
        choose_feature_content <- function(feature='color', content_value) {
            if (content_value == 'WellX') { t_data = t_data %>% mutate(content = substr(Well, 1, 1)) }
            else if (content_value == 'WellY') { t_data = t_data %>% mutate(content = substr(Well, 2, 2)) }
            else if (content_value %in% colnames(t_data)) { t_data = t_data %>% mutate(content = get(content_value)) }
            else if (content_value == '') { t_data = t_data %>% mutate(content='') }
            else { stop(paste0('Invalid ', feature, ' : ', content_value)) }

            if (feature=='shape') { t_data = t_data %>% dplyr::rename(shape=content) }
            else { t_data = t_data %>% dplyr::rename(color=content) }
            return(t_data)
        }
        t_data = choose_feature_content('color', color)
        t_data = choose_feature_content('shape', shape)
        # Generate dose response curve plots
        new_plot = t_data %>% ggplot() + scale_x_log10(breaks=log_breaks) +
                geom_point(aes( Concentration_value, Viability, col=color, shape=shape )) +
                geom_ribbon(aes(Concentration_value, ymin=Vmin, ymax=Vmax), fill="grey", alpha=0.5) +
                geom_point(aes(Concentration_value, Viability, color="controls", alpha=0.5), show.legend=FALSE, data=t_control) +
                geom_line(aes(xx, sigmoid(fits_treatment[[tt]]$coef, log10(xx))), data=dumx) +
                coord_cartesian(ylim=c(-0.5, 1.5)) +
                ggtitle(tt)
        if (nrow(t_data %>% dplyr::filter(Viability < -0.5)) > 0) {
            new_plot = new_plot + geom_text_repel(aes(Concentration_value, -0.5, label=signif(Viability, 2)), data=.%>%dplyr::filter(Viability < -0.5))
        }
        if (nrow(t_data %>% dplyr::filter(Viability > 1.5)) > 0) {
            new_plot = new_plot + geom_text_repel(aes(Concentration_value, 1.5, label=signif(Viability, 2)), data=.%>%dplyr::filter(Viability > 1.5))
        }
        print(new_plot)
        treatment_plots[[tt]] = new_plot
    }
    control_plot = pcontrol %>% dplyr::filter(grepl(paste0("DMSO"), Treatment)) %>% mutate(Concentration_value=1/Concentration_value) %>% ggplot(aes(Concentration_value, Mean_well, color=Analysis_Job)) + geom_point() + scale_x_log10(breaks=log_breaks) + xlab("Carrier dilution") + ggtitle("Controls growth rates")
    print(control_plot)
    return(list(fits=fits_treatment, plots=treatment_plots, controls=control_plot))
}

#' Plot the drug sensitivity with the fit
#'
#' @param fit_result Result of fit_drug_sensitivity
#' @param drug_name Drug for which the plot should be made
#' @export
plot_ic50_fit <- function(fit_result, drug_name='') {
    if (!(drug_name %in% names(fit_result$plots) | drug_name %in% names(fit_result$fits))) {
        stop(paste('Drug', drug_name, 'does not exist in this dataset'))
    }
    log_breaks = 10^-seq(1:10)
    log_breaks = c(log_breaks, 3*log_breaks)
    crange=log10(range(fit_result$plots[[drug_name]]$data$Concentration_value, na.rm=TRUE))
    if (is.na(crange[1])) {crange[1]=-7}
    if (is.na(crange[2])) {crange[2]=-4}
    xx = 10^(seq(crange[1], crange[2], length.out=100))
    dumx=tibble(xx=xx)
    result_plot = fit_result$plots[[drug_name]]$data %>% ggplot() + scale_x_log10(breaks=log_breaks) +
                geom_point(aes( Concentration_value, Viability, col=substr(Well, 1, 1), shape=substr(Well, 2, 2) )) +
                geom_ribbon(aes(Concentration_value, ymin=Vmin, ymax=Vmax), fill="grey", alpha=0.5) +
                #geom_point(aes(Concentration_value, Viability, color="controls", alpha=0.5), show.legend=FALSE, data=t_control) +
                geom_line(aes(xx, sigmoid(fit_result$fits[[drug_name]]$coef, log10(xx))), data=dumx) +
                coord_cartesian(ylim=c(-0.5, 1.5)) +
                ggtitle(drug_name)
    print(result_plot)
}

#' Compute the slope of the confluency in a window
#'
#' Compute the slope of the confluency in a window using the max-min approximation. Does the direct sum, considering the values are log-transformed.
#' @param experiment Dataset as return by read_incu from incucyter
#' @param winsize The half-size of the window in time units
#' @param do_plot Boolean whether the slopes microplate should be plotted
#' @return Invisibly the experiment tibble without the first 'winsize' data points and an extra column 'Slope'.
#' @export
compute_window_slope <- function(experiment, winsize=20, do_plot=TRUE) {
    experiment %>%
        dplyr::filter(Elapsed>winsize) %>%
        group_by(Well, img, Elapsed) %>%
        mutate(Slope=NA) %>%
        group_map(function(.x, .y){
            time = as.numeric(.y$Elapsed)
            experiment %>% dplyr::filter(Elapsed > time-winsize & Elapsed < time+winsize) -> tmp
            max_time = max(tmp$Elapsed)
            min_time = min(tmp$Elapsed)
            .x$Slope = ((experiment %>% dplyr::filter(img==.y$img, Well==.y$Well, Elapsed==max_time) %>% pull(Value))-(experiment %>% dplyr::filter(img==.y$img, Well==.y$Well, Elapsed==min_time) %>% pull(Value)))/(max_time-min_time)
            return(bind_cols(.y, .x)) }) %>%
        bind_rows -> with_slope
    if (do_plot) {
        plot(with_slope %>% mutate(Value=Slope) %>% plot_microplate(plot_name=paste("Windowed slope (", winsize, ")", sep=" ")))
    }

    invisible(with_slope)
}

#' Plot processed growth rates
#'
#' Plot a heatmap of the mean growth rate in each well of a plate.
#' @param pgc A tibble as outputed by process_growth_curve
#' @export
plot_growth_rates <- function(pgc) {
     pgc %>% separate(Well, into=c("Y", "X"), sep=1) %>%
         ggplot(aes(as.factor(as.numeric(X)), reorder(Y,desc(Y)), fill=Viability, label=signif(Viability, 3))) + geom_tile() + geom_text() + facet_wrap(Elapsed~Analysis_Job, switch="y") + xlab('WellX') + ylab('WellY')
}

#' Convert a well ID into X and Y coordinates
#' 
#' Convert a well ID into X and Y coordinates that are factors and probably ordered for facet_wrap
#' @param input_tibble A tibble with column 'Well' matching the regex '[A-Z]\\d\\+' (one letter and a number)
#' export
microplate <- function(input_tibble) {
    input_tibble %>% ungroup %>% separate(Well, into=c("WellX", "WellY"), sep=1, remove=FALSE) %>% mutate(WellY=as.factor(as.numeric(WellY)), WellX=reorder(WellX,desc(WellX)))
}

