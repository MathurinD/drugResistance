# drugResistance
Helper functions to compute drug sensitivity and synergy from incucyte data.

(optional) Load the data using [incucyter](https://github.com/MathurinD/incucyter) and check the growth curves.

```
library(incucyter)
library(drugResistance)
confluency = read_incu(system.file('extdata', 'data.txt', package='drugResistance'), system.file('extdata', 'annotation.csv', package='drugResistance'))
plot_microplate(confluency)
```

The following columns are needed to compute a growth rate:
+ Analysis_job: name of the run/experiment
+ Well: The position of the well as a letter and a number (e.g D5)
+ Treatment: drug_concentration string separated by a "+" for combinations (eg. AZD6244_10nM+Panobinostat_1nM)
+ Ref_T: name of treatment used in the control well(s) (drug + conc, e,g DMSO_0.001)
+ Value: Confluence growth rate
+ Elapsed: Treatment time (unit is free but should be the same in all the column)

The confluency file from read_incu also has the following columns that are not used in the subsequent functions:
- Feature: Type of measurement (e.g. Confluency)
- Unit: of the feature (e.g. %)
- Metric: Type of mesurement
- Reference: position of the control well(s)

Compute growth rates and fit dose-response curves.
Growth rates are normalized for each Analysis_job to the average growth of the control wells (specified for each well in column Ref_T)

```
growth = process_growth_curve(confluency)
plot_growth_rates(growth) # Heatmap of the computed average growth rate in each well
dose_response = fit_drug_sensitivity(growth)
```

process_growth_curve will add the following columns:
- Mean_Well: mean across measurements if the well has several positions (or the only value measured)
- Sd_Well: Standard deviation across measurements if the well has several positions (or zero)
- Viability: normalised data
- Inhibitor: treatement wo/ conc
- Concentration: concentration as a string (e.g 20nM)
- Concentration_value: concentration in micromolar (e.g. 0.02, you can use the helper function `concentration_values` to convert from the Concentration column)

Visualise the dose-response curves and export the fitted values.
```
dose_response$plots
ic50s = data.frame(t(sapply(dose_response, function(xx){c(xx$coef, xx$residual)})))
names(ic50s) = c('slope', 'ic50', 'residual')
write.table(ic50s, system.file('exdata', 'ic50s.csv', package='drugResistance'))
```

You can also compute and visualise synergy scores (note test data are not synergies, see [Dorel2021][https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009515] for example data and code)
```
growth %>% plot_bliss_scores(col_drug="drug1") # Uses the helper function starHeatmap2 which uses ComplexHeatmap to plot growth rate as tile fill color and synergy score as stars
```
## Synergies

If your drug screen has combinations, you can compute synergy:
```
synergy = get_synergy_table(growth, control='DMSO') # Computes multiple synergy statistics
plot_bliss_scores(growth) # Plots only bliss synergy score
write_csv('synergy.csv', get_synergy_table(growth)) # Can be used as input for https://synergyfinder.fimm.fi/ with minor renaming of the columns
```
