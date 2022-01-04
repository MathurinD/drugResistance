# drugResistance
Helper functions to compute drug sensitivity and synergy from incucyte data.

Load the data using [incucyter](https://github.com/MathurinD/incucyter) and check the growth curves.

```
library(incucyter)
library(drugResistance)
confluency = read_incu(system.file('extdata', 'data.txt', package='drugResistance'), system.file('extdata', 'annotation.csv', package='drugResistance'))
plot_microplate(confluency)
```

Compute growth rates and fit dose-response curves.

```
growth = process_growth_curve(confluency)
dose_response = fit_drug_sensitivity(growth)
```

Visualise the dose-response curves and export the fitted values.
```
dose_response$plots
ic50s = data.frame(t(sapply(dose_response, function(xx){c(xx$coef, xx$residual)})))
names(ic50s) = c('slope', 'ic50', 'residual')
write.table(ic50s, system.file('exdata', 'ic50s.csv', package='drugResistance'))
```


