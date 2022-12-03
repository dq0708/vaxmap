# Modeling and presentation of vaccination coverage estimates

Here are some example R code for the analyses and visualizations in "Modeling and presentation of vaccination coverage estimates using data from household surveys" (https://arxiv.org/abs/2004.03127).

- *data_cleaning.R*: code for data cleaning.
- *model_fit.R*: code for fitting varous models using INLA in R. 
- *model_prediction.R*: code for obtaining posterior prediction samples.
- *plot_prediction.R*: code for plotting maps of posterior medians and CI width. 
- *plot_ridgeline.R*: code for creating ridgeline plots of posterior distributions.
- *plot_rank.R*: code for creating histograms of posterior ranking distributions.
- *plot_exceedance.R*: code for creating maps of posterior exceedence probabilities.
- *plot_hatch.R*: code for creating state maps with hatching based on coefficient of variation (CV) of vaccination odds.

- *model_TCP_quantile.R*: code for obtaining quantile intervals, color assignment and true classfication rate (TCP) based on posterior samples. 
- *plot_TCP_quantile.R*: code for plotting maps and histograms with discrete color scales using quantile intervals.

- *model_TCP_threshold.R*:	code for obtaining color assignment and true classfication rate (TCP) given a set of pre-specified intervals based on posterior samples. 
- *plot_TCP_threshold.R*:	code for plotting maps and histograms with discrete color scales using a set of pre-specified intervals.
