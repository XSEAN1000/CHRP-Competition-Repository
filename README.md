# COVID-19 Spread in the U.S. and Temperature: a Statistical Analysis

Thanks for clicking on my project!


**Please see the detailed [summary](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/Presentation1.pdf) of my project, typed and explained neatly for you, with LaTeX.**


### Brief Summary:

Understanding the relationship between disease spread and temperature is important for effectively distributing resources across the country, and for predicting how the disease will effect us over time. We perform a detailed analysis of the correlation between COVID-19 rate of spread in the U.S., and local temperature. Our analysis consists of: 

1. **Estimating the COVID spread rate in each of several communities.** We do this using a first order linear difference approximation (_N'(t)_) for the case totals (i.e. the discrete first derivative of the case totals in a given community); we also use the daily multiplicative factors (_MF(t)_) (i.e. by what factor did the case total multiply from day i to day i+1). The linear difference approximation is a good estimate of the derivative of the Piecewise Exponential approximation to the case total, [pictured here](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/PWE_NY_Apr12.pdf) for NYC, April 12:

2. **Several temperature variables.** Our data from the NOAA consists of daily temperature maximums, minimums, and averages, for several communities. We account for the **variable incubation relate/testing delay** by offsetting temperature and rate data, from 0 days to 2 weeks.

3. **Multiple testing** We have two rate measurements, three temperature variables, 15 possible temperature offsets, and three correlation test statistics (Pearson, Kendall, and Spearman). **We perform one correlation test for each possible combination of rate, temperature, offset, and correlation statistic.** For each choice of rate variable, temperature variable, and correlation statistic, we regard the 15 tests (for 15 offsets) as a multiple test; we evaluate each multiple test based on the standard rule p < 0.05 and on the Benjamini-Hochberg procedure for controlling false discovery rate (q = 0.05). Here is a series of graphs summarizing the test results for:
   * [_N'(t)_ and Spearman's test]
   * [_N'(t)_ and Pearson's test]
   * [_N'(t)_ and Kendall's test]
   * [_MF(t)_ and Spearman's test]
   * [_MF(t)_ and Pearson's test]
   * [_MF(t)_ and Kendall's test]
   
4. **Regression**


Again, thanks for clicking, and hope you enjoy! Constructive criticism and questions are very welcome.




