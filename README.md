# COVID-19 Spread in the U.S. and Temperature: a Statistical Analysis

Thanks for clicking on my project!


**Please see the detailed [summary](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/Presentation1.pdf) of my project, typed and explained neatly for you, with LaTeX.**


### Brief Summary:

Understanding the relationship between disease spread and temperature is important for effectively distributing resources across the country, and for predicting how the disease will effect us over time. We perform a detailed analysis of the correlation between COVID-19 rate of spread in the U.S., and local temperature. Our analysis consists of: 

1. **Estimating the COVID spread rate in each of several communities.** We do this using a first order linear difference approximation (_N'(t)_) for the case totals (i.e. the discrete first derivative of the case totals in a given community); we also use the daily multiplicative factors (_MF(t)_) (i.e. by what factor did the case total multiply from day i to day i+1). The linear difference approximation is a good estimate of the derivative of the Piecewise Exponential approximation to the case total, [pictured here](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/PWE_NY_Apr12.pdf) for NYC, April 12:

2. **Several temperature variables.** Our data from the NOAA consists of daily temperature maximums, minimums, and averages, for [27 United States communities](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/USMAP.png). We account for the **variable incubation rate/testing delay** by offsetting temperature and rate data, from 0 days to 2 weeks.

3. **Multiple testing and Regression.** We have two rate measurements, three temperature variables, 15 possible temperature offsets, and three correlation test statistics (Pearson, Kendall, and Spearman). We pool together pairs (rate, temperature) from each of the 27 cities, and **perform one correlation test for each possible combination of rate, temperature, offset, and correlation statistic.** For each choice of rate variable, temperature variable, and correlation statistic, we regard the 15 tests (for 15 offsets) as a multiple test; we evaluate each multiple test based on the standard rule p < 0.05 and on the Benjamini-Hochberg procedure for controlling false discovery rate (q = 0.05). Here is a series of graphs summarizing the test results for:
   * [_N'(t)_ and Spearman's test](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/result_LD_spearman.png)
   * [_N'(t)_ and Pearson's test](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/result_LD_pearson.png)
   * [_N'(t)_ and Kendall's test](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/result_LD_kendall.png)
   * [_MF(t)_ and Spearman's test](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/result_MF_spearman.png)
   * [_MF(t)_ and Pearson's test](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/result_MF_pearson.png)
   * [_MF(t)_ and Kendall's test](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/result_MF_kendall.png)
   We also include [four plots](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/deriv_mintemp_11_14.pdf) of (linear            difference rate, temperature) offset from 11 days to 14 days, with least squares lines. Colored points correspond to different cities.
   
4. **Summary of results.**


Again, thanks for clicking, and hope you enjoy! Constructive criticism and questions are very welcome.




