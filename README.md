# COVID-19 Spread in the U.S. and Temperature: a Statistical Analysis
Thanks for clicking on my project!


**Please see the detailed [summary](https://github.com/seaneli/CHRP-Competition-Repository/blob/master/PRESENTATION.pdf) of my project, typed and explained neatly, with LaTeX.**


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
   
4. **Summary of results.** In all cases, the test correlations are mild, with absolute value < 0.15. 

   * **Linear Diﬀerence Rate:** Kendall and Spearman tests suggest _N'(t)_ is possibly **positively correlated with the minimum temperature from 1-2 weeks before**, while the (less appropriate) Pearson test suggests _N'(t)_ is possibly **negatively correlated with the maximum temperature, from 0 days-2 weeks before.** All tests are evaluated with the Benjamini-Hochberg procedure, with FDR control at 0.05. Least-squares lines seem to support the above negative Pearson correlations. 
   * **Multiplicative Factor Rate:** All correlations are very mild, with absolute value less than 0.08. Kendall and Spearman tests suggest **no correlation between this rate and any temperature variables.** Pearson’s test suggests **_MF(t)_ is negatively correlated with minimum temperature oﬀset 14 days, and with average temperature, oﬀset 14 days**. The Pearson correlation result is supported by the slope of least-squares lines.
   
 5. **Possible explanation/ future research** Our results suggest the growth rate of COVID-19 cases may be positively correlated with minimum temp. from 1-2 weeks before, and may be negatively correlated with maximum temperature overall. This suggests the highest growth rates appear when the temperature is not too low and not too high (approx. between 50 and 80 degrees). There are many possible explanations for the 2-week delay for minimum temperature effects. **One possible explanation is that the extreme minimum temperatures prolong the COVID-19 incubation period.** To test this hypothesis, we need reliable data on the COVID-19 incubation period (i.e. individual patient data. Date of first exposure to infected person, date of first onset of symptoms. This is not easy data to find, since the date of first exposure will often be impossible to determine.)


Again, thanks for clicking! Constructive criticism and questions are very welcome. 

email: seaneli@rice.edu

EDIT 2021: minor wording edits in presentation.pdf


