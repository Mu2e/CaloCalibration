When a cosmic ray hits one of the calorimeter crystals, it produces two signals (one for each SiPM), that are sampled every 5 ns.
Caloreco.C fits each sigal with a template function over a certain time interval (start and end are relative to the peack time). The timing of the event is then defined as when the fit reaches 20 % of the maximum. Now AnaDriver2_1 can compute the timing difference between the two signals from each crystal over amny events, we take the simga of this timing difference and compute its mean over the crystals, then we plot it against the interval.

Itnervals to be tested are in t_intervals.py
To run the analysis use time_tests.ipynb (preferably on the [EAF](https://mu2ewiki.fnal.gov/wiki/Elastic_Analysis_Facility_(EAF)) for imporved parallelism).
Using v_fit_tries.py sould also work fine.