# pyAPMS
Python implementation of common AP-MS scoring algorithms


In this repo I will upload some of the commonly used AP-MS () scoring algorithms using python3


## WD score (from CompPASS)

Developed by the [Harper lab](https://harper.hms.harvard.edu) this relies on considering both frequency of bait-prey interaction across different experiments and also the prey intensity.
The implementation in mist.py uses as input a 2D adjacency matrix in the form of

|                | Bait1          | Bait2          |
| :------------- | :------------- | :------------- |
| Prey1          | Prey1Bait1     | Prey1Bait2     |
| Prey2          | Prey1Bait1     | Prey2Bait2     |
| Prey3          | Prey3Bait1     | Prey3Bait2     |
| Prey4          | Prey4Bait2     | Prey4Bait2     |


This table needs to be a 2D numpy array. Following calculation of WD score for the provided matrix, a simulated dataset is generated to be able to assess false discovery rate.
Two parameters are used to this end:\n

*iteration* number of generated datasets
*q* quantile to filter the original datasets

To run tests and see input format:

```
import numpy as np
import mist as mist

distr = numpy.random.rand(100,100)
wd_scores = mist.calc_wd_matrix(distr, iteration=1000, q=0.8, plot=True)

# filter original matrix for the wd scores which are significant
distr[wd_scores>0]=0

```
