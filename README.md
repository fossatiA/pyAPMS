# pyAPMS
Python implementation of common AP-MS scoring algorithms


In this repo I will upload some of the commonly used AP-MS () scoring algorithms using python3


## WD score (from CompPASS)

Developed by the [Harper lab](https://harper.hms.harvard.edu) this relies on considering both frequency of bait-prey interaction across different experiments and also the prey intensity.
The implementation in wd_comppass.py uses as input a 2D adjacency matrix in the form of

|                | Bait1          | Bait2          |
| :------------- | :------------- | :------------- |
| Prey1          | Prey1Bait1     | Prey1Bait2     |
| Prey2          | Prey1Bait1     | Prey2Bait2     |
| Prey3          | Prey3Bait1     | Prey3Bait2     |
| Prey4          | Prey4Bait1     | Prey4Bait2     |


This table needs to be a 2D numpy array. Following calculation of WD score for the provided matrix, a simulated dataset is generated to be able to assess false discovery rate.
comppass.py is ported from [this implementation](https://github.com/dnusinow/cRomppass/blob/master/R/comppass.R)
Two parameters are used to this end:\n

*iteration* number of generated datasets
*q* quantile to filter the original datasets

To run tests and see input format:

```
import numpy as np
import wd_comppass as wd

distr = numpy.random.rand(100,100)
wd_scores = wd.calc_wd_matrix(distr, iteration=1000, q=0.8, plot=True)

# filter original matrix for the wd scores which are significant
distr[wd_scores>0]=0

```

## MiST

Developed by the [Krogan lab](https://kroganlab.ucsf.edu) this metric takes into account variability, weighted sum across experiments and uniqueness of a particular bait-prey interaction. mist.py is ported from [this implementation](https://github.com/kroganlab/mist)

The implementation in mist.py uses as input a 2D adjacency matrix in the form

|                | Prey1          | Prey2          |
| :------------- | :------------- | :------------- |
| Bait1          | Prey1Bait1     | Prey2Bait2     |
| Bait1          | Prey1Bait1     | Prey2Bait2     |
| Bait2          | Prey2Bait2     | Prey2Bait2     |
| Bait2          | Prey2Bait2     | Prey2Bait2     |


The input format is a Panda dataframe.
In this example two replicate pulldown with two baits (bait1 and bait2) are performed. A complete example is available in test/mist_test.txt

Two mode can be used for calculation of mist scores

*mode* combination of abundance, reproducibility and specificity for calculating mist score. can be either *_fixed_* or *_PCA_*.
In fixed the product of abundance, reproducibility and specificity is calculated. In PCA mode a PCA is performed and the first component is used as mist score.
Finally, K mean with 2 components is used to separate false and true interaction and another column indicating the cluster is returned

To run tests and see input format:

```
import numpy as np
import mist as mist

distr = mist.create_test_data()
mist = mist.calc_mist_matrix(distr, weights='PCA', plot=True)

```

alternatively test can be run from a terminal

```
python3 mist.py
```
