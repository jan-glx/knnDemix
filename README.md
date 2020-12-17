# knnDemix
R package that provides a non-parametric hypothesis test for the mixture coefficient of a mixture.
Let ![X_0\sim F_0
](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+X_0%5Csim+F_0%0A) and ![X_m\sim \alpha X_0 + (1-\alpha) X_1](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+X_m%5Csim+%5Calpha+X_0+%2B+%281-%5Calpha%29+X_1) 
Then the sum of log p-values obtained by testing for ![\hat{f}_0 \geq \alpha \hat{f}_m](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Chat%7Bf%7D_0+%5Cgeq+%5Calpha+%5Chat%7Bf%7D_m) where ![\hat{f}](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Chat%7Bf%7D) are the k nearest neighbour density estimates at the points in ![X_0](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+X_0) follows a ![\chi^2](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+%5Cchi%5E2) distribution.
This package provides a implementation of a corresponding hypothesis test for $\alpha$. Works. Needs proof.

For example this might be applied to determine the knockout rate from fluorescence cytometry measurement of a biological sample given additionally only measurements from a negative control sample.
