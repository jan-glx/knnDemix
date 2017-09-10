# knnDemix
R package that provides a non-parametric hypothesis test for the mixture coefficient of a mixture.
Let $$X_0\sim F_0$$ and $$X_m\sim \alpha X_0 + (1-\alpha) X_1$$ 
Then the sum of log p-values obtained by testing for $$\hat{f}_0 \geq alpha \hat{f}_m$$ where $$\hat{f}$$ are the k nearest neighbour density estimates at the points in $$X_0$$ follows a $$\Chi^2$$ distribution.
This package provides a implementation of a corresponding hypothesis test for $$\alpha$$. Works. Needs proof.

For example this might be applied to determine the knockout rate from fluorescence cytometry measurement of a biological sample given additionally only measurements from a negative control sample.
