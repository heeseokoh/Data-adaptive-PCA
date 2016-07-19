# Data-adaptive-PCA
Data-adaptive PCA is a dimension reduction technique proposed by Lim & Oh (2016). Conventional PCA uses covariance matrix of the data to extract principal components, so it has limitation that it may not perform well to skewed or asymmetric data. They discovered that linear combinations of asymmetric Huber density functions can represent various types of distribution, such as asymmetric or bimodal ones. Motivated from this, in data-adaptive PCA, they minimizes weighted average of asymmetric Huber norms, instead of square norms. In practice, we can obtain data-adaptive principal components as applying ordinary PCA to pseudo-data which is transfomation of origin data. 

In this branch, there are two code files. 




Reference
[1] Lim, Y., & Oh, H. S. (2015). A data-adaptive principal component analysis: Use of composite asymmetric Huber function. Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2015.1067621, 00.
