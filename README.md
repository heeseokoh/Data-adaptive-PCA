# Data-adaptive-PCA
Data-adaptive PCA is a dimension reduction technique proposed by Lim & Oh (2016). Conventional PCA uses covariance matrix of the data to extract principal components; hence, it may not perform well to skewed or asymmetric data. They discovered that linear combinations of asymmetric Huber density functions can represent various types of distribution, such as asymmetric or bimodal ones. Motivated from this, in data-adaptive PCA, they minimizes weighted average of asymmetric Huber norms, instead of square norms. In practice, we can obtain data-adaptive principal components as applying ordinary PCA to pseudo-data which is transfomation of origin data. 

In this branch, there are two code files. One is source_code.R, and the other is example.R. In soruce_code.R file, there is R function comp_Pseudo which returns loadings of data-adaptive PCA and contribution of each PCs. Before implement the code you should install several R packages, and it is listed on the top of code file.  In example.R file, there are some examples of using data-adaptive PCA using the daily maximum precipitation data in August from CMAP during the period of year 1997 to 2008. The data is available if you load max_prec.RData file. First we only see East Asia region, and grid_finder function finds the index of such region. Next we perform and plot loading matrices of first 4 PCs; result is dpca_loading2 file. Finally we plot 4 PCs as a time series; result is dpca_pcs file.




Reference
[1] Lim, Y., & Oh, H. S. (2016). A data-adaptive principal component analysis: Use of composite asymmetric Huber function. Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2015.1067621, 00.
