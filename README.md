# Data-adaptive-PCA
Data-adaptive PCA (DPCA) is a dimension reduction technique proposed by Lim & Oh (2016). Conventional PCA uses covariance matrix of the data to extract principal components; hence, it may not perform well to skewed or asymmetric data. DPCA considers that linear combinations of asymmetric Huber density functions can represent various types of distribution such as asymmetric or bimodal one. DPCA  minimizes a weighted average of asymmetric Huber norms, instead of square norms. In practice, it can be obtained by applying ordinary PCA to pseudo-data which is a transfomation of the origin data. 

In this branch, there are two code files. One is source_code.R, and the other is example.R. In the soruce_code.R file, there is an R function comp_Pseudo which returns loadings of DPCA and contribution of each PCs. Before implementing the code, you should install several R packages, which are listed in the top line of the code. In the example.R file, there are some examples of using DPCA with daily maximum precipitation data in August during the period of year 1997 to 2008. The data are available with loading max_prec.RData file. First from grid_finder function we find the indices of East Asia region, perform and plot loading matrices of first 4 PCs. We may produce some results, which is dpca_loading2.pdf. Finally, we plot 4 PCs as a time series, which is dpca_pcs.pdf.



Reference
[1] Lim, Y., & Oh, H. S. (2016). A data-adaptive principal component analysis: Use of composite asymmetric Huber function. Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2015.1067621, 00.
