
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **DPCA_source_code** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : DPCA_source_code

Published in : 'Lim, Y., Oh, H. S. (2016): A data-adaptive principal component analysis: Use of
composite asymmetric Huber function. Journal of Computational and Graphical Statistics'

Description : Function comp_Pseudo which returns loadings of DPCA and contribution of each PCs

Keywords : principal component, pca, asymmetric, bimodal, Huber norm, weighted

See also : DPCA_example

Author : OH H. S., Lim Y.

Submitted : 20160728 Burdejova Petra

```


### R Code:
```r
library(rrcov)
library(evd)
library(sn)
library(mvtnorm)


##### Density function ######
#############################

dens = function(x, y) {
    index = density(y)$x
    fit   = density(y)$y[which.min(abs(index - x))]
    return(fit)
}


##### Modify rho_tau,c function #####
#####################################

psi = function(xx, tau = 0.5) {
    cut    = 1.345
    result = matrix(nrow = nrow(xx), ncol = ncol(xx))
    for (i in 1:nrow(xx)) {
        for (j in 1:ncol(xx)) {
            x = xx[i, j]
            if (x < (-cut)) {
                u = (tau - 1)/2
            }
            if (x >= (-cut) & x < 0) {
                u = (1 - tau) * x
            }
            if (x >= 0 & x < cut) {
                u = tau * x
            }
            if (x >= (cut)) {
                u = tau/2
            }
            result[i, j] = u
        }
    }
    return(result)
    
} 	


		
#######################################################################
### Proposed loss function  :weighted average of rho_tau,c function
#######################################################################

com_psi = function(x, opt_w, b, tau_range) {
    A = matrix(0, nrow = nrow(x), ncol = ncol(x))
    for (j in 1:nrow(x)) {
        for (i in 1:length(tau_range)) {
            z = opt_w[i, j] * psi(as.matrix(x[j, ] - b[i, j]), tau_range[i])
            A[j, ] = A[j, ] + z
        }
    }
    return(A)
} 


##########################################################################
########################## Main function #################################
##### Inputs: 
#####   x: data matrix 
#####   q: number of the  component to extract  
##### Outputs: 
#####   result[[1]]: Loading matrix 
#####   result[[2]]: Eigen values 
##########################################################################

comp_Pseudo = function(x, q) {
    n = nrow(x)
    p = ncol(x)
    
    psi <- function(xx, tau = 0.5) {
        cut = 1.345
        result = matrix(nrow = nrow(xx), ncol = ncol(xx))
        for (i in 1:nrow(xx)) {
            for (j in 1:ncol(xx)) {
                x <- xx[i, j]
                if (x < (-cut)) {
                  u = (tau - 1)/2
                }
                if (x >= (-cut) & x < 0) {
                  u = (1 - tau) * x
                }
                if (x >= 0 & x < cut) {
                  u = tau * x
                }
                if (x >= (cut)) {
                  u = tau/2
                }
                result[i, j] <- u
            }
        }
        return(result)      
    }
    
    
    com_psi = function(x, opt_w, b, tau_range) {
        A = matrix(0, nrow = nrow(x), ncol = ncol(x))
        for (j in 1:nrow(x)) {
            for (i in 1:length(tau_range)) {
                z = opt_w[i, j] * psi(as.matrix(x[j, ] - b[i, j]), tau_range[i])
                A[j, ] = A[j, ] + z
            }
        }
        return(A)
    }
    
    tau_range = seq(0.05, 0.95, by = 0.05)
    opt_w = matrix(nrow = length(tau_range), ncol = nrow(x))
    for (i in 1:nrow(x)) {
        a = vector()
        for (k in 1:length(tau_range)) {
            a[k] <- dens(quantile(x[i, ], tau_range[k]), x[i, ])
        }
        
        opt_w[, i] = (a)
        if (max(a) > 1) {
            opt_w[, i] = a/sum(a)
        }
    }
    
    robpca = PcaHubert(x, q)
    B = getLoadings(robpca)[, 1:q]
    
    hatx  = x %*% B %*% t(B)
    b_tau = matrix(nrow = length(tau_range), ncol = (nrow(x)), 0)
    Z     = hatx + com_psi(((x) - (hatx)), opt_w, b_tau, tau_range)
    # Z= hatx +
    # (t(matrix(rep(sd(x-hatx),each=n),p,n)))*com_psi(((x)-(hatx))/(t(matrix(rep(sd(x-hatx),each=n),p,n))),
    # opt_w , b_tau, tau_range)
    if (n > p) {
        BX   = x %*% B
        B    = lm(Z ~ BX - 1)$coeff
        hatx = x %*% t(B) %*% B
        Z    = hatx + com_psi(((x) - (hatx)), opt_w, b_tau, tau_range)
    }
    if (n < p) {
        BX   = x %*% B
        B    = lm(Z ~ BX - 1)$coeff
        hatx = x %*% t(B) %*% B
        Z    = hatx + com_psi(((x) - (hatx)), opt_w, b_tau, tau_range)     
    }
    
    dd  = vector()
    ans = list()
    m   = 0
    while (m < 25) {
        
        m        = m + 1
        ans[[m]] = list()
        
        if (n > p) {
            BX = x %*% t(B)
            B  = lm(Z ~ BX - 1)$coeff
            hatx1 = x %*% t(B) %*% B
            Z = hatx1 + com_psi(((x) - (hatx1)), opt_w, b_tau, tau_range)
             
            BX = x %*% t(B)
            B  = lm(Z ~ BX - 1)$coeff
            hatx2 = x %*% t(B) %*% B
            Z = hatx2 + com_psi(((x) - (hatx2)), opt_w, b_tau, tau_range)              
        }
        
        if (n < p) {
            BX = x %*% t(B)
            B  = lm(Z ~ BX - 1)$coeff
            hatx1 = x %*% t(B) %*% B
            Z = hatx1 + com_psi(((x) - (hatx1)), opt_w, b_tau, tau_range)
             
            BX = x %*% t(B)
            B  = lm(Z ~ BX - 1)$coeff
            hatx2 = x %*% t(B) %*% B
            Z = hatx2 + com_psi(((x) - (hatx2)), opt_w, b_tau, tau_range)          
        }
        
        ans[[m]][[1]] = B
        ans[[m]][[2]] = diag(cov(BX))  #svd(cov(Z))$d
        dd[m] = mean((hatx1 - hatx2)^2)
        # print(mean( ( hatx1-hatx2 ) ^2 ) ) print(diag(cov(BX))[1:q])
    }
    
    res      = list()
    res[[1]] = ans[[which.min(dd)]][[1]]
    res[[2]] = ans[[which.min(dd)]][[2]]
    # print( res[[2]])
    return(res)
}
```
