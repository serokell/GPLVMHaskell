# GPLVMHaskell

## Abstract

Gaussian process latent variable model is a method of machine learning that is widely used in areas such as face recognition, pose estimation, etc. We propose a bit different approach to the representation of these classes of models via functional programming in Haskell and its  advantages.

## Background

### (Probabilistic) principal component analysis (PPCA)

Principal component analysis is a technique for dimensionality reduction. In other words, PCA allows one to reduce the size of data without a loss of general information.

We perform this technique by extracting data and computing eigenvectors and eigenvalues of a covariance matrix to form a feature vector. Eigenvectors are the new axis and eigenvalues show the values of those axes. The larger the variance of data with respect to the axis, the larger the corresponding eigenvalue. Afterwards, we can get rid of axes with low eigenvalues without large losses.

In contrast to PCA, [Probabilistic PCA](http://www.robots.ox.ac.uk/~cvrg/hilary2006/ppca.pdf) is a latent variable model. We have only the observed data set without its preimage. We assume that the given observed data is obtained from some linear function. Every observed vector is represented as a sum of **W** matrix multiplied by the corresponding preimage vector and noise with a zero-mean Gaussian distribution.

### Gaussian process

A [Gaussian process](https://arxiv.org/abs/1505.02965) is a random process with a normal distribution, a mapping from a data space to some space of observed values. We define a Gaussian process by its expected value and a covariance matrix computed by some kernel function.

In other words, a function is said to be a Gaussian process (GP) if and only if for all input vectors, their joint distribution is normal with some mean and a covariance matrix. As we have told before, we obtain a covariance matrix using a kernel function. A covariance matrix (symmetric and positive) in turn allows us to define that a mapping between the data set and the observable values of output preserves "proximity" between the points of this data set.
Note that we shouldn't find the parameters of those functions because only the hyperparameters of a kernel function matter. Thus, a Gaussian process provides a non-parametric way to solve a regression problem.

That is, a GP is a prior over functions and prior depends on a mean and a kernel function. We infer a distribution over possible functions to predict new inputs. It is clear that we don't need all the functions, that would be too much. So, we put constraints on desirable functions, and the distribution is defined on some subset of functions most of the time. After that, we use this distribution to get a posterior over functions that belong to this subset.

## GPLVM summary

Gaussian process latent variable model is a version of probabilistic principal component analysis [initially proposed](https://papers.nips.cc/paper/2540-gaussian-process-latent-variable-models-for-visualisation-of-high-dimensional-data.pdf) by N. D. Lawrence. GPLVM was made for visualization of higher dimensional data via dimensional reduction. In a matter of fact, GPLVM is a kind of combination of PPCA and a Gaussian process.

A Gaussian process is a non-parametric approach that allows one to find a distribution over some class of functions with respect to observed data.
In GPLVM, we assume that we have generated the observed data from a low dimensional data, and this generation process is described by some function `f` up to Gaussian distributed noise `e`. `f` has a GP prior, i.e. GP prior with a zero mean and some kernel function k with some hyperparameters. Note that these functions are not necessarily linear.

Currently GPLVM is still not implemented.

## Haskell implementation

### Why Haskell?

We may describe the general case of GPLVM representing the basic notions of this model as datatypes. This approach allows us to consider Haskell of GPLVM closer to the theory.
In contrast to Python, we may use Haskell type system to propose a more universal consideration of GPLVM independently of the subject matter. That is, we propose a general GPLVM as a set of abstract datatypes and related polymorphic functions.
Thus, our description of this model consists of a set of primitive idioms that allow one to implement a much more transparent code that might be easily structured into the system of modules.

Examples of proposed definitions (TODO: rewrite examples via actual code updates):

```haskell
newtype ObservedData a = ObservedData { _unObservedData :: Matrix a }
```

```haskell
data GaussianProcess a = GaussianProcess {
      _pointsGP :: Matrix a
    , _outputGP :: ObservedData a
    , _kernelGP :: Matrix a -> Matrix a -> Matrix a
    , _distrGP ::  Distribution a
    }
```

```haskell
data GaussianProcessLatentVariable a = GaussianProcessLatentVariable {
      _GPLVMObserved :: ObservedData a
    , _kerGPLVM :: Matrix a -> Matrix a -> Matrix a
    , _distrGPLVM ::  Distribution a
    }
```

One could test PPCA algorithm via the application: 

```$ GPLVMHaskell-exe test.pca.txt 1000 False```  

The first parameter is a file with rows as observation points. 
The second parameter is number of iterations. 
The third one should be False for non-typesafe version and vice versa. 

If there are "NaNs" in a file, then the algorithm version which could handle 
lost or untrusted data would be applied. In this case the allication will return not only ***W*** matrix and log-likelihood but also restored original matrix.  