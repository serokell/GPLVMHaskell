# GPLVMHaskell

## Background

### (Probabilistic) principal component analysis (PPCA)

Principal component analysis is a technique for dimensionality reduction.

### Gaussian process

Gaussian process is a random process with normal distribution, a mapping from a latent space to some space of observed values. Gaussian process is defined by its excepted value and covariance matrix, that we obtain using some kernel function.

In other words, a function is said to be a gaussian process if and only if for all input vectors, their joint distribution is normal with some excepted value from an arbitrary mean function.

As we told before, we obtain a covariance matrix using a kernel function. Note that, covariance matrix is a symmetric and positive.

That is, where GP is a prior over functions and prior dependes on a mean and this covariance function.

## GPLVM summary


Gaussian process latent variable model is a version of probabilistic principal component analysis [initially proposed](https://papers.nips.cc/paper/2540-gaussian-process-latent-variable-models-for-visualisation-of-high-dimensional-data.pdf) by N. D. Lawrence. GPLVM was provided for visualisation of higher dimensional data.
In a matter of fact, GPLVM is a kind of combination of PPCA and gaussian process.

Gaussian process is a non-parametric approach that allows one to find a function with respect to observed data.

In GPLVM we assume that we generated the observed data from a low dimensional data, and this generation process is described by some non-linear function `f` up to gaussian distributed noise `e`. `f` has a GP prior, i.e. GP prior with a zero mean and some kernel function $k$ with some hyperparameters.
The example of a kernel function is the RBF-kernel.

TODO: fill the explanation

## Applications

TODO: some words on computer vision: pose estimation, face recognition.


## Haskell implementation

### Why Haskell?

TODO: write down arguments


Examples of definitions:

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

## What we need to do for GPLVM implementation.

Final goal: formalization of GPLVM via Haskell type system to propose a general framework to GPLVM.

1. Consider the general case of PCA and PPCA and formalize in Haskell. Provide the way to work with PCA and PPCA;
2. Similarly with Gaussian process;
3. Research the current abilities of probabilistic programming, machine learning, and Bayes inference in Haskell;
4. Similarly with linear algebra in Haskell, the most of observed libraries are quite old;
5. Define the list of required mathematical constructions and operations from algebra, probability theory and numeric analysis.
6. Formalize GPLVM using the background above:
6.1. General GPLVM
6.2. Demo application.
