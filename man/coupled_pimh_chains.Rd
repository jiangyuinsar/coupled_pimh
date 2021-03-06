% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/debiasing_algorithm.R
\name{coupled_pimh_chains}
\alias{coupled_pimh_chains}
\title{Coupled PIMH chains}
\usage{
coupled_pimh_chains(single_kernel, coupled_kernel, pimh_init, N, K = 1,
  max_iterations = Inf, preallocate = 10, verbose = FALSE)
}
\arguments{
\item{single_kernel}{function taking a state (in a vector) and an iteration, and returning
a list with a key named \code{chain_state} and containing the next state.}

\item{coupled_kernel}{function taking two states (in two vectors) and an iteration,
and returning a list with keys \code{chain_state1} and \code{chain_state2}.}

\item{K}{number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{K},
default to 1).}

\item{max_iterations}{number of iterations at which the function stops if it is still running  (default to Inf).}

\item{preallocate}{expected number of iterations, used to pre-allocate memory (default to 10).}

\item{rinit}{function taking no arguments are returning an initial state for a Markov chain.}
}
\description{
Sample two PIMH chains, each following \code{single_kernel} marginally,
and \code{coupled_kernel} jointly, until min(max(tau, K), max_iterations), where tau
is the first time at which the two chains meet (i.e. take the same value exactly).
Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
are initialized from the distribution provided in \code{rinit}.

 See \code{\link{get_hmc_kernel}}
for an example of function returning the appropriate kernels.
}
