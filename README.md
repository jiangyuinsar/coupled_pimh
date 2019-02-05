# Coupled PIMH

Code to produce figures using in Unbiased Smoothing using Particle Independent Metropolis–Hastings.
Package can be installed through downloading the repository, opening coupled_pimh.Rproj in RStudio and running devtools::document().
Three models provided - linear-Gaussian state space model, stochastic volatility and Markov jump process - as well as some example code for SMC samplers. 

It is recommended to begin by running 

'inst/reproduce/stochastic_kinetic/stochastic_kinetic_plot_state.R'

to gain unbiased estimates of smoothing distributions for the Markov jump process. 
