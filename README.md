# CentralSpinChain

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lovemy569.github.io/CentralSpinChain.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lovemy569.github.io/CentralSpinChain.jl/dev/)
[![Build Status](https://github.com/lovemy569/CentralSpinChain.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lovemy569/CentralSpinChain.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/lovemy569/CentralSpinChain.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/lovemy569/CentralSpinChain.jl)

The model involves a central spin of interest and a controllable spin chain environment.

$$\begin{aligned}
\hat{H}& =\hat{H}_S+\hat{H}_B+\hat{H}_I \\
\hat{H}_{S}& =\Omega\hat{\tau}^z  \\
\hat{H}_{B}& =\sum_{i=1}^L(-F i+\alpha i^2/L^2)\hat{\sigma}_i^z+\frac{1}2\sum_{i=1}^Lh_i\hat{\sigma}_i^z+\frac14\sum_{i=1}^LJ_z(i,t,\omega_b)\hat{\sigma}_i^z\hat{\sigma}_{i+1}^z+J_x(i,t,\omega_b)\hat{\sigma}_i^x\hat{\sigma}_{i+1}^x+J_y(i,t,\omega_b)\hat{\sigma}_i^y\hat{\sigma}_{i+1}^y  \\
\hat{H}_I& =\frac14\frac\gamma L\sum_{i=1}^LJ_z(t,\omega_s)\hat{\sigma}_i^z\hat{\tau}^z+J_x(t,\omega_s)\hat{\sigma}_i^x\hat{\tau}^x+J_y(t,\omega_s)\hat{\sigma}_i^y\hat{\tau}^y.
\end{aligned}$$

$$Z=\sum_{\sigma_1\sigma_2\sigma_3\ldots}e^{-E(\sigma_1,\sigma_2,\sigma_3,\ldots)/T}$$

Change the relevant parameters in the source code to perform the corresponding calculations.
