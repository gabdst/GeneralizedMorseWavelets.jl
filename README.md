# Julia Library for Generalized Morse Wavelets

A Julia package for generating Generalized Morse Wavelets (GMW) (see
e.g. (Lilly and Olhede 2012)) and their derivatives.

## Theory

The first order Generalized Morse Wavelets are defined in frequency by:

$$
\widehat{\psi}_{a,u,\beta,\gamma}(w) = C_{a,\beta,\gamma}\, (aw)^\beta \exp \left( -(aw)^\gamma -iwu \right)
$$ where $C_{a,\beta,\gamma}$ is a normalization constant, usually a L2
normalization or a frequency peak normalization, $a$ and $u$ are scaling
and translation parameters and $\beta$ and $\gamma$ are shape
parameters.

In practice, Generalized Morse Wavelest can be used to compute
continuous wavelet transforms. Given a signal $x$, its wavelet
coefficients $\tilde{x}_{a,u}$ are computed at each scale $a$ and time
$u$ with fixed shape parameters $\beta$ and $\gamma$ (generally it is
recommended to choose $\beta$=1 and $\gamma$=3 as it gives the wavelet a
good localization in time-frequency space, almost a gaussian-like
localisation, see (Lilly and Olhede 2012) and references therein).

The wavelet coefficient of a signal $x$ at scale $a$ and time $u$ is
then simply computed using: $$
\tilde{x}[a,u] = \langle x,\psi_{a,u,\beta,\gamma} \rangle
$$ where $\psi_{a,u,\beta,\gamma}$ is the wavelet in time.

This package also provide Generalized Morse Wavelets at higher order,
although their use is more limited as their properties in frequency are
not yet well understood for varying shape parameters $\beta$ and
$\gamma$[^1]. For more information, see (Lilly and Olhede 2012).

## Usage

To generate the analytic part of a first order `k=0` Generalized Morse
Wavelet, normalized in energy `normalization=:L2`, with `N=1024` time
samples, with shape parameters `\beta=1`,`\gamma=3`, and at scale `a=5`
and time index `u=512`, use:

``` julia
using GMW
using FFTW
k=0
a=5
u=512
β=1
γ=3
N=1024
normalization=:L2
g_fft=gmw(k,a,u,β,γ,N,normalization) # Analytic part of the wavelet of size div(N,2)+1
g_time=irfft(g_fft,N) # real part of the corresponding complex analytic generalized morse wavelet
```

In order to compute continuous wavelet tranforms, we also provide a
handy function `gmw_grid` to initialize a bank of Generalized Morse
Wavelets whose frequency peaks are logarithmically positionned in
frequency.

With the scales $a_i$ defined by $$
a_i=a_02^{i/Q},\, \forall i=0\dots JQ-1
$$ where $J$ is the number of octaves and $Q$ the number of
inter-octaves, $a_0$ the initial lowest scale such that the highest
frequency peak is positionned at $w_\mathrm{max}$.

The parameters of the bank are then computed using:

``` julia
J=8
Q=4
wmin=0 # Minimum frequency peak allowed
wmax=pi # Maximum frequency peak allowed
normalization=:peak # Wavelet frequency peak normalized to 1
g_params=gmw_grid(β,γ,J,Q,wmin,wmax) # Get the parameters of GMW bank, returns params in the form [a,u,β,γ]
g=gmw(0,g_params[1]...,N,normalization) # Get the first wavelet of the bank (g_params starts from the lowest scale)
```

The resulting bank of wavelets have their frequency peaks $w_i$
positionned at: $$
w_i=w_\mathrm{max}2^{-i/Q},\, \forall i=0\dots JQ-1
$$ If a minimum frequency threshold $w_\mathrm{min}$ is specified, all
wavelets with frequency peaks below this threshold are removed from the
bank of filters resulting in a bank of size less than `JQ`.

Now given our bank of filters we can compute the continuous wavelet
transform of a given signal $x$:

``` julia
conv(x,g) = irfft(rfft(x) .* g,N) # For simplicity here we only compute the real part of the continuous wavelet transform
get_gmw(g_p) = gmw(0,g_p...,N,normalization)
x = randn(N)
x_cwt = [ conv(x,get_gmw(g_p)) for g_p in g_params] # Continous Wavelet transform
```

Let’s visualize the transform:

``` julia
using Plots
heatmap(hcat(x_cwt...)')
```

## TODO

- documentation, add math definitions of GMWs
- add self-dual filter bank as practical test
- first-order derivative for higher order GMWs
- second-derivatives of GMWs

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Lilly2012" class="csl-entry">

Lilly, J. M., and S. C. Olhede. 2012. “Generalized Morse Wavelets as a
Superfamily of Analytic Wavelets.” *IEEE Transactions on Signal
Processing* 60 (11): 6036–41.
<https://doi.org/10.1109/tsp.2012.2210890>.

</div>

</div>

[^1]: The reason is that it depends on the localization of the zeros of
    high order Laguerre polynomials
