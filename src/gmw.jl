using SpecialFunctions
using LinearAlgebra
# See "Generalized Morse Wavelets", Sofia C. Olhede, Andrew T. Walden, IEEE Transactions on signal processing, 2002

"""
  wdomain(a::Real,N::Integer,M::Integer=N)

Provide a pulsation domain scaled by a factor `a` from normalized frequencies 0 to `(M-1)/N` on a domain of size `N`"
"""
function wdomain(a::Real,N::Integer,M::Integer=N)
  domain=Array{ComplexF64}(undef,M)
  wdomain!(domain,a,N)
  return domain
end



"""
  wdomain!(domain,a::Real,N::Integer,M::Integer=N)

In-place version of `wdomain`.
"""
function wdomain!(domain::AbstractVector,a::Real,N::Integer)
  M=length(domain)
  for i in 1:M
    domain[i]=(2*pi*a)*(i-1)/N
  end
  return domain
end

"The alternative *r* parameter"
function r_βγ(β,γ)
    return (2*β+1)/γ
end

"Normalization factor, care must be taken, it depends on *gamma()* which is valid only for reasonable inputs"
function a_β_γ(k::Integer,β::Real,γ::Real,r = r_βγ(β,γ))
    return sqrt(pi * γ * 
                2^r * 
                gamma(k+1) / gamma(k+r))
end

A_l2(k,a,β,γ)=sqrt(a)*a_β_γ(k,β,γ)
function A_peak(k::Integer,β::Real,γ::Real)
  (k==0) || throw(error("No peak normalization for k != 0"))
  A=exp( ( 1 - (log(β)-log(γ)) )*β/γ)
  return A
end

function core_exp!(w::AbstractVector,a::Real,u::Real,β::Real,γ::Real,A::Real)
  for i in axes(w,1)
    x=w[i]
    θ=cis(-u*x) # Phase
    x=a*x
    w[i]=exp(β*log(x)-exp(γ*log(x)))*θ*A
  end
  return w
end

function lpoly_wave(w::AbstractVector,k::Integer,β::Real,γ::Real,r=r_βγ(β,γ),w_γ=2*exp.(γ .* log.(w)))
    return laguerrepoly(k,r-1,w_γ)
end

function _calc_wavelet!(w::AbstractVector,k::Integer,a::Real,u::Real,β::Real,γ::Real,A::Real,r=r_βγ(β,γ))
  if k !=0
    throw(error("Not implemented yet"))
#  L=lpoly_wave(w,k,β,γ,r,2*w_γ)
#  core_f = core_exp(w,β,γ,w_γ) .* L
  else
    core_exp!(w,a,u,β,γ,A)
  end
end


function gmw!(domain,k::Integer,a::Real,u::Real,β::Real,γ::Real,N::Integer,normalization::Symbol)
  (a > 0) || throw(error("Bad range need a > 0, got a=$a"))
  (β > 0) || throw(error("Bad range need β > 0, got β=$β"))
  (γ > 0) || throw(error("Bad range need γ > 0, got γ=$γ"))
  M = div(N,2)+1
  (length(domain) == M) || throw(error("Domain of size $(length(domain)), expected $M."))
  if normalization == :L2
    A=A_l2(k,a,β,γ)*sqrt(2)
  elseif normalization == :peak
    A=A_peak(k,β,γ)
  else
    throw(error("Invalid Normalization $normalization"))
  end
  wdomain!(domain,1,N)
  _calc_wavelet!(view(domain,2:M),k,a,u,β,γ,A)
  return domain
end

"""
  gmw(k::Integer,a::Real,u::Real,β::Real,γ::Real,N::Integer,normalization::Symbol)

Generalized Morse Wavelets parametrized with order `k`, scaling `a`, translation `u`, and first and second shape parameters `β` and`γ`. Valid only for `k>=0`, `a,β,γ > 0`.
"""
function gmw(k::Integer,a::Real,u::Real,β::Real,γ::Real,N::Integer,normalization::Symbol)
  M = div(N,2)+1
  domain=Array{ComplexF64}(undef,M)
  gmw!(domain,k,a,u,β,γ,N,normalization)
end

w_psi(β,γ) = exp(log(β/γ)/γ)
a_for_wpsi(wpsi,β,γ) = exp((log(β)-log(γ))/γ - log(wpsi))
peak_w(a,β,γ)= w_psi(β,γ)/a
logpeak_w(a,β,γ)= (log(β)-log(γ))/γ - log(a)
peak_n(a,β,γ,n) = peak_w(a,β,γ)*n/(2*pi)
duration(β,γ) = sqrt(β*γ)


"""
  gmw_grid(β, γ, J, Q; wmin=0,wmax=pi,phase=0.)

Creates a grid of generalized morse wavelets. The grid starts by initiliazing the first wavelet at the lowest scalce with frequency peak `wmax`, it then generates subsequent wavelets by scaling it by a factor `2^(1/JQ)`. `J`is the number of octaves and `Q` the number of inter-octaves. It results in a grid of `JQ` wavelets.
"""
function gmw_grid(β,γ,J,Q,wmin=0,wmax = pi,phase=0.)
  logpsi = log(w_psi(β,γ))
  nscale = J*Q
  logw = log(wmax) .- log(2)*(0:(nscale-1))/Q
  if !iszero(wmin)
    mask = logw .> log(wmin)
    logw=logw[mask]
  end
  a = exp.(logpsi .- logw )
  nscale = length(a)
  return [ [a[i],phase,β,γ] for i in 1:nscale ]
end



"""
  mappeakgrid(grid,n)

map a grid of generalized morse wavelets parameters to their peak in indices.
"""
mappeakgrid(grid,n)=map(x->peak_n(x[1],x[3],x[4],n),grid)

"""
  mappeakgrid(grid)

map a grid of generalized morse wavelets parameters to their peak in frequency pulsation
"""
mappeakgrid(grid)=map(x->peak_w(x[1],x[3],x[4]),grid)

#### CHECKING GMW wavelets with generalized gamma function
# using SpecialFunctions
# use SpecialFunctions.gamma_inc(α,x) -> (Γ(α)^-1)∫_0^x t^(α-1)exp(-t)dt

