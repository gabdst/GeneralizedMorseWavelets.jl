# First Derivatives for k > 0 is WIP
function logA_peak_βγ_deriv(β::Real,γ::Real)
  ∂β=(-1/γ)*(log(β)-log(γ))
  ∂γ=-(β/γ)*∂β
  return (∂β,∂γ)
end

"Jacobian of the transformation (β,γ) -> (r,γ)"
function _∇zg_A(β::Real,γ::Real,r=r_βγ(β,γ)) # Where z=(β,γ)
    return [ 2/γ -r/γ ;  0  1 ]
end

"Jacobian of the normalization factor"
function _∇gA(k,β,γ,r=r_βγ(β,γ),A = a_β_γ(k,β,γ))
    return (A/2).*[ (log(2)-digamma(k+r)) 1/γ ] 
end

"Jacobian of (β,γ) -> (r-1,2w^γ)"
function _∇zg_L(w::AbstractArray,β::Real,γ::Real,r=r_βγ(β,γ),w_γ=2*exp.(γ .* log.(w)))
  return [ 2/γ -r/γ ; zeros(length(w)) log.(w).*w_γ]
end

function a_β_γ_deriv(k::Integer,β::Real,γ::Real,r=r_βγ(β,γ),A = a_β_γ(k,β,γ))
  k ==0 || throw(error("Not Implemented"))
  ∇gA = _∇gA(k,β,γ,r,A)
  ∇zg = _∇zg_A(β,γ,r)
  out = ∇gA*∇zg
  ∂β = out[1,1]
  ∂γ = out[1,2]
  return (∂β,∂γ)
end

function core_exp_βγderiv(w::AbstractArray,β::Real,γ::Real,w_γ::AbstractArray=exp.(γ .* log.(w)))
  glog = core_exp(w,β,γ,w_γ) .* log.(w)
  ∂β = glog
  ∂γ = -w_γ.*glog
  return (∂β,∂γ)
end

function lpoly_wave_βγderiv(w::AbstractArray,k::Integer,β::Real,γ::Real,r=r_βγ(β,γ),w_γ=2*exp.(γ .* log.(w))) # here w_γ is doubled, OK
  ∂α_l = α_derivative_laguerre(k,r-1,w_γ)
  ∂x_l = x_derivative_laguerre(k,r-1,w_γ)
  ∇zg = _∇zg_L(w,β,γ,r,w_γ)
  ∂β = ∂α_l*∇zg[1,1] #+ ∂x_l*J[2:end,1] J[2:end,1] is zero
  ∂γ = ∂α_l*∇zg[1,2] + ∂x_l .* ∇zg[2:end,2] # Could be done with a tensor product instead by using a larger jacobian
  #∂β=(2/γ)*∂α_l
  #∂γ= γ*log.(w) .* w_γ .* x_derivative_laguerre(k,r-1,w_γ) - r*∂α_l # the old way was wrong I think, missing an inverse gamma somewhere
  return (∂β,∂γ)
end

function _calc_wavelet_βγderiv(w::AbstractArray,k::Integer,β::Real,γ::Real,r=r_βγ(β,γ),  w_γ = exp.(γ .* log.(w)))
    A = a_β_γ(k,β,γ)
    ∂A = a_β_γ_deriv(k,β,γ,r,A)
    core= core_exp(w,β,γ,w_γ)
    ∂core  = core_exp_βγderiv(w,β,γ,w_γ)
    if k > 0 
        L=lpoly_wave(w,k,β,γ,r,2*w_γ)
        ∂L=lpoly_wave_βγderiv(w,k,β,γ,r,2*w_γ)
        Ccore=L.*A
        CL=A.*core
        CA=core.*L
        ∂β=sqrt(2)*(Ccore.*∂core[1] + CL .* ∂L[1] + CA .* ∂A[1])
        ∂γ=sqrt(2)*(Ccore.*∂core[2] + CL .* ∂L[2] + CA .* ∂A[2])
    else
        Ccore=A
        CA=core
        ∂β=sqrt(2)*(Ccore.*∂core[1] +  CA .* ∂A[1])
        ∂γ=sqrt(2)*(Ccore.*∂core[2] +  CA .* ∂A[2])
    end
    return (∂β,∂γ)
end

function u_derivative(g,k,a,u,β,γ,N,normalization)
  M=length(g)
  return -1im*wdomain(1,N,M).*g
end

function w_derivative(g,k,a,u,β,γ,N,normalization)
  M=length(g)
  r = r_βγ(β,γ)
  ∂w = wdomain(a,N,M)
  if k != 0
    throw(error("Not Implemeted Yet"))
    Agprev = (k+r-1)*γ*a_β_γ(k,β,γ)/a_β_γ(k-1,β,γ)
    gprev = Agprev*gmw(k-1,a,u,β,γ,N,normalization)
    for i in 2:M
      Ag = ( β  - γ*exp(γ*log(∂w[i]))  + γ*k)
      ∂w[i] =  (Ag*g[i] - gprev[i])/∂w[i]
    end
  else
    for i in 2:M
      Ag = ( β  - γ*exp(γ*log(∂w[i])) )
      ∂w[i] = Ag*g[i]/∂w[i]
    end
    return ∂w
  end
end

function a_derivative(g,k,a,u,β,γ,N,normalization) # Phase is conserved
  M=length(g)
  if k == 0 # Faster than the following
    ∂a = wdomain(a,N,M)
    if normalization==:L2
      for i in 2:M
        ∂a[i] = (1/a)*((2*β+1)/2 - γ*exp(γ*log(∂a[i]))) * g[i]
      end
      return ∂a
    end
    if normalization==:peak
      for i in 2:M
        ∂a[i] = (1/a)*(β - γ*exp(γ*log(∂a[i]))) * g[i]
      end
      return ∂a
    end
  else
    throw(error("Not implemented yet"))
    w = wdomain(1,N,M)
    return w_derivative(g,k,a,u,β,γ,N,normalization).*w .+ (0.5/a)*g
  end
end

function γ_derivative(g,k,a,u,β,γ,N,normalization)
  M=length(g)
  if k == 0 # FASTER
    ∂γ=wdomain(a,N,M)
    if normalization == :L2
      r = r_βγ(β,γ)
      ddig = log(2) - digamma(r)
      for i in 2:M
        lw=log(∂γ[i])
        ∂γ[i] = ( (1/(2*γ))*(1-r*ddig) - lw*exp(γ*lw)) * g[i]
      end
      return ∂γ
    end
    if normalization == :peak
      ∂logA∂β,∂logA∂γ=logA_peak_βγ_deriv(β,γ)
      for i in 2:M
        lw=log(∂γ[i])
        ∂γ[i] = ( ∂logA∂γ - lw*exp(γ*lw)) * g[i]
      end
      return ∂γ
    end
  else
    throw(error("Not implemented yet"))
    w = wdomain(a,N)[1:M]
    Cphase = tau_expo(u,N)*sqrt(a) # Keep the phase
    partials=_calc_wavelet_βγderiv(w[2:M],k,β,γ)
    ∂γ = vcat(0,Cphase[2:M] .* partials[2])
    return ∂γ
  end
end

function β_derivative(g,k,a,u,β,γ,N,normalization)
  M=length(g)
  if k == 0 # FASTER
    ∂β=wdomain(a,N,M)
    if normalization == :L2
      r = r_βγ(β,γ)
      ddig = log(2) - digamma(r)
      for i in 2:M
        lw=log(∂β[i])
        ∂β[i] = ( (1/γ)*ddig + lw)*g[i]
      end
      return ∂β
    end
    if normalization == :peak
      ∂logA∂β,∂logA∂γ=logA_peak_βγ_deriv(β,γ)
      for i in 2:M
        lw=log(∂β[i])
        ∂β[i] = ( ∂logA∂β + lw)*g[i]
      end
      return ∂β
    end
  else
    throw(error("Not implemented yet"))
    w = wdomain(a,N)[1:M]
    Cphase = tau_expo(u,N)*sqrt(a)# Keep phase
    partials=_calc_wavelet_βγderiv(w[2:M],k,β,γ) 
    ∂β = vcat(0,Cphase[2:M] .* partials[1])
    return ∂β
  end
end
