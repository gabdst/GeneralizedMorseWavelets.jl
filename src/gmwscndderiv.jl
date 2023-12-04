### Everything below is WIP

function _∇∇zg_A(β::Real,γ::Real,r=r_βγ(β,γ))
    return [ 0 -2/γ^2 ;
            0 0;
            -2/γ^2 2*r/γ^2 ;
            0 0 ]
end

# ∇zg jacobian for the transformation z=(β,γ) -> (r,γ)=g(β,γ)=g(z)
# Using dvec(∇za) = (I kron ∇ga)dvec(∇zg) + (∇zg' kron I)dvec(∇ga)
function a_β_γ_2deriv(k::Integer,β::Real,γ::Real,r=r_βγ(β,γ),A=a_β_γ(k,β,γ))
    ∇gA = _∇gA(k,β,γ,r,A)
    ∇∇gA = _∇∇gA(k,β,γ,r,A,∇gA)
    ∇zg = _∇zg_A(β,γ,r)
    ∇∇zg = _∇∇zg_A(β,γ,r)
    ∇∇zA=∇∇gA*∇zg
    out = (∇zg')*∇∇zA + kron(I(2),∇gA)*∇∇zg

#    ∂β2 = (∇gA*∇∇zg[1:2,1])[1] + (∇zg')[1,1:2]'*∇∇gA[1:2,1]
#    ∂βγ = (∇gA*∇∇zg[1:2,2])[1] + (∇zg')[1,1:2]'*∇∇gA[1:2,2]
#    ∂γ2 = (∇gA*∇∇zg[3:4,2])[1] + (∇zg')[2,1:2]'*∇∇gA[1:2,2]
    ∂β2 = out[1,1]
    ∂βγ = out[1,2]
    ∂γ2 = out[2,2]
    return (∂β2,∂βγ,∂γ2)
end

function core_exp_∂βγ∂βγderiv(w::AbstractArray,β::Real,γ::Real,w_γ::AbstractArray=exp.(γ .* log.(w)))
    glog2 = core_exp(w,β,γ,w_γ) .* (log.(w).^2)
    ∂β2 = glog2
    ∂βγ = -glog2 .* w_γ
    ∂γ2 = glog2 .* w_γ.*(w_γ .- 1)
    return (∂β2,∂βγ,∂γ2)
end

function _∇∇zg_L(w::AbstractArray,β::Real,γ::Real,r=r_βγ(β,γ),w_γ=2*exp.(γ .* log.(w)))
    N = length(w_γ)
    return [ 0 -2/γ^2; 
            zeros(N) zeros(N) ; 
            -2/γ^2 2*r/γ^2 ; 
            zeros(N) ((log.(w).^2).*w_γ) ]
end # To check

# Here we want to use ∇∇zL=dvec(∇gL∇zg)/dz = (I kron ∇gL) ∇∇zg + (∇zg' kron I) ∇∇gL∇zg
function lpoly_wave_∂βγ∂βγderiv(w::AbstractArray,k::Integer,β::Real,γ::Real,r=r_βγ(β,γ),w_γ=2*exp.(γ .* log.(w))) # here w_γ is doubled, OK
    N=length(w_γ)
    ∂α_l = α_derivative_laguerre(k,r-1,w_γ)
    ∂x_l = x_derivative_laguerre(k,r-1,w_γ)
    ∇gL = [ ∂α_l ∂x_l ]
    ∂α2,∂αx,∂x2 = x2_αx_α2_derivative_laguerre(k,r-1,w_γ)
    ∇∇gL = [ ∂α2 ∂αx ; ∂αx ∂x2 ]

    ∇zg= _∇zg_L(w,β,γ,r,w_γ) # (N+1) x 2
    ∇∇zg = _∇∇zg_L(w,β,γ,r,w_γ) # (2(N+1) x 2 )

    # We do the calculation by hand since we actually are deriving for each w[i], ideally we could use of a tensor product implem.

    # (I kron ∇gL) ∇∇zg
    out = [ (∇gL[:,1]*∇∇zg[1,1]+∇gL[:,2].*∇∇zg[2:N+1,1]) (∇gL[:,1]*∇∇zg[1,2]+∇gL[:,2].*∇∇zg[2:N+1,2]) ; 
           (∇gL[:,1]*∇∇zg[N+2,1]+∇gL[:,2].*∇∇zg[N+3:end,1])  (∇gL[:,1]*∇∇zg[N+2,2]+∇gL[:,2].*∇∇zg[N+3:end,2]) ]

    # (∇zg' kron I) ∇∇gL∇zg
    JJ = [ (∇∇gL[1:N,1] * ∇zg[1,1]+∇∇gL[1:N,2].*∇zg[2:end,1]) (∇∇gL[1:N,1] * ∇zg[1,2]+∇∇gL[1:N,2].*∇zg[2:end,2]) ;
          (∇∇gL[N+1:end,1] * ∇zg[1,1]+∇∇gL[N+1:end,2].*∇zg[2:end,1]) (∇∇gL[N+1:end,1] * ∇zg[1,2]+∇∇gL[N+1:end,2].*∇zg[2:end,2])]
    # Warning: with ∇zg transposed !
    out += [ (∇zg[1,1]*JJ[1:N,1]+∇zg[2:end,1].*JJ[N+1:end,1]) (∇zg[1,1]*JJ[1:N,2]+∇zg[2:end,1].*JJ[N+1:end,2]) ;
             (∇zg[1,2]*JJ[1:N,1]+∇zg[2:end,2].*JJ[N+1:end,1]) (∇zg[1,2]*JJ[1:N,2]+∇zg[2:end,2].*JJ[N+1:end,2])]
    ∂β2=out[1:N,1]
    ∂βγ=out[N+1:end,1]
    ∂γ2=out[N+1:end,2]
    return (∂β2, ∂βγ, ∂γ2)
end

function _calc_wavelet_∂βγ∂βγderiv(w::AbstractArray,k::Integer,β::Real,γ::Real,r=r_βγ(β,γ),  w_γ = exp.(γ .* log.(w)))
    A = a_β_γ(k,β,γ)
    ∂A = a_β_γ_deriv(k,β,γ,r,A)
    ∂A2 = a_β_γ_2deriv(k,β,γ,r,A)
    G= core_exp(w,β,γ,w_γ)
    ∂G  = core_exp_βγderiv(w,β,γ,w_γ) # (n,n)
    ∂G2  =  core_exp_∂βγ∂βγderiv(w,β,γ,w_γ) # (n,n,n)
    if k > 0 # Broken somewhere, check gmwgrad.jl
        L=lpoly_wave(w,k,β,γ,r,2*w_γ)
        ∂L=lpoly_wave_βγderiv(w,k,β,γ,r,2*w_γ)
        ∂L2=lpoly_wave_∂βγ∂βγderiv(w,k,β,γ,r,2*w_γ)
        CG=L.*A
        CL=A.*G
        CA=G.*L

        # Common vectors to ∂β2 and ∂βγ
        c_βγA= ( L .* ∂G[1] + G .* ∂L[1] ) 
        c_βγL= ( A .* ∂G[1] + G .* ∂A[1] )
        c_βγG= ( A .* ∂L[1] + L .* ∂A[1] )

        ∂β2 = sqrt(2)*(∂G2[1] .* CG + ∂L2[1] .* CL + ∂A2[1] .* CA  + ∂A[1] .* c_βγA + ∂L[1] .* c_βγL +  ∂G[1] .* c_βγG)
        ∂βγ = sqrt(2)*(∂G2[2] .* CG + ∂L2[2] .* CL + ∂A2[2] .* CA  + ∂A[2] .* c_βγA + ∂L[2] .* c_βγL +  ∂G[2] .* c_βγG)

        # Vectors for ∂γ2
        c_βγA= ( L .* ∂G[2] + G .* ∂L[2] ) 
        c_βγL= ( A .* ∂G[2] + G .* ∂A[2] )
        c_βγG= ( A .* ∂L[2] + L .* ∂A[2] )
        ∂γ2 = sqrt(2)*(∂G2[3] .* CG + ∂L2[3] .* CL + ∂A2[3] .* CA + ∂A[2] .* c_βγA + ∂L[2] .* c_βγL +  ∂G[2] .* c_βγG)
        return (∂β2,∂βγ,∂γ2)
    else # OK
        CG=A
        CA=G
        ∂β=sqrt(2)*(CG.*∂G[1] +  CA .* ∂A[1])
        ∂β2=sqrt(2)*(∂A[1].*∂G[1] + CG.*∂G2[1] + ∂G[1].*∂A[1] + CA.*∂A2[1])
        ∂βγ=sqrt(2)*(∂A[2].*∂G[1] + CG.*∂G2[2] + ∂G[2].*∂A[1] + CA.*∂A2[2])
        ∂γ2=sqrt(2)*(∂A[2].*∂G[2] + CG.*∂G2[3] + ∂G[2].*∂A[2] + CA.*∂A2[3])
        return (∂β2,∂βγ,∂γ2)
    end
end

function ∂w2_derivative(g,k,a,u,β,γ,N,M;normalization)
    r = r_βγ(β,γ)
    w = wdomain(a,N)
    out = zeros(ComplexF64,M)
    w_g = w_derivative(g,k,a,u,β,γ,N,M,normalization)
    if k != 0
        Agprev = (k+r-1)*γ*a_β_γ(k,β,γ)/a_β_γ(k-1,β,γ)
        gprev=gmw(k-1,a,u,β,γ,N;normalization=normalization)
        w_gprev = w_derivative(gprev,k-1,a,u,β,γ,N,M,normalization)

        gprev = Agprev*gprev
        for i in 2:M
            Ag = ( β  - γ*exp(γ*log(w[i]))  + γ*k)
            ∇Ag = -(γ^2)*exp((γ-1)*log(w[i]))
            out[i] =  -(Ag*g[i] - gprev[i])/(w[i]^2) + 
            (∇Ag*g[i] + w_g[i]*Ag - w_gprev[i])/w[i]
        end
    else
        for i in 2:M
            Ag = ( β  - γ*exp(γ*log(w[i])) )/w[i]
            ∇Ag = -Ag/w[i] - (γ^2) * exp((γ-2)*log(w[i]))
            out[i] = Ag*w_g[i] + ∇Ag*g[i]
        end
    end
    return out
end

function ∂βγ∂w_derivative(g,k,a,u,β,γ,N,M=div(N,2)+1;normalization)
    r = r_βγ(β,γ)
    w = wdomain(a,N)
    ∂β = zeros(ComplexF64,M)
    ∂γ = zeros(ComplexF64,M)
    ∂βg,∂γg = βγ_derivative(g,k,a,u,β,γ,N,M,normalization)
    if k != 0
        ∇zg = _∇zg_A(β,γ,r)
        Ak  = a_β_γ(k,β,γ,r)
        ∇Ak = _∇gA(k,β,γ,r,Ak)
        Akm = a_β_γ(k-1,β,γ,r)
        ∇Akm    = _∇gA(k-1,β,γ,r,Akm)
        Agprev  = (k+r-1)*γ*Ak/Akm
        ∇Agprev = [ γ*Ak/Akm (k+r-1)*Ak/Akm]
        ∇Agprev  += (k+r-1)*γ*( ∇Ak/Akm .- (Ak/(Akm^2))*∇Akm ) 
        ∇Agprev  = ∇Agprev*∇zg
        gprev    = gmw(k-1,a,u,β,γ,N)
        ∂βγgprev = βγ_derivative(gprev,k-1,a,u,β,γ,N,M,normalization)
        ∂βgprev = ∇Agprev[1]*gprev + Agprev*∂βγgprev[1]
        ∂γgprev = ∇Agprev[2]*gprev + Agprev*∂βγgprev[2]
        for i in 2:M
            Ag    = ( β  - γ*exp(γ*log(w[i]))  + γ*k)
            ∂βAg  = 1
            ∂γAg  = k - (1+γ*log(w[i]))*exp(γ*log(w[i]))
            ∂β[i] = (∂βAg*g[i] + Ag*∂βg[i] - ∂βgprev[i])/w[i]
            ∂γ[i] = (∂γAg*g[i] + Ag*∂γg[i] - ∂γgprev[i])/w[i]
        end
    else
        for i in 2:M
            Ag    = ( β  - γ*exp(γ*log(w[i])) )/w[i]
            ∂βAg  = 1/w[i]
            ∂γAg  = -(1+γ*log(w[i]))*exp((γ-1)*log(w[i]))
            ∂β[i] = Ag*∂βg[i] + ∂βAg*g[i]
            ∂γ[i] = Ag*∂γg[i] + ∂γAg*g[i]
        end
    end
    return (∂β,∂γ)
end

function ∂a∂w_derivative(g,k,a,u,β,γ,N,M;normalization)
    r = r_βγ(β,γ)
    w = wdomain(a,N)
    ∂aw = wdomain(1,N)
    out = zeros(ComplexF64,M)
    ∂ag = a_derivative(g,k,a,u,β,γ,N,M,normalization)
    if k != 0
        Agprev = (k+r-1)*γ*a_β_γ(k,β,γ)/a_β_γ(k-1,β,γ)
        gprev = gmw(k-1,a,u,β,γ,N;normalization=normalization)
        ∂agprev = Agprev.*a_derivative(gprev,k-1,a,u,β,γ,N,M,normalization)
        gprev=Agprev.*gprev
        for i in 2:M
            Ag = ( β  - γ*exp(γ*log(w[i]))  + γ*k)
            ∂aAg = -∂aw[i]*(γ^2)*exp((γ-1)*log(w[i]))
            out[i] =  (∂aAg*g[i] + Ag*∂ag[i] - ∂agprev[i])/w[i] - ∂aw[i]*(Ag*g[i]-gprev[i])/(w[i]^2) 
        end
    else # OK
        for i in 2:M
        Ag = ( β  - γ*exp(γ*log(w[i])) )
        ∂aAg = -∂aw[i]*(γ^2)*exp((γ-1)*log(w[i]))
        out[i] = Ag*( ∂ag[i]/w[i] - ∂aw[i]*g[i]/(w[i]^2)) + ∂aAg*g[i]/w[i]
        end
    end
    return out
end

function ∂a2_derivative(g,k,a,u,β,γ,N,M;normalization)
  w = wdomain(1,N)[1:M]
  return ∂a∂w_derivative(g,k,a,u,β,γ,N,M;normalization).*w - (0.5/(g.a^2))*g + (0.5/(a))*a_derivative(g,k,a,u,β,γ,N,M,normalization)
end

function ∂a∂βγ_derivative(g,k,a,u,β,γ,N,M,normalization)
    w = wdomain(1,N)[1:M]
    g∂βγ = βγ_derivative(g,k,a,u,β,γ,N,M,normalization)
    ∂βγw = ∂βγ∂w_derivative(g,k,a,u,β,γ,N,M,normalization)
    ∂a∂β = w.*∂βγw[1] + (0.5/a) * g∂βγ[1]
    ∂a∂γ = w.*∂βγw[2] + (0.5/a) * g∂βγ[2]
    return (∂a∂β,∂a∂γ)
end

function ∂u∂auβγ_derivative(g,k,a,u,β,γ,N,M,normalization)
  d_k = -Complex.(0,wdomain(1,N,M))
  ∂uu = d_k .* u_derivative(g,k,a,u,β,γ,N,M,normalization)
  ∂ua = d_k .* a_derivative(g,k,a,u,β,γ,N,M,normalization)
  g∂βγ = βγ_derivative(g,k,a,u,β,γ,N,M,normalization)
  ∂uβ = d_k .* g∂βγ[1]
  ∂uγ = d_k .* g∂βγ[2]
  return (∂ua, ∂uu ,∂uβ ,∂uγ)
end

function ∂βγ∂βγ_derivative(g,k,a,u,β,γ,N,M,normalization)
    w = wdomain(a,N)[1:M]
    Cphase = tau_expo(u,N)*sqrt(a) # Keep phase
    partials=_calc_wavelet_∂βγ∂βγderiv(w[2:M],k,β,γ)
    ∂β2 = vcat(0,Cphase[2:M] .* partials[1])
    ∂βγ = vcat(0,Cphase[2:M] .* partials[2])
    ∂γ2 = vcat(0,Cphase[2:M] .* partials[3])
    return (∂β2,∂βγ,∂γ2)
end

#function ∇jacobian(g::GMF) # Derivative of vectorized jacobian
#    ∂uag,∂u2g,∂uβg,∂uγg = ∂u∂auβγ_derivative(g)
#    ∂u = vcat(∂uag,∂u2g,∂uβg,∂uγg)
#
#    ∂a2g = ∂a2_derivative(g)
#    ∂aβg, ∂aγg = ∂a∂βγ_derivative(g)
#    ∂a = vcat(∂a2g,∂uag,∂aβg,∂aγg)
#
#    ∂β2g,∂βγg,∂γ2g = ∂βγ∂βγ_derivative(g)
#    ∂β = vcat(∂aβg,∂uβg,∂β2g,∂βγg)
#    ∂γ = vcat(∂aγg,∂uγg,∂βγg,∂γ2g)
#    return (∂a,∂u,∂β,∂γ)
#
#end
#
#function ∇logabsdetjacobian(g::GMF)
#    # See matrix cook book 4.1.3 for complex derivative involving determinants
#    # dot product perform a conjugate transposition
#    # we have ∂det(X^HX)/∂X = det(X^H*X)(X^H*X)^-1X^H)^T
#    # or on ee.ic.ac.uk/hp/staff/dmb/matrix/calculus.html:
#    #  d(ln(det(X^HX))) = vec(conj(X) (X^Tconj(X)^(-1)))^T * vec(dX) + vec(X(X^HX)^(-1))^T  * vec(conj(dX))
#    J = jacobian(g)
#    ∇J = ∇jacobian(g)
##    JC1=conj(J)*(inv(transpose(J)*conj(J)))
##    JC2=J*(inv(J'*J))
##    ∂a = transpose(vec(JC1))*∇J[1] + transpose(vec(JC2))*conj(∇J[1])
##    ∂u = transpose(vec(JC1))*∇J[2] + transpose(vec(JC2))*conj(∇J[2])
##    ∂β = transpose(vec(JC1))*∇J[3] + transpose(vec(JC2))*conj(∇J[3])
##    ∂γ = transpose(vec(JC1))*∇J[4] + transpose(vec(JC2))*conj(∇J[4])
#    JC2=J*(inv(J'*J))
#    ∂a,∂u,∂β,∂γ = real(vec(JC2)'*hcat(∇J...))
#    # So Ok dont forget that we are computing log(sqrt(det(J'J)))=0.5log(det(J'J)) -> ∂ =0.5*2*vec(X^+)*conj(hcat(∇J...))
#    return (∂a,∂u,∂β,∂γ)
#end
