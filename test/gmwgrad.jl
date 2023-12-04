function apply∇jac(H,dv,N)
    M=div(N,2)+1
    out = zeros(ComplexF64,M)
    for i in 2:M
        out[i]=dv'*H[i:M:end,:]*dv
    end
    return out
end

N=2048
k=0
a=rand()*5+3
u=rand(0:2047)
β=rand()+1
γ=rand()+1
x=range(0,10,length=1024)
α=1.7
amp(x) = sqrt(sum(abs2,x))
M=div(N,2)
r=r_βγ(β,γ)
w = wdomain(a,N)
w_γ=2*exp.(γ .* log.(w))

# FIRST ORDER DIFFERENTIATION

remind_∇zg(dβ,dγ) = begin # OK
    gA = (β,γ) -> [r_βγ(β,γ);γ]
    dg = gA(β+dβ,γ+dγ)
    g = gA(β,γ)
    ∇zg = _∇zg_A(β,γ)
    dz = [dβ;dγ]
    return amp(dg-g-∇zg*dz)
end

remind_∇∇zg(dβ,dγ) = begin # OK
    d∇zg = _∇zg_A(β+dβ,γ+dγ)
    ∇zg = _∇zg_A(β,γ)
    ∇∇zg = _∇∇zg_A(β,γ)
    dz = [dβ;dγ]
    return amp(vec(d∇zg) - vec(∇zg) - ∇∇zg*dz)
end

remind_∇gA(dβ,dγ)= begin # OK
    dz=[dβ;dγ]
    ∇zg = _∇zg_A(β,γ)
    dg=∇zg*dz
    dA=a_β_γ(k,β+dβ,γ+dγ)
    A=a_β_γ(k,β,γ)
    ∇gA = _∇gA(k,β,γ)
    return dA - A - (∇gA*dg)[1]
end

remind_∇∇gA(dβ,dγ) = begin # OK
    dz=[dβ;dγ]
    ∇zg = _∇zg_A(β,γ)
    dg=∇zg*dz
    d∇gA = _∇gA(k,β+dβ,γ+dγ)
    ∇gA = _∇gA(k,β,γ)
    ∇∇gA = _∇∇gA(k,β,γ)
    return vec(d∇gA) - vec(∇gA) - ∇∇gA*dg
end

remind_aβγ_deriv(dβ,dγ) = begin # OK
    adβγ=a_β_γ(k,β+dβ,γ+dγ)
    aβγ = a_β_γ(k,β,γ)
    aβγ_deriv = a_β_γ_deriv(k,β,γ)
    return adβγ - aβγ - dβ*aβγ_deriv[1] - dγ*aβγ_deriv[2]
end

remind_a_β_γ_∂βγ∂βγderiv(dβ,dγ)= begin # ok
    d∂βγ = a_β_γ_deriv(k,β+dβ,γ+dγ)
    ∂βγ  = a_β_γ_deriv(k,β,γ)
    ∂∂βγ = a_β_γ_2deriv(k,β,γ)
    r_β = d∂βγ[1] - ∂βγ[1] - ∂∂βγ[1]*dβ - ∂∂βγ[2]*dγ
    r_γ = d∂βγ[2] - ∂βγ[2] - ∂∂βγ[2]*dβ - ∂∂βγ[3]*dγ
    return (amp(r_β),amp(r_γ))
end

remind_a_β_γ_2deriv(dβ,dγ) = begin # OK
    adβγ=a_β_γ(k,β+dβ,γ+dγ)
    aβγ = a_β_γ(k,β,γ)
    ∂β,∂γ= a_β_γ_deriv(k,β,γ)
    ∂β2,∂βγ,∂γ2=a_β_γ_2deriv(k,β,γ)
    return adβγ - aβγ - dβ*∂β - dγ*∂γ - 0.5*((dβ^2)*∂β2 + dγ^2*∂γ2) - dβ*dγ*∂βγ
end

remind_laguerre_deriv(dα,dx) = begin #OK
    ldαdx=laguerrepoly(k,α+dα,x.+dx)
    l=laguerrepoly(k,α,x)
    a_deriv = α_derivative_laguerre(k,α,x)
    x_deriv = x_derivative_laguerre(k,α,x)
    return ldαdx - l - dα*a_deriv - dx*x_deriv
end

remind_derivative_laguerre(dα,dx) = begin # OK
    α_deriv(dα,dx) = α_derivative_laguerre(k,α+dα,x.+dx)
    x_deriv(dα,dx) = x_derivative_laguerre(k,α+dα,x.+dx)
    x2_deriv,αx_deriv,α2_deriv = x2_αx_α2_derivative_laguerre(k,α,x)
    r_α = α_deriv(dα,dx) - α_deriv(0,0) - dα*α2_deriv - dx*αx_deriv # ok
    r_x = x_deriv(dα,dx) - x_deriv(0,0) - dx*x2_deriv - dα*αx_deriv # ok
    return (r_α,r_x)
end

remind_laguerre_2deriv(dα,dx) = begin # OK 
    ldαdx=laguerrepoly(k,α+dα,x.+dx)
    l=laguerrepoly(k,α,x)
    α_deriv = α_derivative_laguerre(k,α,x)
    x_deriv = x_derivative_laguerre(k,α,x)
    α2_deriv,αx_deriv,x2_deriv = x2_αx_α2_derivative_laguerre(k,α,x)
    return ldαdx - l - dα*α_deriv - dx*x_deriv -0.5*((dα^2)*α2_deriv + (dx^2)*x2_deriv) - dα*dx*αx_deriv
end

remind_lpoly_wave_deriv(dβ,dγ) = begin # OK
    dL=lpoly_wave(w[2:M],k,β+dβ,γ+dγ)
    L=lpoly_wave(w[2:M],k,β,γ)
    ∂L=lpoly_wave_βγderiv(w[2:M],k,β,γ)
    dv = [dβ dγ]
    return amp(dL - L - ∂L[1]*dv[1] - ∂L[2]*dv[2])
end

# Jacobian matrix testing

remind_∇zgL(dβ,dγ) = begin # OK
    gL= (β,γ) -> (r_βγ(β,γ),2*exp.(γ .* log.(w[2:end])))
    dg = gL(β+dβ,γ+dγ)
    g = gL(β,γ)
    ∇zg = _∇zg_L(w[2:end],β,γ)
    dz = [dβ dγ]
    r_β=dg[1] - g[1] - ((∇zg[1,1:2]')*dz')[1]
    r_γ=dg[2] - g[2] - (∇zg[2:end,1:2])*dz'
    return (amp(r_β),amp(r_γ))
end

remind_∇∇zgL(dβ,dγ) = begin # OK
    ∇zg = _∇zg_L(w[2:end],β,γ)[:]
    d∇zg= _∇zg_L(w[2:end],β+dβ,γ+dγ)[:]
    ∇∇zg = _∇∇zg_L(w[2:end],β,γ)
    dz = [dβ dγ]
    r=d∇zg - ∇zg - (∇∇zg*dz')
    return amp(r)
end

remind_lpoly_wave_∂βγ∂βγderiv(dβ,dγ) = begin # OK for k<=1
    d∂L = lpoly_wave_βγderiv(w[2:M],k,β+dβ,γ+dγ)
    ∂L = lpoly_wave_βγderiv(w[2:M],k,β,γ)
    ∂β2,∂βγ,∂γ2 = lpoly_wave_∂βγ∂βγderiv(w[2:M],k,β,γ)
    dv = [ dβ dγ ]
    r_β=d∂L[1] - ∂L[1] - ∂β2*dβ - ∂βγ*dγ
    r_γ=d∂L[2] - ∂L[2] - ∂γ2*dγ - ∂βγ*dβ 
    return (amp(r_β),amp(r_γ))
    # For K<=2 OK
    # FOR K>2 OK~, the approximation is not that good O(amp([dβ,dγ]))
    # BUT: See next test with full approx of lpoly_wave
    # Bad approx with ∂γ2 for k>2
end

remind_lpoly_wave_2deriv(dβ,dγ) = begin # BIG OK, the approximation is better than wave_deriv and for k=0,1,2,3 ( so for big k also ?)
    # It confirms that lpoly_wave_∂βγ∂βγderiv works also for k>2
    dv = [ dβ dγ]
    wh=w[2:end]
    dl = lpoly_wave(wh,k,β+dβ,γ+dγ)
    l = lpoly_wave(wh,k,β,γ)
    ∇l = lpoly_wave_βγderiv(wh,k,β,γ)
    ∇∇l = lpoly_wave_∂βγ∂βγderiv(wh,k,β,γ)
    r = dl - l - ∇l[1]*dβ - ∇l[2]*dγ - 0.5*∇∇l[1]*dβ^2 - 0.5*∇∇l[3]*dγ^2 - ∇∇l[2]*dβ*dγ
    return amp(r)
end 

remind_core_exp_deriv(dβ,dγ) = begin # OK
    wh=w[2:end]
    dcore =core_exp(wh,β+dβ,γ+dγ)
    core= core_exp(wh,β,γ)
    ∂core  = core_exp_βγderiv(wh,β,γ)
    return amp(dcore - core -∂core[1]*dβ -∂core[2]*dγ)
end

remind_core_exp_∂βγ∂βγderiv(dβ,dγ) = begin # OK
    wh = w[2:end]
    d∂core  = core_exp_βγderiv(wh,β+dβ,γ+dγ)
    ∂core = core_exp_βγderiv(wh,β,γ)
    ∂∂core = core_exp_∂βγ∂βγderiv(wh,β,γ)
    r_β=d∂core[1]-∂core[1]-∂∂core[1]*dβ - ∂∂core[2]*dγ
    r_γ=d∂core[2]-∂core[2]-∂∂core[2]*dβ - ∂∂core[3]*dγ
    return (amp(r_β),amp(r_γ))
end

remind_core_exp_2deriv(dβ,dγ) = begin # Ok, 100 times better
    wh=w[2:end]
    dcore=core_exp(wh,β+dβ,γ+dγ)
    core=core_exp(wh,β,γ)
    ∂core = core_exp_βγderiv(wh,β,γ)
    ∂∂core = core_exp_∂βγ∂βγderiv(wh,β,γ)
    return amp(dcore - core -∂core[1]*dβ -∂core[2]*dγ - 0.5*(∂∂core[1]*dβ^2+∂∂core[3]*dγ^2) - dβ*dγ*∂∂core[2])
end

remind_calc_wavelet_deriv(dβ,dγ) = begin # OK
    wh = w[2:end]
    dc_w = _calc_wavelet(wh,k,β+dβ,γ+dγ)
    c_w = _calc_wavelet(wh,k,β,γ)
    ∂c_w = _calc_wavelet_βγderiv(wh,k,β,γ)
    return amp(dc_w-c_w -∂c_w[1]*dβ-∂c_w[2]*dγ)
end

remind_calc_wavelet_∂βγ∂βγderiv(dβ,dγ) = begin # Not Ok for k>0
    wh = w[2:end]
    d∂c_w = _calc_wavelet_βγderiv(wh,k,β+dβ,γ+dγ)
    ∂c_w = _calc_wavelet_βγderiv(wh,k,β,γ)
    ∂∂c_w = _calc_wavelet_∂βγ∂βγderiv(wh,k,β,γ)
    r_β=d∂c_w[1]-∂c_w[1]-∂∂c_w[1]*dβ-∂∂c_w[2]*dγ
    r_γ=d∂c_w[2]-∂c_w[2]-∂∂c_w[2]*dβ-∂∂c_w[3]*dγ
    return (r_β,r_γ,amp(r_β),amp(r_γ))
end

remind_calc_wavelet_2deriv(dβ,dγ) = begin # OK
    wh = w[2:end]
    dc_w = _calc_wavelet(wh,k,β+dβ,γ+dγ)
    c_w = _calc_wavelet(wh,k,β,γ)
    ∂c_w = _calc_wavelet_βγderiv(wh,k,β,γ)
    ∂∂c_w = _calc_wavelet_∂βγ∂βγderiv(wh,k,β,γ)
    return amp(dc_w-c_w -∂c_w[1]*dβ-∂c_w[2]*dγ - 0.5*(∂∂c_w[1]*dβ^2 +∂∂c_w[3]*dγ^2) - dβ*dγ*∂∂c_w[2])
end

remind_βγderiv(dβ,dγ) = begin # OK
    gdβγ=GMF(k,a,u,β+dβ,γ+dγ,N)
    g=GMF(k,a,u,β,γ,N)
    ∂β,∂γ = βγ_derivative(g)
    return amp(gdβγ.gmw - (g.gmw + dβ*∂β + dγ*∂γ))
end

remind_∂βγ∂βγderiv(dβ,dγ) = begin # Not Ok for k>0
    dg=GMF(k,a,u,β+dβ,γ+dγ,N)
    g=GMF(k,a,u,β,γ,N)
    d∂βγ = βγ_derivative(dg)
    ∂βγ = βγ_derivative(g)
    ∂∂βγ = ∂βγ∂βγ_derivative(g)
    r_β = d∂βγ[1] - ∂βγ[1] - ∂∂βγ[1]*dβ - ∂∂βγ[2]*dγ
    r_γ = d∂βγ[2] - ∂βγ[2] - ∂∂βγ[2]*dβ - ∂∂βγ[3]*dγ
    return (amp(r_β),amp(r_γ))
end

remind_βγ2deriv(dβ,dγ) = begin # OK
    gdβγ=GMF(k,a,u,β+dβ,γ+dγ,N)
    g=GMF(k,a,u,β,γ,N)
    dv = [dβ;dγ]
    ∂β,∂γ = βγ_derivative(g)
    ∇g=[∂β ∂γ]
    ∂β2, ∂βγ, ∂γ2 = ∂βγ∂βγ_derivative(g)
    return amp(gdβγ.gmw - g.gmw - ∇g*dv - 0.5*(∂β2*dβ^2+∂γ2*dγ^2) - dβ*dγ*∂βγ)
end

remind_uderiv(du) = begin # OK
    gdu=GMF(k,a,u+du,β,γ,N)
    g=GMF(k,a,u,β,γ,N)
    return amp(gdu.gmw - (g.gmw + du * u_derivative(g)))
end

remind_∂u∂auβγ_derivative(du,da,dβ,dγ) = begin # OK
    g   = GMF(k,a,u,β,γ,N)
    d∂u = u_derivative(GMF(k,a+da,u+du,β+dβ,γ+dγ,N))
    ∂u = u_derivative(GMF(k,a,u,β,γ,N))
    ∂ua, ∂uu ,∂uβ ,∂uγ =∂u∂auβγ_derivative(g::GMF)
    J = [ ∂ua ∂uu ∂uβ ∂uγ ] 
    dv = [ da du dβ dγ]
    return amp(d∂u - ∂u - J*dv')
end

remind_aderiv(da) = begin # OK
    gda=GMF(k,a+da,u,β,γ,N)
    g=GMF(k,a,u,β,γ,N)
    return amp(gda.gmw - (g.gmw + da * a_derivative(g)))
end

remind_a2deriv(da) = begin  # OK
    dg=GMF(k,a+da,u,β,γ,N)
    g=GMF(k,a,u,β,γ,N)
    d∂a = a_derivative(dg)
    ∂a = a_derivative(g)
    ∂a2 = ∂a2_derivative(g)
    r∂a=amp(d∂a - ∂a - ∂a2*da) # Not very good
    rag=amp(dg.gmw - g.gmw -∂a*da - 0.5*(da^2)*∂a2) # Better than remind_aderiv
    return (r∂a,rag)
end

remind_∂a∂βγ_derivative(dβ,dγ) = begin # Ok for k=0
    dg=GMF(k,a,u,β+dβ,γ+dγ,N)
    g=GMF(k,a,u,β,γ,N)
    d∂a=a_derivative(dg)
    ∂a=a_derivative(g)
    ∂β,∂γ=∂a∂βγ_derivative(g)
    dv = [dβ dγ]'
    return amp(d∂a - ∂a - [∂β ∂γ]*dv)
end

remind_jacobian(da,du,dβ,dγ) = begin # OK for k=0
    dg = GMF(k,a+da,u+du,β+dβ,γ+dγ,N)
    g  = GMF(k,a,u,β,γ,N)
    jac =jacobian(g)
    dv = [ da ,du ,dβ ,dγ]
    return amp(dg.gmw -g.gmw - jac*dv)
end

remind_∇jacobian(da,du,dβ,dγ) = begin # OK for k=0
    dg=GMF(k,a+da,u+du,β+dβ,γ+dγ,N)
    g=GMF(k,a,u,β,γ,N)
    djac=jacobian(dg)
    jac=jacobian(g)
    hes=∇jacobian(g)
    dv=[da,du,dβ,dγ]
    return vec(djac) - vec(jac) - hcat(hes...)*dv
end

remind_hessian(da,du,dβ,dγ)= begin # OK for k=0
    dg = GMF(k,a+da,u+du,β+dβ,γ+dγ,N)
    g  = GMF(k,a,u,β,γ,N)
    jac =jacobian(g)
    hes = ∇jacobian(g)
    dv = [ da ,du ,dβ ,dγ]
    return amp(dg.gmw -g.gmw - jac*dv -0.5*apply∇jac(hcat(hes...),dv,N))
end

remind_logabsdet(da,du,dβ,dγ) = begin # Ok for k=0, even though I dont know why
    dg=GMF(k,a+da,u+du,β+dβ,γ+dγ,N)
    g=GMF(k,a,u,β,γ,N)
    logd = logabsdetjacobian(g)
    dlogd=logabsdetjacobian(dg)
    ∂logd = [ ∇logabsdetjacobian(g)...]'
    dv=[da,du,dβ,dγ]
    return amp(dlogd-logd-∂logd*dv)
end 

