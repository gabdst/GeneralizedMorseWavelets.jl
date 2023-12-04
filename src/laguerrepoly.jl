import Base:size,getindex,setindex!

function laguerrepoly(domain::AbstractArray,k::Integer,α::Real)
    (k<0) && throw(error("Generalized Laguerre polynomials are defined for k >0"))
    if k==0
        poly = ones(Float64,size(domain))
    elseif k==1
        poly = Array{Float64}((1+α .- domain))
    else
        poly = _generate_k_laguerre(k,α,domain)
    end
    return poly
end

function _next_laguerre(domain::AbstractArray,L_km::Array{Float64},L_k::Array{Float64},α::Real,k::Integer)
	L_knext = ((2*k+1+α .- domain).*L_k .- (k+α).*L_km)./(k+1)
	return L_knext
end

function _generate_k_laguerre(domain::AbstractArray,k,α)
    L_km = laguerrepoly(0,α,domain)
    L_k  = laguerrepoly(1,α,domain)
    i=1
    while i<k
        buf = _next_laguerre(domain,L_km,L_k,α,i)
        L_km = L_k
        L_k = buf
        i+=1
    end
    return L_k
end

function x_derivative_laguerre(domain::AbstractArray,k::Integer,α::Real)
    out = zeros(length(domain))
    if k == 0
        return out
    else
        return -laguerrepoly(domain,k-1,α+1)
    end
end

function x2_derivative_laguerre(domain::AbstractArray,k::Integer,α::Real)
    out = zeros(length(domain))
    if k > 1
        return laguerrepoly(domain,k-2,α+2,)
    else 
        return out
    end
end

function x2_αx_α2_derivative_laguerre(domain::AbstractArray,k::Integer,α::Real)
    if k == 0 || k == 1
        N=length(domain)
        return (zeros(N),zeros(N),zeros(N))
    end
    ∂x2 = x2_derivative_laguerre(domain,k,α)
    ∂xα = -α_derivative_laguerre(domain,k-1,α+1)
    ∂α2 = α2_derivative_laguerre(domain,k,α)
    return (∂x2,∂xα,∂α2)
end

function α_derivative_laguerre(domain::AbstractArray,k::Integer,α::Real)
    out = zeros(length(domain))
    if k == 0
        return out
    end
    L_km = laguerrepoly(domain,0,α)
    out = L_km/k
    if k == 1
        return out
    end
    L_k  = laguerrepoly(domain,1,α)
    out += L_k/(k-1)
    i=1
    while i<k-1
        buf = _next_laguerre(domain,L_km,L_k,α,i)
        L_km = L_k
        L_k = buf
        i+=1
        out += L_k/(k-i)
    end
    return out
end

_cumul(k,j) = begin
    c=0
    for i in (j+1):(k-1)
        c+= 1/((k-i)*(i-j))
    end
    return c
end

function α2_derivative_laguerre(domain::AbstractArray,k::Integer,α::Real)
    # ∂α{∑(i=0 to k-1) L(α,i)/(k-i) } = ∑(i=1 to k-1) ∑(j=0 to i-1) L(α,j)/( (k-i) * (i-j) )
    #                                 = ∑(j=0 to k-2) L(α,j) ∑(i=j+1 to k-1) 1/( (k-i)*(i-j) )
    #                                 = ∑(j=0 to k-2) L(α,j) * C(j)
    
    out = zeros(length(domain))
    if k == 0 || k == 1
        return out
    end
    L_km = laguerrepoly(domain,0,α)
    C = _cumul(k,0)
    out = L_km*C
    if k == 2
        return out
    end
    L_k  = laguerrepoly(domain,1,α)
    out+=L_k*_cumul(k,1)
    i=1
    while i<k-2
        buf = _next_laguerre(domain,L_km,L_k,α,i)
        L_km = L_k
        L_k = buf
        C=_cumul(k,i)
        i+=1
        out += L_k*C
    end
    return out
end
