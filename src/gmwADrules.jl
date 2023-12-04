using ChainRulesCore

function ChainRulesCore.rrule(::typeof(gmw),k::Integer,a::Real,u::Real,β::Real,γ::Real,N::Integer,normalization::Symbol)
  g = gmw(k,a,u,β,γ,N,normalization)
  function gmw_pullback(Δg) 
    ∂k = NoTangent()
    ∂a = @thunk( begin 
      ∂a = a_derivative(g,k,a,u,β,γ,N,normalization)
      ∂a = real(dot(Δg,∂a)) #/sqrt(dot(∂a,∂a))) 
    end)

    ∂u = @thunk( begin 
      ∂u = u_derivative(g,k,a,u,β,γ,N,normalization)
      ∂u = real(dot(Δg,∂u)) #/sqrt(dot(∂u,∂u)))
    end)

    ∂β=@thunk( begin 
      ∂β = β_derivative(g,k,a,u,β,γ,N,normalization)
      ∂β = real(dot(Δg,∂β)) #/sqrt(dot(∂β,∂β)))
    end)

    ∂γ = @thunk( begin 
      ∂γ = γ_derivative(g,k,a,u,β,γ,N,normalization)
      ∂γ = real(dot(Δg,∂γ)) #/sqrt(dot(∂γ,∂γ)))
    end)
    return (NoTangent(),∂k,∂a,∂u,∂β,∂γ,NoTangent())
  end
  return g, gmw_pullback
end


function ChainRulesCore.frule(deltas,::typeof(gmw),k::Integer,a::Real,u::Real,β::Real,γ::Real,N::Integer,normalization::Symbol)
  _,_,Δa,Δu,Δβ,Δγ,_,_ =deltas
  g = gmw(k,a,u,β,γ,N,normalization)
  ∂a = a_derivative(g,k,a,u,β,γ,N,normalization)
  ∂u = u_derivative(g,k,a,u,β,γ,N,normalization)
  ∂β = β_derivative(g,k,a,u,β,γ,N,normalization)
  ∂γ = γ_derivative(g,k,a,u,β,γ,N,normalization)
  Δg = Δa*∂a + Δu*∂u + Δβ*∂β + Δγ*∂γ
  return (g,Δg) 
end
## frule !?
