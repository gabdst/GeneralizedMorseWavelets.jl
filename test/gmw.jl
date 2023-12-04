push!(LOAD_PATH,"../")
using GMW
using FFTW
using LinearAlgebra
using Test

k=0
a=1
u=0
β=1
γ=3
N=1024

@testset for normalization in [:L2,:peak]

  g=gmw(k,a,u,β,γ,N,normalization)
  h=copy(g)
  gmw!(g,k,a,u,β,γ,N,normalization)
  @test all(h .== g)

  g=gmw(k,a,u,β,γ,N,normalization)
  h=copy(g)
  gmw!(g,k,a,u,β,β,N,normalization) # β != γ
  @test any(h .!= g)
end

@testset "Normalization Test" begin
  normalization=:L2
  g=gmw(k,a,u,β,γ,N,normalization)
  @test isapprox(norm(irfft(g,N)),sqrt(2))

  normalization=:peak
  g=gmw(k,a,u,β,γ,N,normalization)
  w_psi=peak_n(a,β,γ,N)
  w_psi_g = argmax(abs.(g))-1
  @test isapprox(w_psi,w_psi_g,atol=1)
end

@testset "Translation Test" for u in [0,128,512,-511]
  normalization=:peak
  g=gmw(k,a,u,β,γ,N,normalization)
  θ=g./abs.(g)
  θ[1]=1
  x=vcat(1,zeros(N-1))
  x_shift_θ=irfft(rfft(x) .* θ,N)
  @test isapprox(circshift(x,u),irfft(rfft(x).*θ,N))
end
