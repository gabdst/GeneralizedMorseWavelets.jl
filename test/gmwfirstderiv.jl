push!(LOAD_PATH,"../")
using GMW
import GMW: a_derivative,β_derivative,γ_derivative,u_derivative,logA_peak_βγ_deriv,A_peak
using Test

N=2048
k=0
L=5
grid=gmw_grid(1,3,Int(log2(N)-1),1)
ϵ=1e-3
c=10

function remind_derivative(f,∂f,x,ϵ)
  df=f(x+ϵ)-f(x)
  r=ϵ*∂f(x)
  return (df,r)
end

# We only pass the ϵ test, should get the ϵ^2/2 tolerance test.
get_gmw(a,u,β,γ)=gmw(k,a,u,β,γ,N,:peak)
@testset "Derivatives w peak normalization" for pg in grid
  a,u,β,γ=pg
  p=(k,pg...,N,:peak)
  r_β=remind_derivative(β->log(A_peak(k,β,γ)), β->logA_peak_βγ_deriv(β,γ)[1], β, ϵ)
  r_γ=remind_derivative(γ->log(A_peak(k,β,γ)), γ->logA_peak_βγ_deriv(β,γ)[2], γ, ϵ)
  @test isapprox(r_β...,atol=c*ϵ^2)
  @test isapprox(r_γ...,atol=c*ϵ^2)
  r_a=remind_derivative(a->get_gmw(a,u,β,γ),a->a_derivative(get_gmw(a,u,β,γ),p...),a,ϵ)
  r_u=remind_derivative(u->get_gmw(a,u,β,γ),u->u_derivative(get_gmw(a,u,β,γ),p...),u,ϵ) ### BUG when u=0 and u+du !=0
  r_β=remind_derivative(β->get_gmw(a,u,β,γ),β->β_derivative(get_gmw(a,u,β,γ),p...),β,ϵ)
  r_γ=remind_derivative(γ->get_gmw(a,u,β,γ),γ->γ_derivative(get_gmw(a,u,β,γ),p...),γ,ϵ)
  @test isapprox(r_a...,atol=c*ϵ^2) 
  @test isapprox(r_β...,atol=c*ϵ^2)
  @test isapprox(r_γ...,atol=c*ϵ^2)
  @test isapprox(r_u...,atol=c*ϵ^2)
end

get_gmw(a,u,β,γ)=gmw(k,a,u,β,γ,N,:L2)
@testset "Derivatives w\\ L2 normalization" for pg in grid
  a,u,β,γ=pg
  p=(k,pg...,N,:L2)
  r_a=remind_derivative(a->get_gmw(a,u,β,γ),a->a_derivative(get_gmw(a,u,β,γ),p...),a,ϵ)
  r_u=remind_derivative(u->get_gmw(a,u,β,γ),u->u_derivative(get_gmw(a,u,β,γ),p...),u,ϵ)
  r_β=remind_derivative(β->get_gmw(a,u,β,γ),β->β_derivative(get_gmw(a,u,β,γ),p...),β,ϵ)
  r_γ=remind_derivative(γ->get_gmw(a,u,β,γ),γ->γ_derivative(get_gmw(a,u,β,γ),p...),γ,ϵ)
  @test isapprox(r_a...,atol=c*ϵ^2) 
  @test isapprox(r_β...,atol=c*ϵ^2)
  @test isapprox(r_γ...,atol=c*ϵ^2)
  @test isapprox(r_u...,atol=c*ϵ^2)
end
