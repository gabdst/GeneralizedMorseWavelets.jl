module GeneralizedMorseWavelets
    include("laguerrepoly.jl")
    include("gmw.jl")
    include("gmwfirstderiv.jl")
    include("gmwscndderiv.jl")
    include("gmwADrules.jl")
    export w_psi,peak_w,peak_n,a_for_wpsi,mappeakgrid,duration,gmw_grid,gmw,gmw!
end
