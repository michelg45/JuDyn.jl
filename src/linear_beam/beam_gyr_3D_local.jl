"""
    beam_gyr_3D_local

        function constructing the gyroscopic mass kernel of a beam element, based on a linear interpolation of the deflection.

            Input:

                l::Float64      element length.
                m::Float64      mass per unit length.

            Output:

                S_el::Vector{Matrix}  set of (12x12) matrices forming the gyroscopic matrix of the element: G_el = sum(omega[i]*S_el[i]) 


            Calling sequence: 
                S_el = beam_gyr_3D_local(l,m)

"""
function beam_gyr_3D_local(l::Float64,m::Float64,degree::String)
    mass_kernel = bar_mass_2D(l,m,degree)
    if degree == "quadratic"
        ndim = 18
        nd = 3
        iloc = [1; 7; 13; 2; 8; 14; 3; 9; 15]
    else
        ndim = 12
        nd = 2
        iloc = [1; 7; 2; 8; 3; 9]
    end
    M_gyr = zeros(ndim,ndim)
     
    M_gyr[iloc,iloc] = cat(mass_kernel, mass_kernel, mass_kernel; dims=(1,2))
    S_el = Vector{Matrix{Float64}}(undef,3)
    locel = Vector{Vector{Int}}(undef,3)
    for i = 1:3
        S_el[i] = zeros(ndim,ndim)
        locel[i] = zeros(Int,nd)
        locel[i][1] = i 
        locel[i][2] = i+6
        nd == 3 && (locel[i][nd] = i+12)
    end
    for i = 1:3
        j = i+1
        j > 3 && (j = mod(j,3))
        k = i+2
        k > 3 && (k = mod(k,3))
        S_el[i][locel[j],locel[k]] = M_gyr[locel[i],locel[i]]
        S_el[i][locel[k],locel[j]] = - M_gyr[locel[i],locel[i]]
    end
            

    return S_el
end

function beam_gyr_3D_local(l::Float64,m::Float64)
    degree = "linear"
    return S_el = beam_gyr_3D_local(l,m,degree)
end