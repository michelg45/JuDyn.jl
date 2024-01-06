"""
    waterfall_plot(...)
        function to create a 3-D waterfall plot
        calling sequence : 
        waterfall_plot(yrange::Matrix{Float64}, x::Vector{Float64}, z::Array{Float64,2},
        labels::Vector{String},color::String,perm_yz::Bool)
"""
function waterfall_plot(yrange::Matrix{Float64}, x::Vector{Float64}, z::Array{Float64,2}, labels::Vector{String},color::String,perm_yz::Bool)
    

    
    
    ny = size(yrange,1)
    nx = size(x,1)
    nz = size(z,2)
    nx != nz &&(error("waterfall plot: x and z[1,:] data must have same size"))
    ny != size(z,1) &&(error("waterfall plot: y and z[:,1] data must have same size"))
    y = zeros(nz)
    
    wplot = plot(xlabel = labels[1], ylabel = labels[2], zlabel = labels[3], legend =:none)
    for i = 1:ny
        y .= yrange[i]
        wplot = plot!(x,y,z[i,:],color=color)
    end

    return wplot
end

