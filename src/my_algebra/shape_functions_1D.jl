
"""
     shape_functions_1D
"""
function shape_functions_1D(Nnodes,xi)
    if Nnodes == 2
         F = [0.5*(1.0-xi); 0.5*(1.0+xi)]
         DF = [-0.5; 0.5]
    elseif Nnodes == 3
         F = [0.5*xi*(xi-1.0); 1.0-xi*xi; 0.5*xi*(1.0+xi)]
         DF = [xi-0.5; -2.0*xi; xi+0.5]
    elseif Nnodes == 8
         eta = xi[1]
         mu = xi[2]
         F = 0.5625*[-(1.0-xi)*(1.0/3.0+xi)*(1.0/3.0-xi);
         3.0*(1.0-xi)*(1.0+xi)*(1.0/3.0-xi); 3.0(1.0-xi)*(1.0+xi)*(1.0/3.0+xi);
         -(1.0+xi)*(1.0/3.0+xi)*(1.0/3.0-xi)]
         DF = [-1.6875*xi*xi+1.125*xi+0.0625; 5.0625*xi*xi-1.125*xi-1.6875; 
         -5.0625*xi*xi-1.125*xi+1.6875; 1.6875*xi*xi+1.125*xi-0.0625]
    end
    return F, DF
 end