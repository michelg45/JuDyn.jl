function gauss_points(ndim::Int, npoints::Int)
    if ndim == 1
        if npoints == 1
            x = [0.0]
            w = [2.0]
        elseif npoints == 2
            x = 1.0/sqrt(3.0)*[-1.0; 1.0]
            w = ones(2)
        elseif npoints == 3
            x = sqrt(3.0/5.0)*[0.0; -1.0; 1.0]
            w = 1.0/9.0*[8.0; 5.0; 5.0]
        elseif npoints == 4
            x = [-sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0));
                  sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0));
                -sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0));
                sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0))]
            w = 1.0/36.0*(18.0*ones(4)+sqrt(30.0)*[1.0; 1.0; -1.0; -1.0])
        end
        elseif ndim == 2
            if npoints == 1
                x = [[0.0; 0.0 ]]
                w = [4.0]
            elseif npoints == 4;
                x = 1.0/sqrt(3.0)*[[-ones(2)]; [[1.0; -1.0;]]; [[-1.0; 1.0;]]; [ones(2)]] 
                w = ones(4)
            elseif npoints == 7

                """w = [1.0/40.0, 1.0/15.0, 1.0/40.0, 1.0/15.0, 1.0/40.0, 1.0/15.0, 9.0/40.0]
                x = [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [0.5, 0.5], [0.0, 1.0], [0.0, 0.5], [1.0/3.0, 1.0/3.0]]"""

                w = zeros(7)
                w[1] = 0.225
                w[2:4] .= 0.132394152788506
                w[5:7] .= 0.125939180544827

                w = 0.5*w

                alpha_1 = 0.059715871789770
                beta_1 =  0.470142064105115
                alpha_2 = 0.797426985353087
                beta_2 =  0.101286507323456
                
                tr_coord =  Vector{Vector}(undef,7)

                tr_coord[1] = [1.0/3.0, 1.0/3.0,1.0/3.0]
                tr_coord[2]  = [alpha_1,beta_1,beta_1]
                tr_coord[3] = [beta_1,alpha_1,beta_1]
                tr_coord[4] = [beta_1,beta_1,alpha_1]
                tr_coord[5]  = [alpha_2,beta_2,beta_2]
                tr_coord[6] = [beta_2,alpha_2,beta_2]
                tr_coord[7] = [beta_2,beta_2,alpha_2]

                vertices = zeros(2,3)
                vertices[1,2] = 1.0
                vertices[2,3] = 1.0

                x = Vector{Vector}(undef,7)
                for i = 1:7
                    x[i] = vertices*tr_coord[i]
                end


            elseif npoints == 9
                x = sqrt(3.0/5.0)*[[-ones(2)]; [[0.0; -1.0]];[[1.0; -1.0]];
                                    [[-1.0; 0.0]]; [zeros(2)]; [[1.0; 0.0]];
                                    [[-1.0; 1.0]]; [[0.0; 1.0]]; [ones(2)]]
                w = 1/81.0*[25.0; 40.0; 25.0; 40.0; 64.0; 40.0; 25.0; 40.0; 25.0]
        end
    end
    return x, w
end