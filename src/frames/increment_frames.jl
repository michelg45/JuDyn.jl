"""
    increment_frames

        function incrementing theta node reference configurations at end of step
        in array current_frames:  x := x + Dx, psi := RV3(psi,Dpsi)


"""
    function increment_frames(Dy::Vector,y::Vector,Nnodes::Int)



            cf = Main.Frames.current_frames



        for inode=1:Nnodes

            inv_loc=copy(cf[inode].inv_loc)

            global inv_loc_x = inv_loc[1:3]
            global inv_loc_theta = inv_loc[4:6]

            global min_theta = findmin(inv_loc_theta)[1]
            global max_theta = findmax(inv_loc_theta)[1]
            global min_x = findmin(inv_loc_x)[1]
            global max_x = findmax(inv_loc_x)[1]

            if  cf[inode].type == "frame" 
                if  min_theta > 0
                    Dpsi = RV3(Dy[inv_loc_theta])
                else
                    Dpsi = RV3()
                    if max_theta > 0
                        for ic=1:3
                            inv_loc_theta[ic] > 0 &&  (Dpsi[ic] = Dy[inv_loc_theta[ic]])
                        end               
                    end
                end
                cf[inode].psi = RV3(cf[inode].psi, Dpsi)

                if  min_theta > 0
                    y[inv_loc_theta] = cf[inode].psi.v
                else
                    if max_theta > 0
                        for ic=1:3
                            inv_loc_theta[ic] > 0 &&  (y[inv_loc_theta[ic]]  = cf[inode].psi[ic])
                        end               
                    end
                end
                
            end
                
            if  min_x > 0
                Dx = Vec3(Dy[inv_loc_x])
            else
                Dx = Vec3()
                max_x = findmax(inv_loc_x)[1]
                if max_x > 0
                    for ic=1:3
                        inv_loc_x[ic] > 0 &&  (Dx[ic] = Dy[inv_loc_x[ic]])
                    end               
                end
            end

            cf[inode].x =  cf[inode].x + Dx

            if  min_x > 0
                y[inv_loc_x] = cf[inode].x.v
            else
                if max_x > 0
                    for ic=1:3
                            inv_loc_x[ic] > 0 &&  (y[inv_loc_x[ic]]  = cf[inode].x[ic])
                    end               
                end
            end

        end
end     # end update frames