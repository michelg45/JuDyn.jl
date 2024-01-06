function update_frames(y::Vector)

    """
        function update_frames

        updates node reference configurations at end of step
        in current_frames:  x := x + Dx, psi := RV3(psi,Dpsi)


    """

        Nnodes=Main.model_container.Nodes
        cf = Main.Frames.current_frames



        for inode=1:Nnodes

            inv_loc=copy(cf[inode].inv_loc)
            inv_loc_x = inv_loc[1:3]
            inv_loc_r = inv_loc[4:6]

            # update position

            cf[inode].xn = cf[inode].x
            cf[inode].Dx=Vec3(copy(zeros(3)))
            cf[inode].xdot = Vec3(zeros(3))

            if findmin(inv_loc_x)[1] > 0
                y[inv_loc_x] = cf[inode].x.v

            else
                for ic=1:3
                    inv_loc_x[ic] > 0 &&  (y[inv_loc_x[ic]] = cf[inode].x[ic])
                end
            end

            # update rotation for "frame" nodes

            if cf[inode].type == "frame"

                cf[inode].psin =  cf[inode].psi
                cf[inode].Dpsi=RV3(copy(zeros(3)))
                cf[inode].psidot=Vec3(copy(zeros(3)))
                cf[inode].T=Mat3(copy(eye(3)))
                cf[inode].W=Vec3(copy(zeros(3)))
                cf[inode].DW=Mat3(copy(zeros(3,3)))

                if findmin(inv_loc_r)[1] > 0
                    y[inv_loc_r] = cf[inode].psi.v

                else
                    for ic=1:3
                        inv_loc_r[ic] > 0 && ( y[inv_loc_r[ic]] = cf[inode].psi[ic] )
                    end
                end

            end

        end

end # update_frames
