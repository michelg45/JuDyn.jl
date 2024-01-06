function correct_frames(dy::Vector{Float64},theta_p::Float64,Nnodes::Int)


      cf = Main.Frames.current_frames

      for inode=1:Nnodes

          inv_loc=copy(cf[inode].inv_loc)

          inv_loc_x = inv_loc[1:3]

          if findmin(inv_loc_x)[1] > 0
              cf[inode].Dx =cf[inode].Dx + Vec3(dy[inv_loc_x])
              cf[inode].xdot=cf[inode].xdot + Vec3(theta_p*dy[inv_loc_x])
          else
              for ic=1:3
                      inv_loc_x[ic] > 0 &&
                      (
                        cf[inode].Dx[ic]   = cf[inode].Dx[ic] + dy[inv_loc_x[ic]];
                        cf[inode].xdot[ic] = cf[inode].xdot[ic] + theta_p*dy[inv_loc_x[ic]]
                      )
                              end
          end

#          cf[inode].x = cf[inode].xn + cf[inode].Dx

          if cf[inode].type == "frame"
              inv_loc_r = inv_loc[4:6]
              if findmin(inv_loc_r)[1] > 0
              Dpsi = cf[inode].Dpsi
              dpsi = RV3(dy[inv_loc_r])
              cf[inode].Dpsi   = RV3(Dpsi,dpsi)
              cf[inode].psidot = cf[inode].psidot + Vec3(theta_p*dy[inv_loc_r])
              else
                  Dpsi = RV3()
                  dpsi = RV3()
                  for ic=1:3
                      inv_loc_r[ic] > 0 &&
                       (
                          Dpsi[ic] = cf[inode].Dpsi[ic];
                          dpsi[ic] = dy[inv_loc_r[ic]];
                          cf[inode].psidot[ic] = cf[inode].psidot[ic] + theta_p*dy[inv_loc_r[ic]]
                        )
                  end
                  cf[inode].Dpsi=RV3(Dpsi,dpsi)
              end


#               cf[inode].psi = RV3(cf[inode].psin,  cf[inode].Dpsi)
              cf[inode].T = tang(cf[inode].Dpsi)
              cf[inode].W = cf[inode].T*cf[inode].psidot
              cf[inode].DW = Dtang(cf[inode].psidot,cf[inode].Dpsi)
          end
      end

end # correct_frames
