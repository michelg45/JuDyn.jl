function shape_functions_2D(Nnodes,xi)

    if Nnodes == 3

          F =  [1.0-xi[1]-xi[2]; xi[1]; xi[2]]

          DF = [-1.0 -1.0; 1.0  0.0; 0.0  1.0]

    elseif Nnodes == 4

          F =  [0.25*(1.0-xi[1])*(1.0-xi[2]);
               0.25*(1.0+xi[1])*(1.0-xi[2]);
               0.25*(1.0+xi[1])*(1.0+xi[2]);
               0.25*(1.0-xi[1])*(1.0+xi[2])]

          DF = [-0.25*(1.0-xi[2])   -0.25*(1.0-xi[1]);
               0.25*(1.0-xi[2])   -0.25*(1.0+xi[1]);
               0.25*(1.0+xi[2])    0.25*(1.0+xi[1]);
              -0.25*(1.0+xi[2])    0.25*(1.0-xi[1])]

    elseif Nnodes == 5

        f_E = (1.0-xi[1])*(1.0-xi[1])*(1.0-xi[2])*(1.0-xi[2])

        F = [0.25*(1.0-xi[1])*(1.0-xi[2])-0.25*f_E;
        0.25*(1.0+xi[1])*(1.0-xi[2])-0.25*f_E;
        0.25*(1.0+xi[1])*(1.0+xi[2])-0.25*f_E;
        0.25*(1.0-xi[1])*(1.0+xi[2])-0.25*f_E;
        f_E]

        d1_f_E =-2.0*xi[1]*(1.0-xi[2])*(1.0-xi[2])
        d2_f_E =-2.0*xi[2]*(1.0-xi[1])*(1.0-xi[1])

        DF = [-0.25*(1.0-xi[2])-0.25*d1_f_E  (-0.25)*(1.0-xi[1])-0.25*d2_f_E;
        0.25*(1.0-xi[2])-0.25*d1_f_E    (-0.25)*(1.0+xi[1])-0.25*d2_f_E;
        0.25*(1.0+xi[2])-0.25*d1_f_E   0.25*(1.0+xi[1])-0.25*d2_f_E;
        -0.25*(1.0+xi[2])-0.25*d1_f_E  0.25*(1.0-xi[1])-0.25*d2_f_E;
        d1_f_E  d2_f_E]

    elseif Nnodes == 6
        
        F = [1.0 - 3.0*xi[1] - 3.0*xi[2] + 2.0*xi[1]^2 +  4.0*xi[1]*xi[2] + 2.0*xi[2]^2; 
             -1.0*xi[1]  + 2.0*xi[1]^2; 
             -1.0*xi[2]  + 2.0*xi[2]^2 ;
             4.0*xi[1] - 4.0*xi[1]^2 - 4.0*xi[1]*xi[2]; 
             4.0*xi[1]*xi[2];
             4.0*xi[2] - 4.0*xi[2]^2 - 4.0*xi[1]*xi[2]]

        DF=[-3.0+4.0*xi[1]+4.0*xi[2] -3.0+4.0*xi[1]+4.0*xi[2];
             -1.0+4.0*xi[1] 0.0;
             0.0 -1.0+4.0*xi[2];
             4.0-8.0*xi[1]-4.0*xi[2] -4.0*xi[1];
             4.0*xi[2] 4.0*xi[1];
             -4.0*xi[2] 4.0-8.0*xi[2]-4.0*xi[1]]

    elseif Nnodes == 8

         eta = xi[1]
         mu = xi[2]

          F = [-0.25*(1.0-eta)*(1.0-mu)*(1.0+eta+mu);
              -0.25*(1.0+eta)*(1.0-mu)*(1.0-eta+mu);
              -0.25*(1.0+eta)*(1.0+mu)*(1.0-eta-mu);
              -0.25*(1.0-eta)*(1.0+mu)*(1.0+eta-mu);
              0.5*(1.0+eta)*(1.0-eta)*(1.0-mu);
              0.5*(1.0+eta)*(1.0+mu)*(1.0-mu);
              0.5*(1.0+eta)*(1.0-eta)*(1.0+mu);
              0.5*(1.0-eta)*(1.0+mu)*(1.0-mu)]

          DF = [-0.25*(1.0-mu)*(-mu-2.0*eta)  -0.25*(1.0-eta)*(-eta-2.0*mu);
              -0.25*(1.0-mu)*(mu-2.0*eta)   -0.25*(1.0+eta)*(eta-2.0*mu);
              -0.25*(1.0+mu)*(-mu-2.0*eta)  -0.25*(1.0+eta)*(-eta-2.0*mu);
              -0.25*(1.0+mu)*(mu-2.0*eta)   -0.25*(1.0-eta)*(eta-2.0*mu);
              -(1.0-mu)*eta                 -0.5*(1.0-eta)*(1.0+eta);
               0.5*(1.0-mu)*(1.0+mu)        -(1.0+eta)*mu;
              -(1.0+mu)*eta                  0.5*(1.0-eta)*(1.0+eta);
              -0.5*(1.0-mu)*(1.0+mu)        -(1.0-eta)*mu]
    elseif Nnodes == 9
                eta = xi[1]
                mu = xi[2]
               
                O_eta = 1.0 - eta
                O_mu = 1.0 - mu
                Opeta = 1.0 + eta
                Opmu = 1.0 + mu
               
                etamu = eta*mu
                O_eta2 = 1.0 - eta*eta
                O_mu2 = 1.0 - mu*mu

                O_2eta = 1.0 - 2.0*eta
                O_2mu = 1.0 - 2.0*mu
                Op2eta = 1.0 + 2.0*eta
                Op2mu = 1.0 + 2.0*mu              
               
                F =   [0.25*O_eta*O_mu*etamu;
                       -0.25*Opeta*O_mu*etamu;
                       0.25*Opeta*Opmu*etamu;
                       -0.25*O_eta*Opmu*etamu;
                       -0.5*mu*O_mu*O_eta2;
                       0.5*eta*Opeta*O_mu2;
                       0.5*mu*Opmu*O_eta2;
                       -0.5*eta*O_eta*O_mu2;
                       O_eta2*O_mu2]

                DF =  [0.25*mu*O_mu*O_2eta  0.25*O_eta*eta*O_2mu;
                      -0.25*mu*O_mu*Op2eta  -0.25*Opeta*eta*O_2mu;
                      0.25*mu*Opmu*Op2eta  0.25*Opeta*eta*Op2mu;
                      -0.25*mu*Opmu*O_2eta  -0.25*O_eta*eta*Op2mu;
                      mu*O_mu*eta  -0.5*O_eta2*O_2mu;
                      0.5*O_mu2*Op2eta  -eta*Opeta*mu;
                      -mu*Opmu*eta  0.5*O_eta2*Op2mu;
                     -0.5*O_mu2*O_2eta  eta*O_eta*mu;
                      -2.0*eta*O_mu2  -2.0*mu*O_eta2]

    end
    return F, DF
 end