"""
    lower_lim

        function calculating the frame number i1 such that  h*i1 <= t < h*(i1+1) 

            Input: 

                t       time of evaluation
                T       simulation time interval
                Npas    number of time steps

            Outfut:

                i1      frame number such that  0 < h*i1 <= t < h*(i1+1) <= Npas

            Calling sequence:

                i1 = lower_lim(t,T,Npas)


"""
function lower_lim(t::Float64,T::Float64,Npas::Int)
    h = T/(Npas-1)
    i1 = Int(floor(t/h))+1
    return i1
    end
    