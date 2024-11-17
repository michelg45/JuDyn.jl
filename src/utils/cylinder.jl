function  cylinder(R::Float64, h::Float64, N::Int,theta_min::Float64,theta_max::Float64)
    dtheta = theta_max-theta_min
    theta =  [i for i in 0:N]*dtheta/N .+ theta_min
    x = broadcast(x -> R*cos(x),theta)
    y = broadcast(x -> R*sin(x),theta)
    z = [i for i in 0:N]*h/N
    Npoints = (N+1)^2
    vertices = zeros(3,Npoints)
    ipoint = 0
    for i = 1:N+1
        for k = 1:N+1
            ipoint += 1
            vertices[1,ipoint] = x[i]
            vertices[2,ipoint] = y[i]
            vertices[3,ipoint] = z[k]
        end
    end
    Nfaces = 2*N^2
    faces = zeros(Int, 3, Nfaces)
    iface = 0
    for i = 1:N
        for k = 1:N
            i1 = (i-1)*(N+1)+k
            i2 = i*(N+1) +k
            i3 = i*(N+1) +k+1
            i4 = (i-1)*(N+1) +k+1
            iface += 1
            faces[:,iface] = [i1, i2, i3]
            iface += 1
            faces[:,iface] = [i1, i3, i4]
        end
    end
    faces = faces'

    return vertices, faces

    end

function  cylinder(R::Float64, h::Float64, N::Int)
    theta = [i for i in 0:N]*2.0*pi/N
    x = broadcast(x -> R*cos(x),theta)
    y = broadcast(x -> R*sin(x),theta)
    z = [i for i in 0:N]*h/N
    Npoints = (N+1)^2
    vertices = zeros(3,Npoints)
    ipoint = 0
    for i = 1:N+1
        for k = 1:N+1
            ipoint += 1
            vertices[1,ipoint] = x[i]
            vertices[2,ipoint] = y[i]
            vertices[3,ipoint] = z[k]
        end
    end
    Nfaces = 2*N^2
    faces = zeros(Int, 3, Nfaces)
    iface = 0
    for i = 1:N
        for k = 1:N
            i1 = (i-1)*(N+1)+k
            i2 = i*(N+1) +k
            i3 = i*(N+1) +k+1
            i4 = (i-1)*(N+1) +k+1
            iface += 1
            faces[:,iface] = [i1, i2, i3]
            iface += 1
            faces[:,iface] = [i1, i3, i4]
        end
    end
    faces = faces'
    
    return vertices, faces
    
    end