"""
   record_animation

    function to create a video (.mp4) file  froma a mesh record_animation

        Input:

            node_positions::Array{Float64,3}   3 x Nnodes x Npas array collecting node_positions versus time. 
            triangles::Array{Int,2}            3 x Nt array collecting a set of triangles supported by the node mesh.
            T                                  simulation time interval.
            Npas                               number of time steps of the simulation.
            "name"                             name of the video file.
            nframes                            number of video fremes.
            framerate                          numer of frames/s.
            t_init                             simulation time at which the visualization is initialized (first frame must be 3D).

        Output:
            "name.mp4" file.

        Calling sequence: 
            record_animation(node_positions,triangles,T,Npas,name,nframes, framerate)

"""
function record_animation(node_positions::Array{Float64,3},triangles::Array{Int,2},
    T::Float64,Npas::Int,mp4_file::String,nframes::Int, framerate::Int,t_init::Float64)

   
    frames = 1:1:nframes

    hfin = T/(nframes-1)

    Xmin = findmin(node_positions[1,:,:])[1]
    Xmax = findmax(node_positions[1,:,:])[1]
    Ymin = findmin(node_positions[2,:,:])[1]
    Ymax = findmax(node_positions[2,:,:])[1]
    Zmin = findmin(node_positions[3,:,:])[1]
    Zmax = findmax(node_positions[3,:,:])[1]

    lims = ([Xmin,Ymin,Zmin], [Xmax, Ymax, Zmax])

    # lims = HyperRectangle(Vec3f0(Xmin,Ymin,Zmin),Vec3f0(Xmax, Ymax, Zmax))

    # scene = Scene(limits = lims,aspect = (1,1,1))

    scene = Scene()
    
    # update_limits!(scene,lims)


    shape = Observable(shape_interpol(node_positions,t_init,T,Npas))

    cylinder_radius = 0.05
    cylinder_height = 0.1
    n_segments = 20

    cyl_vertices, cyl_faces = cylinder(cylinder_radius, cylinder_height, n_segments)
    # mesh!(scene, shape, triangles, color = (:cyan,0.8), transparency = true)
    # mesh!(scene, cyl_vertices', cyl_faces, color = (:red,0.2),transparency = true)

    f = mesh(shape,triangles,color =:cyan, transparency = true )
    mesh!(scene, cyl_vertices', cyl_faces, color=:red)
    display(f)
    

    record(f,mp4_file*".mp4",frames; framerate ) do i
    t= (i-1)*hfin
    shape[] = shape_interpol(node_positions,t,T,Npas)
    end
end