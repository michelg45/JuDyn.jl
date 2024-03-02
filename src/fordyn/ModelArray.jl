"""
    ModelArray

Data structure containing the general data of the model. Its full content can be consulted in the source code.

"""
mutable struct ModelArray

    Name::String
    Nodes::Int
    Elements::Int
    Nodal_forces::Int
    Nodal_imposed_displacements::Int
    Nodal_torques::Int
    SuperElements::Int
    Rigid_bodies::Int
    Rigid_masses::Int
    Node_links::Int
    Frame_links::Int
    Frame_springs::Int
    Ground_hinges::Int
    Hinges::Int
    Prismatic_joints::Int
    Spherical_joints::Int
    Beams::Int
    Shells::Int
    SuperBeams::Int
    Ground_spring_dampers::Int
    Inequalities::Int
    Linear_constraints::Int
    Ndofs::Int
    Ndofs_x::Int
    Ndofs_int::Int
    Ndofs_q::Int
    Ndofs_v::Int
    Ndofs_mult::Int
    max_x::Int
    max_int::Int
    max_v::Int
    max_mult::Int
    end_of_nodes::Bool
    end_of_elements::Bool
    end_of_assembly::Bool
    end_of_initial_conditions::Bool
    struc_loc_q::Vector{Int}
    struc_loc_v::Vector{Int}
    modified_locs::Vector{Vector{Int}}
    initial_shape::Bool
    gravity::Vec3
    uniform_rotation::Bool
    rotation_speed::Vec3
    matrix_update::Bool
    eigvals::Bool
    init_file::String


"""
    ModelArray()

Function creating an array of ModelArray type.

Calling sequence (once by the create_model function) : 
      
            global model_container = ModelArray() 
"""
function ModelArray(Name)

        Nodes=0
        Elements=0
        Nodal_forces=0
        Nodal_torques=0
        Nodal_imposed_displacements=0
        SuperElements=0
        Rigid_bodies=0
        Rigid_masses=0
        Node_links=0
        Frame_links=0
        Frame_springs=0
        Ground_hinges = 0
        Hinges = 0
        Prismatic_joints = 0
        Spherical_joints = 0
        Beams=0
        Shells=0
        SuperBeams=0
        Ground_spring_dampers = 0
        Inequalities=0
        Linear_constraints = 0
        Ndofs=0
        Ndofs_x=0
        Ndofs_int=0
        Ndofs_q=0
        Ndofs_v=0
        Ndofs_mult=0
        max_x=0
        max_int=0
        max_mult=0
        max_v=0
        end_of_nodes=false
        end_of_elements=false
        end_of_assembly=false
        end_of_initial_conditions=false
        struc_loc_q=Vector{Int}[]
        struc_loc_v=Vector{Int}[]
        modified_locs = []
        initial_shape = false
        gravity = Vec3(0.0, 0.0, 0.0)
        uniform_rotation = false
        rotation_speed = Vec3(0.0, 0.0, 0.0)
        matrix_update = true
        eigvals = false
        init_file = ""
        
        
        return new(Name,Nodes,Elements,Nodal_forces,Nodal_torques,Nodal_imposed_displacements,
                   SuperElements,Rigid_bodies,Rigid_masses,
                   Node_links,Frame_links,Frame_springs,Ground_hinges,
                   Hinges,Prismatic_joints,Spherical_joints,Beams,Shells,SuperBeams,Ground_spring_dampers,Inequalities,Linear_constraints,Ndofs,Ndofs_x,Ndofs_int,Ndofs_q,
                   Ndofs_v,Ndofs_mult,max_x,max_int,max_v,max_mult,
                   end_of_nodes,end_of_elements,end_of_assembly,
                   end_of_initial_conditions,struc_loc_q,struc_loc_v,modified_locs,initial_shape,gravity,
                   uniform_rotation,rotation_speed,matrix_update,eigvals,init_file)

    end

end

"""
ModelArray()

function creating the "Main.model_container".

Calling sequence (once by the create_model function)

    global model_container = ModelArray(name::String)

    or 

    global model_container = ModelArray()
    in which case, name = "new_problem"

"""
function ModelArray()

    name="new_problem"
    return ModelArray(name)

end