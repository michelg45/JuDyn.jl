"""
    select_results_nodes

        function returning the vectors "ipos", "irot"
        containing the locations in the structural set
        of the translation and rotation components
        of the set "nodes" of frame nodes. 
        calling sequence: ipos, irot = select_results_nodes(nodes::Vector{Int})

        Example

        nodes = [2,3]
        2-element Vector{Int64}:
            2
            3
        ipos, irot = select_results_nodes(nodes)
        ([[1, 2, 3], [6, 7, 8]], [[4, 5, 0], [9, 10, 0]])    

"""
function select_results_nodes(nodes::Vector{Int})

Nnodes = max(size(nodes,1),size(nodes,2))
ipos = Vector{Vector{Int}}(undef,Nnodes)
irot  = Vector{Vector{Int}}(undef,Nnodes)

 for i =1:Nnodes
       ipos[i]=[find_node_components(nodes[i],1); find_node_components(nodes[i],2); find_node_components(nodes[i],3)]
       irot[i]=[find_node_components(nodes[i],4); find_node_components(nodes[i],5); find_node_components(nodes[i],6)]
 end

 return ipos, irot

end