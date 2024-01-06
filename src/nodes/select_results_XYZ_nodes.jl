"""
    select_results_XYZ_nodes
    
        function returning the locations (idx,idy,idz)
        in the structural set of the translation components
        of the set "nodes" of frame nodes. 
        calling sequence: idx,idy,idz = select_results_XYZ_nodes(nodes::Vector{Int})

        Example

        nodes = [2,3]
        2-element Vector{Int64}:
            2
            3
        idx, idy, idz  = select_results_XYZ_nodes(nodes)
        ([1, 6], [2, 7], [3, 8])

"""
function select_results_XYZ_nodes(nodes::Vector{Int})

Nnodes = max(size(nodes,1),size(nodes,2))
idx = Vector{Int}(zeros(Nnodes))
idy = Vector{Int}(zeros(Nnodes))
idz = Vector{Int}(zeros(Nnodes))

 for i =1:Nnodes
       idx[i]=find_node_components(nodes[i],1)
       idy[i]=find_node_components(nodes[i],2)
       idz[i]=find_node_components(nodes[i],3)
 end


 return idx,idy,idz

end
