"""
    print_node
        function to print the initial data at node "inode"
        as collected the "Main.node_container" array
        calling sequence: print_node(inode::Int)

        Example

        print_node(2)
        node number: 2
        type: frame
        parent node: 0
        localization: [7, 8, 9, 10, 11, 12]
        inverse localization: [1, 2, 3, 4, 5, 0]
        initial position: Vec3([1.000000e+00, 0.000000e+00, 0.000000e+00])
        initial orientation: RV3([0.000000e+00, 0.000000e+00, 0.000000e+00])

"""
function print_node(inode::Int)
    nc = Main.node_container
    println("node number: ",nc.node_numbers[inode])
    println("type: ", nc.types[inode])
    println("parent node: ",nc.parent_nodes[inode])
    println("localization: ", nc.locs[inode])
    println("inverse localization: ",nc.inv_locs[inode])
    println("initial position: ",nc.init_positions[inode])
    println("initial orientation: ", nc.init_orientations[inode])
end
