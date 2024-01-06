
"""
    create_model()

        function initiating the JUDYN model. 
            It creates the 3 general containers of the JUDYN model: 'JuDyn.FORDYN.model_container', 
            'JuDyn.FORDYN.element_container' and 'JuDyn.FORDYN.node_container'.

            calling sequence: create_model()
            
"""
function create_model()
    global node_container = NodeArray()
    global model_container = ModelArray()
    global element_container= ElementArray()
    println("model_container created")
    println("node_container created")
    println("element_container created")
    return model_container, node_container, element_container
end