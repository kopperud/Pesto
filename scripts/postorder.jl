function postorder(model::SSEconstant, data::SSEdata)
    ## Compute the extinction probability through time
    n = length(data.state_space)
    i_not_js = [setdiff(1:n, i) for i in 1:n]
    pE = [model.λ, model.μ, model.η, i_not_js, n]
    alg = RK4()

    tree_height = maximum(data.node_depth)
    tspan = (0.0, tree_height)
    E0 = [1.0 - data.ρ, 1.0 - data.ρ]
    pr = ODEProblem(extinction_prob, E0, tspan, pE)
    E = solve(pr, alg)

    ## Postorder traversal: computing the branch probabilities through time
    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)

    ## Storing the Ds
    D_ends = zeros(typeof(model.λ[1]), nrows, 2)
    ## Storing the scaling factors
    sf = zeros(typeof(model.λ[1]), nrows)

    pD = [model.λ, model.μ, model.η, i_not_js, n, E]

    #for i in 1:nrows
    for i in data.po
        anc, dec = data.edges[i,:]

        #if is tip
        if dec <= Ntip
            species = data.tiplab[dec]
            trait_value = data.trait_data[species]

            if trait_value == "?" ## If we don't know or or didn't observe the trait
                D = [1.0, 1.0] .* data.ρ
            else ## If we observed the trait and measured it 
                trait_idx = convert.(Float64, trait_value .== data.state_space)
                D = trait_idx .* data.ρ
            end

            u0 = typeof(model.λ[1]).(D)

            node_age = data.node_depth[dec]
            parent_node = parental_node(dec, data)
            parent_node_age = data.node_depth[parent_node]
            tspan = (node_age, parent_node_age)

            #prob = remake(pr, u0 = u0, tspan = tspan)
            prob = ODEProblem(backward_prob, u0, tspan, pD)
            sol = solve(prob, alg, save_everystep = false)[end]

            k = sum(sol)
            sol = sol ./ k
            D_ends[i,:] = sol
            logk = log(k)
            sf[i] = logk
        else
            left_node, right_node = descendant_nodes(dec, data)
            left_edge = findall(data.edges[:,2] .== left_node)[1]
            right_edge = findall(data.edges[:,2] .== right_node)[1]
            node_age = data.node_depth[dec]

            D_left = D_ends[left_edge,:]
            D_right = D_ends[right_edge,:]

            D = D_left .* D_right .* model.λ
            u0 = D

            parent_node_age = data.node_depth[anc]
            tspan = (node_age, parent_node_age)

            #prob = remake(pr, u0 = u0, tspan = tspan)
            prob = ODEProblem(backward_prob, u0, tspan, pD)
            sol = solve(prob, alg, save_everystep = false)[end]
            k = sum(sol)
            sol = sol ./ k
            D_ends[i,:] = sol
            if k < 0.0
                logk = 0.0
            else
                logk = log(k)
            end
            sf[i] += logk
        end
    end
    return(D_ends, sf, E)
end
