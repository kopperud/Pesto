function descendant_nodes(node, data)
    desc_edge_idxs = findall(data.edges[:,1] .== node)
    desc = data.edges[desc_edge_idxs,:]
    res = desc[:,2]
end

function parental_node(node, data)
    parental_edge_idx = findall(data.edges[:,2] .== node)
    parent_node = data.edges[parental_edge_idx,1][1]
    return(parent_node)
end

function readtree(treefile)
    @rput treefile
    R"""
    library(ape)

    phy <- read.nexus(treefile)
    nde <- node.depth.edgelength(phy)
    node_depths <- max(nde) - nde
    phy$node_depths <- node_depths
    phy$branching_times <- branching.times(phy)

    po <- postorder(phy)
    phy$po <- po
    """
    @rget phy
    return(phy)
end

function make_SSEdata(phy, datafile, ρ; include_traits = true)
    df = CSV.File(datafile)
   
    if include_traits
        trait_data = Dict(taxon => string(state) for (taxon, state) in zip(df[:Taxon], df[:state]))
    else
        trait_data = Dict(taxon => "?" for taxon in df[:Taxon])
    end
    node_depth = phy[:node_depths]
    tiplab = phy[:tip_label]
    branching_times = phy[:branching_times]

    state_space = sort(unique(values(trait_data)))
#    state_space = unique(string.(convert.(Int64, phy[:tip_state])))
#    tip_state = string.(convert.(Int64, phy[:tip_state]))

 #   trait_data = Dict{String,String}()

  #  if include_traits
  #      for (species, val) in zip(tiplab, tip_state)
  #          if val == "1"
  #              trait_data[species] = state_space[1]
  #          elseif val == "0"
  #              trait_data[species] = state_space[2]
  #          else
  #              trait_data[species] = "?"
  #          end
  #      end
  #  else
  #      trait_data = Dict(sp => "?" for sp in tiplab)
  #  end

    edges = convert.(Int64, phy[:edge])
    el = phy[:edge_length]
    po = phy[:po]

    data = SSEdata(state_space, trait_data, edges, tiplab, node_depth, ρ, el, branching_times, po)
    return(data)
end
