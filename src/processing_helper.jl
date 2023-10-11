function match_strains(geno::DataFrame, pheno_cropped::DataFrame, genome_start_index::Int)
    pheno_strains = pheno_cropped[:,"id"]
    gInfo = geno[:,1:(genome_start_index-1)]
    geno_strains = names(geno)[genome_start_index:end]
    shared_strains = intersect(pheno_strains, geno_strains)
    genome_df = geno[:,shared_strains]
    pheno_cropped = filter(row -> row["id"] in shared_strains, pheno_cropped)
    return gInfo, genome_df, pheno_cropped
end

function remove_duplicates(pheno::DataFrame)
    duplicate_list = []
    for trait in names(pheno)
        if length(unique(pheno[:,trait])) == 1
            push!(duplicate_list, trait)
        end
    end
    return pheno[:, setdiff(names(pheno), duplicate_list)]
end

function apply_filter_by_value(pheno::DataFrame, filter::Array)
    select = filter[1]
    trait = filter[3]
    val_type = eltype(pheno[:,trait])
    value = parse(val_type, filter[2])
    if select == "Keep"
        filter!(row -> row[trait] == value, pheno)
    elseif select == "Remove"
        filter!(row -> row[trait] != value, pheno)
    end
end

function apply_filter_by_index(pheno::DataFrame, filter::Array)
    select = filter[1]
    index_input = filter[2]
    if contains(index_input, "-")
        split_str = split(index_input, "-")
        start_ind = parse(Int64, split_str[1])
        end_ind = parse(Int64, split_str[2])
        if select == "Remove"
            deleteat!(pheno, start_ind:end_ind)
        else
            deleteat!(pheno_nomissing, (end_ind + 1):(size(pheno_nomissing)[1]))
            deleteat!(pheno, 1:(start_ind-1))
        end
    else
        if select == "Remove"
            deleteat!(pheno, parse(Int64,index_input))
        end
    end
end

function remove_outliers(pheno::DataFrame, trait::String)
    mean = Statistics.mean(pheno[:, trait])
    std = Statistics.std(pheno[:,trait])
    zscores = (pheno[:,trait].-mean)./std
    for i in length(zscores):-1:1
        if zscores[i] > 3 || zscores[i] < -3
            deleteat!(pheno, i)
        end
    end
end

function process_data(geno::DataFrame, pheno::DataFrame, traitID::String; genome_start_index=5, geno_missing_value="U", transformation="None", filters=nothing, remove_outliers="False")
    geno_nomissing = subset(geno, All() .=> ByRow(!=(geno_missing_value)))
    pheno_nomissing = dropmissing(pheno[:,["id",traitID]])
    if filters !== nothing
        for filter in filters
            if (filter[1] == "filter by value")
                apply_filter_by_value(pheno_nomissing, filter)
            elseif (filter[1] == "filter by index")
                apply_filter_by_index(pheno_nomissing, filter)
            end
        end
    end
    if remove_outliers == "True"
        pheno_nomissing = remove_outliers(pheno_nomissing, traitID)
    end
    if transformation == "log2"
        pheno_nomissing[:,traitID] = log2.(pheno_nomissing[:, traitID])
    end
    gInfo, genome_df, pheno_filtered = match_strains(geno_nomissing, pheno_nomissing, genome_start_index)
    geno_processed = Matrix{Float64}(permutedims(genome_df))
    pheno_y = reshape(pheno_filtered[:,traitID],:,1)
    final_strains = pheno_filtered[:,"id"]
    return gInfo, geno_processed, pheno_y, final_strains
end

function process_data_loco(geno::DataFrame, pheno::DataFrame, traitID::String; genome_start_index=5, geno_missing_value="U", transformation="None", filters=nothing, remove_outliers="False")
    geno_nomissing = subset(geno, All() .=> ByRow(!=(geno_missing_value)))
    pheno_nomissing = dropmissing(pheno[:,["id",traitID]])
    if filters !== nothing
        for filter in filters
            if (filter[1] == "filter by value")
                apply_filter_by_value(pheno_nomissing, filter)
            elseif (filter[1] == "filter by index")
                apply_filter_by_index(pheno_nomissing, filter)
            end
        end
    end
    if remove_outliers == "True"
        pheno_nomissing = remove_outliers(pheno_nomissing, traitID)
    end
    if transformation == "log2"
        pheno_nomissing[:,traitID] = log2.(pheno_nomissing[:, traitID])
    end
    gInfo, genome_df, pheno_filtered = match_strains(geno_nomissing, pheno_nomissing, genome_start_index)
    loco_geno = get_loco_geno(hcat(gInfo, genome_df))
    pheno_y = reshape(pheno_filtered[:,traitID],:,1)
    final_strains = pheno_filtered[:,"id"]
    return gInfo, loco_geno, pheno_y, final_strains
end


function process_data(geno::DataFrame, pheno::DataFrame, traitID::String, covar_traitID::String; genome_start_index=5, geno_missing_value="U", transformation="None", filters=nothing)
    geno_nomissing = subset(geno, All() .=> ByRow(!=(geno_missing_value)))
    pheno_nomissing = dropmissing(pheno[:,["id",traitID, covar_traitID]])
    if filters !== nothing
        for filter in filters
            if (filter[1] == "filter by value")
                apply_filter_by_value(pheno_nomissing, filter)
            elseif (filter[1] == "filter by index")
                apply_filter_by_index(pheno_nomissing, filter)
            end
        end
    end

    if remove_outliers == "True"
        pheno_nomissing = remove_outliers(pheno_nomissing, traitID)
    end

    if transformation == "log2"
        pheno_nomissing[:,traitID] = log2.(pheno_nomissing[:, traitID])
    end

    gInfo, genome_df, pheno_filtered = match_strains(geno_nomissing, pheno_nomissing, genome_start_index)
    geno_processed = Matrix{Float64}(permutedims(genome_df))
    pheno_y = reshape(pheno_filtered[:,traitID],:,1)
    covar = reshape(pheno_filtered[:,covar_traitID],:,1)
    final_strains = pheno_filtered[:,"id"]
    
    return gInfo, geno_processed, pheno_y, covar, final_strains
end



function process_data_loco(geno::DataFrame, pheno::DataFrame, traitID::String, covar_traitID::String; genome_start_index=5, geno_missing_value="U", transformation="None", filters=nothing, remove_outliers="False")
    geno_nomissing = subset(geno, All() .=> ByRow(!=(geno_missing_value)))
    pheno_nomissing = dropmissing(pheno[:,["id",traitID, covar_traitID]])
    if filters !== nothing
        for filter in filters
            if (filter[1] == "filter by value")
                apply_filter_by_value(pheno_nomissing, filter)
            elseif (filter[1] == "filter by index")
                apply_filter_by_index(pheno_nomissing, filter)
            end
        end
    end
    if remove_outliers == "True"
        pheno_nomissing = remove_outliers(pheno_nomissing, traitID)
    end
    if transformation == "log2"
        pheno_nomissing[:,traitID] = log2.(pheno_nomissing[:, traitID])
    end
    gInfo, genome_df, pheno_filtered = match_strains(geno_nomissing, pheno_nomissing, genome_start_index)
    loco_geno = get_loco_geno(hcat(gInfo, genome_df))
    pheno_y = reshape(pheno_filtered[:,traitID],:,1)
    covar = reshape(pheno_filtered[:,covar_traitID],:,1)
    final_strains = pheno_filtered[:,"id"]
    return gInfo, loco_geno, pheno_y, covar, final_strains
end



function process_data(geno::DataFrame, pheno::DataFrame; genome_start_index=5, geno_missing_value="U", transformation="None")
    geno_processed = subset(geno, All() .=> ByRow(!=(geno_missing_value)))
    pheno = remove_duplicates(dropmissing!(pheno))
    if transformation == "log2"
        pheno[:, 2:end] = log2.(remove_duplicates(dropmissing!(pheno))[:, 2:end])
    end 
    gInfo, geno_processed, pheno_processed = match_strains(geno_processed, pheno, genome_start_index)
    geno_processed = Matrix{Float64}(permutedims(geno_processed))
    final_strains = pheno_processed[:, "id"]
    pheno_processed = Matrix{Float64}(pheno_processed[:,2:end])
    return gInfo, geno_processed, pheno_processed, final_strains
end



#LOCO Functions taken from BigRiverQTL package to enable loco scans:

function calcLocoKinship(G::Vector{Matrix{Float64}})

	N = length(G)
	K = Vector{Matrix{Float64}}(undef, N)

	Q = calcKinship.(G)
	KK = sum(Q)

	for i in 1:N
		K[i] = (KK - Q[i]) ./ (N - 1)
	end

	return K
end

function get_loco_geno(dfG::DataFrame;
	chromosome_colname::String = "Chr",
	idx_start::Int = 5, # TODO is it better to indicate geno column names?
	kwargs...,
)

	gdf = groupby(dfG, chromosome_colname)

	N = length(gdf)

	G = Vector{Matrix{Float64}}(undef, N)

	for i in 1:N
		G[i] = Matrix{Float64}(permutedims(gdf[i][:, idx_start:end]))
	end

	return G
end

function get_loco_geno_info(dfG::DataFrame;
	chromosome_colname = "Chr",
	idx_info = collect(1:4),
	kwargs...,
)

	select(dfG, idx_info) |>
	x -> groupby(x, chromosome_colname) |>
		 x -> DataFrame.(collect(x))
end


function loco_scan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, K::Vector{Matrix{Float64}};
	kwargs...)

    N = length(G)
	results_loco = map((g, k) -> scan(y, g, k; kwargs...), G, K)

	results_keys = keys(results_loco[1]);

	if :L_perms in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
			L_perms = reduce(vcat, ([results_loco[i].L_perms for i in 1:N])),
		)
	elseif :h2_each_marker in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			h2_each_marker = reduce(vcat, ([results_loco[i].h2_each_marker for i in 1:N])),
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)
	else
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)		
	end
end

function loco_scan(y::Matrix{Float64}, G::Vector{Matrix{Float64}}, covar::Matrix{Float64},
	K::Vector{Matrix{Float64}}; kwargs...)

    N = length(G)
	results_loco = map((g, k) -> scan(y, g, covar, k; kwargs...), G, K)

	results_keys = keys(results_loco[1]);

	if :L_perms in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
			L_perms = reduce(vcat, ([results_loco[i].L_perms for i in 1:N])),
		)
	elseif :h2_each_marker in results_keys
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			h2_each_marker = reduce(vcat, ([results_loco[i].h2_each_marker for i in 1:N])),
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)
	else
		return (sigma2_e = [results_loco[i].sigma2_e for i in 1:N],
			h2_null = [results_loco[i].h2_null for i in 1:N],
			lod = reduce(vcat, ([results_loco[i].lod for i in 1:N])),
		)		
	end
end