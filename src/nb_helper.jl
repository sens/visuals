function get_specieslist()
	specieslist = []
	global dfspecies = list_species()
	for row in eachrow(dfspecies)
		name = row[3]
		push!(specieslist, name)
	end
	return specieslist
end


function get_grouplist(species)
	displayname_grouplist = []
	global dfgroups = list_groups(species)
	for row in eachrow(dfgroups)
		displayname = row[1]
		push!(displayname_grouplist, displayname)
	end
	return displayname_grouplist
end


function get_group_shortname(group_displayname)
	for row in eachrow(dfgroups)
		if (cmp(row[1], group_displayname) == 0)
			return row[6]
		end
	end
end

function get_datasetlist(group_displayname)
	datasetlist = []
	group_shortname = get_group_shortname(group_displayname)
	global dfdatasets = list_datasets(group_shortname)
	for row in eachrow(dfdatasets)
		fullname = row[4]
		push!(datasetlist, fullname)
	end
	return datasetlist
end

function get_dataset_shortabbrev(dataset_fullname)
	for row in eachrow(dfdatasets)
		if (cmp(row[4], dataset_fullname) == 0)
			return row[9]
		end
	end
end


function remove_missing(geno)
	filtered = subset(geno, All() .=> ByRow(!=("U")))
	return filtered
end

function transpose_geno(geno)
	genoT = permutedims(geno)
	return genoT
end

function crop_geno(geno)
	if ("cM" in names(geno)) && ("Mb" in names(geno))
		return geno[:,5:end]
	elseif ("cM" in names(geno)) || ("Mb" in names(geno))
		return geno[:,4:end]
	else
		println("unexpected column names")
	end
end

function download_data(species_name, group_shortname, dataset_shortabbrev, data_dir_path)
	data_raw_geno_filepath = joinpath(data_dir_path, "dataRaw/genoData/genoRawFile")
	data_raw_pheno_filepath = joinpath(data_dir_path, "dataRaw/phenoData/phenoRawFile")
	download_geno(group_shortname, path=data_raw_geno_filepath)
	if dataset == "Clinical Phenotypic Data"
		download_pheno(group_shortname, path=data_raw_pheno_filepath)
	else
		download_omics(dataset_shortabbrev, path=data_raw_pheno_filepath)
	end
end

function create_dictionary(control_filepath)
	new_dict = CSV.File(control_filepath) |> Dict
	return new_dict
end

function filter_geno(control_dictionary, geno)
	new_geno = filter(geno, control_dictionary)
	return new_geno
end

function count_missing_in_pheno(phenodf)
	missing_cnts = DataFrame(id = String[], Missing_Count = Int[])
		
	for col in names(phenodf)[1:end]
		cnt = 0
		for i in 1:(size(phenodf)[1])
			if typeof(phenodf[i,col]) == Missing
				cnt += 1
			end
		end
		push!(missing_cnts, [col, cnt])
	end
	return missing_cnts
end

function get_bxdtrait_description(trait_str, trait_description_df)
	index = findfirst(isequal(trait_str),trait_description_df.trait)
	description = trait_description_df[index, "description"]
	return description
end

function get_chr_indices(gInfo)
	chr_indices = []
	chromosomes = gInfo[:,1]
	first = 1
	prev_chr = chromosomes[1]
	for i in 1:(size(chromosomes)[1])
		curr_chr = chromosomes[i]
		if curr_chr != prev_chr
			range = [first, i-1]
			push!(chr_indices, range)
			first = i
			prev_chr = curr_chr
		elseif i == (size(chromosomes)[1])
			range = [first,i]
			push!(chr_indices,range)
		end
	end
	return chr_indices
end

function get_loco_kinships(geno_processed, chr_indices)
	kinships_array = []
	for i in 1:(size(chr_indices)[1])
		curr = chr_indices[i]
		if i == 1
			loco_geno = geno_processed[:,(curr[2]+1):end]
		elseif i == (size(chr_indices)[1])
			loco_geno = geno_processed[:,1:(curr[1]-1)]
		else
			before = geno_processed[:,1:(curr[1]-1)]
			after = geno_processed[:,(curr[2]+1):end]
			loco_geno = hcat(before, after)
		end
		kinship = calcKinship(loco_geno)
		push!(kinships_array, kinship)
	end
	return kinships_array
end

function get_significant_loci(thresh::Float64, lod_scores, gInfo)
    idx_thresh = findall(lod_scores .>= thresh)
	df_results = gInfo[idx_thresh,:]
	df_results.LOD = lod_scores[idx_thresh]
    return df_results
end

function download_significant_loci(filepath::String, thresh::Float64, lod_scores, gInfo)
    idx_thresh = findall(lod_scores .>= thresh)
	df_results = gInfo[idx_thresh,:]
	df_results.LOD = lod_scores[idx_thresh]
    CSV.write(filepath, df_results)
end
