

"""
intersect_strain(data::PhenoGenoData; pheno_strain_col="StrainName")

Create a new `PhenoGenoData` object by combining phenotype and genotype dataframes
with shared strain names.

# Arguments
- `data`: the original dataset containing phenotype and genotype information.
- `pheno_strain_col`: the name of the column in the phenotype dataframe that 
  holds the strain names.

"""
function intersect_strain(gn_df::DataFrame, phn_df::DataFrame; pheno_strain_col = "StrainName")
    colnames_geno = names(gn_df)
    
    geno_strains = colnames_geno[5:end]
    pheno_strains = phn_df[:, pheno_strain_col]

    shared_strains = intersect(pheno_strains, geno_strains)

    return (
        geno = select(gn_df, vcat(colnames_geno[1:4], shared_strains)),
		pheno = filter(row ->row[pheno_strain_col] in shared_strains, phn_df)
        )
end



