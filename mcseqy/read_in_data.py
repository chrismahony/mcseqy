def load_10x_h5_to_anndata(file_paths, sample_names=None, combine=False):
    """
    Load 10x Genomics .h5 files into AnnData objects and optionally combine them.
    
    Parameters:
    - file_paths: List of strings or Path objects, each path pointing to a 10x .h5 file.
    - sample_names: List of sample names, one per file path. Defaults to the file name if not provided.
    - combine: Boolean, if True, combines all matrices into a single AnnData object.
    
    Returns:
    - If combine is False, returns a list of AnnData objects.
    - If combine is True, returns a single combined AnnData object with sample names.
    """
    
    if sample_names is None:
        sample_names = [Path(fp).stem for fp in file_paths]
    
    if len(sample_names) != len(file_paths):
        raise ValueError("The length of sample_names must match the number of file_paths.")

    adata_list = []
    
    for file_path, sample_name in zip(file_paths, sample_names):
        # Read the .h5 file into an AnnData object
        adata = sc.read_10x_h5(file_path)
        
        # Ensure gene names are unique
        adata.var_names_make_unique()
        
        # Add sample name as metadata
        adata.obs['sample'] = sample_name
        
        # Append AnnData object to list
        adata_list.append(adata)
    
    # Combine into a single AnnData object if specified
    if combine:
        adata_combined = sc.concat(adata_list, join='outer', label="sample", keys=sample_names)
        return adata_combined
    
    return adata_list
