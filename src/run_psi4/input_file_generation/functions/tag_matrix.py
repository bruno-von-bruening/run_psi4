def tensor_tag_from_indices(indices):
    """Return the tensor tag string for a given index combination (e.g. [0,1,2] -> 'xyz')."""
    mapper = ["x", "y", "z"]
    return ''.join(mapper[i] for i in indices)

def tag_matrix(M, field=None):
    def get_val(ind):
        return float(M[tuple(ind)])
    new_M={}
    M = np.asarray(M)
    base_tag=f"@{field}" if not field is None else ''

    dims=M.shape
    # get all indices of M
    all_indices = list(np.ndindex(M.shape))
    for ind in all_indices:
        tag=tensor_tag_from_indices(ind)
        new_M.update({tag+base_tag: get_val(ind)})

    return new_M