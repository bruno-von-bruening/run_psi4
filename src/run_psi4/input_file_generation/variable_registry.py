
var_mem=dict()
def get_var(key):
    assert key in var_mem.keys(), f"The requested \'{key}\' was not set in variable memory" 
    return var_mem[key]