from . import *

@val_call
def make_storage_variable(job:psi4_job) -> List[str]:
    # Storage creation
    storage_lines=[
        f"# Store variables here (eg IP)",
        f"from collections import defaultdict # allows for creating key paths beyond one level: dic=defaultdict(dict); dic['a']['b']=c", 
        f"{job.var_names.storage}=defaultdict(dict)",
        f"{job.var_names.storage}.update( files={{}}, results={{}}, settings={{}}, about_the_run={{ 'timings': {{}} }} )  # optional defaultdict will initialise anyway",
        f"{job.var_names.storage_file_name}=\'{job.files.storage}\' # will drop here at the end of the run",
    ]
    return storage_lines

def dump_storage(storage_var, storage_file_var):
    # Storage dump
    storage_drop_lines=[
        f"# Time the complete run",
        f"{storage_var}['about_the_run']['timings'].update({{ 'total_run_time': time.time()-run_started_at,"+
            f" 'unaccounted_time': time.time()-run_started_at - sum(list({storage_var}['about_the_run']['timings'].values()))}})",
        f"# Drop the storage as yaml file",
        f"def drop_yaml(storage, storage_file):",
        f"    with open(storage_file, 'w') as wr:",
        f"        yaml.dump({storage_var},wr)",
        f"    print(f\"Dropped stored variables under {{storage_file}}\")",
            f"drop_yaml( dict({storage_var}), {storage_file_var})",
        ]
    return storage_drop_lines



def update_storage(
        job: psi4_job, 
        kwargs:dict=Field(min_length=1), # Do not pass empty dic!
):
    lines=[]
    store_var=job.var_names.storage

    # loop through the dict
    @val_call
    def loop(
            dic:dict, 
            origin:List[str]=[], # which keys led here
    )-> List[ Tuple[List[str],str] ]:
        """Get the direction and allways return string! """
        res=[]
        for k,v in dic.items():
            path=origin+[k]
            if isinstance(v, dict):
                res+= loop(v, path) 
            elif isinstance(v, str):
                res+=[ (path,v) ]
            else: raise Exception(f"Expected string at branch of dictionary, got: {origin}")
        return res
    
    entries=loop(kwargs)
    print(entries)

    lines+=[
        f"try:",]+[
            ind+store_var+''.join( f"[\'{k}\']" for k in keys) + f" = {v}"
            for keys,v in entries
        ]+[
        f"except Exception as ex:",
        ind+f"print_out(f\"WARNING could not save: {{str(ex)}}\\n\") # probably not necessary as long as defaultdict is used"
    ]
    return lines


@val_call
def store_timing(job:psi4_job, key:str, start:str):
    lines=[
        f"try:",
        ind+f"{job.var_names.storage}['about_the_run']['timings'].update({{ \'{key}\':time.time()-{start} }})",
        f"except Exception as ex:",
        ind+f"print_out(f\"WARNING: could not save timing for \'{key}\'\")",
    ]
    return lines