from . import *

@val_call
def make_storage_variable(storage_var, storage_file_var, storage_file) -> Tuple[List[str], str]:
    # Storage creation
    storage_lines=[
        f"# Store variables here (eg IP)",
        f"{storage_var}={{'files':{{}},  'results':{{}}, 'settings':{{}}, 'about_the_run':{{ 'timings': {{}} }} }}",
        f"{storage_file_var}=\'{storage_file}\'",
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
            f"drop_yaml({storage_var}, {storage_file_var})",
        ]
    return storage_drop_lines