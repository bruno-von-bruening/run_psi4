from . import *

@val_call
def set_hardware(memory='4GB', num_threads=4, storage_var=None):
    hardware_settings=dict(
        memory=p4ut.embed_value(memory),
        num_threads=p4ut.embed_value(num_threads),
    )
    content=[]
    content+=[f"# Hardware settings "]
    hardware_lines=[ 4*" "+f"{k}={v}," for k, v in hardware_settings.items() ]
    content+=["hardware_settings=dict("]+ hardware_lines +[ ")"]
    content+=[f"set_{k}(hardware_settings['{k}'])" for k ,v in hardware_settings.items()]
    
    if not isinstance(storage_var, type(None)):
        content+=[ f"{storage_var}['about_the_run'].update(hardware_settings)"]
    content='\n'.join([ f"{line}" for line in content])
    return content

def machine_settings_wrapper(v):
    if isinstance(v, machine_settings):
        return v
    elif isinstance(v, dict):
        return machine_settings(**v)
    else:
        raise Exception(f"Cannot wrap machine_settings from {v} of type {type(v)}")
    
class machine_settings(BaseModel):
    num_threads     : int = 4    
    memory          : str = '4GB'
    @val_call
    def hardware_lines(self, storage_var:str):
        return set_hardware(memory=self.memory, num_threads=self.num_threads, storage_var=storage_var).split('\n')
    def __init__(self, **kwargs):
        if 'memory' in kwargs.keys():
            mem=kwargs['memory']
            try:
                float(mem)
                mem=f"{mem.strip()}GB"
            except:
                assert any([ mem.upper().endswith(x) for x in ['MB','GB'] ]), f"wrong memory {mem}"
                number=mem.upper().removesuffix('GB').removesuffix('MB').rstrip('_')
                try:
                    float(number)
                except:
                    raise Exception(f"wrong memory {mem}. Not a number: {number}")
                mem=f"{number} {mem[-2:].upper()}"
            kwargs['memory']=mem
        super().__init__(**kwargs)


my_machine_settings=Annotated[machine_settings, BeforeValidator(machine_settings_wrapper)]