from . import *

from qcp_global_utils.shell_processes.execution import run_shell_command

@val_call
def execute_psi4_shell(input_file:pdtc_file, psi4_executable:str|None=None, verbose:int=0 ):
    
    import subprocess
    if psi4_executable is None:
        psi4_executable=os.path.join( os.path.dirname(os.path.realpath(sys.executable)),  'psi4' )
        if not os.path.isfile(psi4_executable): raise Exception(f"Cannot find 'psi4' executable in same folder as  {sys.executable}")
    #else:
    #    assert os.path.isfile(conda_exc), f"The conda executable does not exist: {conda_exc}"
    #    psi4_executable=f"{conda_exc} run -n {conda_env} psi4"

    output_file=input_file.split('.')[0]+'.psi4out'
    cmd=f"{psi4_executable} {input_file} {output_file}" # the \" are important for esacping e.g. parenthesis in ccsd(t)
    print(f"Excuting \'{cmd}\' in shell")



    @val_call
    def find_clean_file(pwd=os.getcwd()):
        """ Find files in the current working directory that are likely to be scratch files created by psi4 (e.g. files with .tmp extension) and delete them (these files are not needed and can take up a lot of space) """
        # Find files with .tmp extension
        the_format=r"psi[4]?.[0-9]+.clean"
        found_files = [ x for x in os.listdir(pwd) if re.fullmatch(the_format, x) ]
        if len(found_files)==0:
            raise Exception(f"No clean files found in {pwd}")
        elif len(found_files) > 1:
            raise Exception(f"Found multiple clean files in {pwd} (these files are likely to be scratch files created by psi4, but this is unexpected, because there should only be one clean file per run):\n   {found_files}\n   Will not delete any files to be safe!")
        return found_files[0]


    @val_call
    def get_checkpoint_file(id):
        scratch=os.environ['PSI_SCRATCH']
        escaped_prefix_path = r".*"+re.escape(f".default.{id}.")+ r"[0-9]+\.npy$"
        # try to find this file file match on system
        return [ os.path.join(scratch,x) for x in os.listdir(scratch) if re.fullmatch(escaped_prefix_path, x) ]

    @val_call
    def clean_up(file:pdtc_file):
        """ Excpects files in which every line poiints to absolute path of file
        commonly produced by psi4
        """
        with open(file,'r') as rd:
            lines=rd.readlines()
        # Raise warning if invalid files are found (make a list of valid and invalid files)
        valid_files, invalid_files = [], []
        for x in [ x.strip() for x in lines]:
            if os.path.isfile(x):
                valid_files.append(x)
            else:
                invalid_files.append(x)

        # Also .npy files are apparently dropped that are not included
        try:
            if len(valid_files)>0:
                valid_files+=get_checkpoint_file(valid_files[0].split('.')[1])
        except Exception as ex:
            warn(f"Problem in finding .npy files in scratch folder: {ex}")

        if len(invalid_files) > 0:
            warn(f"Found {len(invalid_files)} invalid files in {file} (these files do not exist):\n   {invalid_files}\n   Will not delete these files!")
        
        # Continously print progress (clean line and print again , like progress bar (no new line, reset cursor to beginning of line) and remove the valid files)
        for i, x in enumerate(valid_files):
            os.remove(x)
            print(f"Removed {x} ({i+1}/{len(valid_files)})", end='\r')

        # move the clean file to back.
        print(f"Moving clean file {file} to back.{os.path.basename(file)}")
        os.rename(file, os.path.join(os.path.dirname(file), 'back.' + os.path.basename(file)))

    def clean():
        print(f"Will manually clean up scratch files", end='\r')
        clean_file=find_clean_file()
        print(f"Found clean file: {clean_file}. Cleaning up scratch files...")
        clean_up(clean_file)
    
    
    def periodic_check(p, data={}):
        if len(data)==0:
            data=dict(
                disk_max_in_GB=0,
                disk_avr_in_GB=0,
                n=1
            )
        try:
            clean_file=find_clean_file()
        except Exception as ex:
            if verbose>0:
                warn(f"Could not find clean file in periodic check: {ex}")
            clean_file=None
        if clean_file:
            with open(clean_file,'r') as rd:
                files=[ x.strip() for x in rd.readlines()]
            if len(files)>0:
                files+=get_checkpoint_file( files[0].split('.')[1] )
            # aggregate the size of these files
            cur_data=0
            for file in files:
                try:
                    cur_data += os.path.getsize(file) / (1024**3)
                except Exception as ex:
                    pass # pretty normal since files may have been deleted after reading them
                    # warn(f"Problem in getting size of file {file}: {ex}")
            if cur_data > data['disk_max_in_GB']:
                data['disk_max_in_GB']=cur_data

            data['disk_avr_in_GB']=data['disk_avr_in_GB']*(data['n']-1)/data['n'] + cur_data/data['n']
        data['n']+=1
        return data 
    try:
        start=time.time()
        di = run_shell_command(cmd, uncritical=True, monitor_resources=True, resources_refresh=2.0,
                               periodic_check=periodic_check, refresh=2)
        took_time=time.time()-start
    except KeyboardInterrupt as ex:
        clean()
        raise KeyboardInterrupt(f"KeyboardInterrupt: Cleaned scratch files") from ex

    return_code, pid, stdout, stderr = di['returncode'], di['pid'], di['stdout'], di['stderr']
    resources=di.get('resources', None)
    scratch_info=di.get( getattr(periodic_check, "__name__" ).split(".")[-1], None)

    if resources is None: 
        warn(f"Did not get resource usage information from the shell command execution (this is unexpected, but not critical, so I will continue without this information)")
        resources={}
    resources['disk_usage']=scratch_info

    debug=False
    def interprete_code(code):
        if code==139:
            return "Insufficient memory caught inside psi4"
        else:
            return f"Unknown error code {code}"
    def format_output(s:List[str]):
        return '\n  |'.join(['']+s)

    if return_code != 0:
        if not debug:
            print(f"Will clean up scratch files before quitting")
            try:
                clean()
            except Exception as ex:
                print(f"Cleaning failed with exception: {ex}")
        quit(f"Psi4 finished with return code {return_code}" +
            ( f"(hint=\'{interprete_code(return_code)}\')" if interprete_code(return_code) else "") +
             f" (after {took_time:.3f} seconds):\n   stdout={format_output(stdout)}  error={format_output(stderr)}")
    else:
        print(f"Psi4 finished successfully (total real time {took_time:.3f} seconds)")

    return {
        'psi4_output_file': output_file,
        'stdout': stdout,
        'stderr': stderr,
        'return_code': return_code,
        'resources': resources,
    }