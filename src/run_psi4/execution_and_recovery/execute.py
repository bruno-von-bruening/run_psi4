from . import *

def execute_psi4_shell(input_file, psi4_executable=None ):
    
    import subprocess
    if psi4_executable is None:
        psi4_executable=os.path.join( os.path.dirname(os.path.realpath(sys.executable)),  'psi4' )
        if not os.path.isfile(psi4_executable): raise Exception(f"Cannot find 'psi4' executable in same folder as  {sys.executable}")
    #else:
    #    assert os.path.isfile(conda_exc), f"The conda executable does not exist: {conda_exc}"
    #    psi4_executable=f"{conda_exc} run -n {conda_env} psi4"

    output_file=input_file.split('.')[0]+'.psi4out'
    cmd=f"{psi4_executable} \"{input_file}\" \"{output_file}\"" # the \" are important for esacping e.g. parenthesis in ccsd(t)
    print(f"Excuting \'{cmd}\' in shell")
    start=time.time()
    psi4_run=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = psi4_run.communicate()
    return_code=psi4_run.returncode
    if return_code != 0:
        quit(f"Psi4 finished with return code {return_code}:\n   stdout={stdout.decode()}  error={stderr.decode()}")
    else:
        print(f"Psi4 finished successfully (after {time.time()-start:.3f} seconds)")
    
    return {
        'psi4_output_file': output_file,
        'stdout': stdout.decode(),
        'stderr': stderr.decode(),
        'return_code': return_code,
    }