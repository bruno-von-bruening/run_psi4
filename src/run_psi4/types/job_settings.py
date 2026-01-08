from run_psi4.util.import_helper import *

class pol_ff_options(controller_model):
    #cc_convergence: float = 1.e-12
    do_moments: bool = True
    delta_dip: float = 1.e-3
    delta_qad: float = 1.e-4
    e_convergence: float|None=None

    @classmethod
    def moments_possible(cls, scf_type, method):
        return ( 
            scf_type.lower()=='df' and any(method.lower().startswith(x) for x in ['cc','mp'])
        )

    def __init__(self, **data):
        super().__init__(**data)
        if self.e_convergence is None:
            self.e_convergence=self.delta_qad**2/100 # Should be 0.01 exact