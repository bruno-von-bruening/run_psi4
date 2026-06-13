from run_psi4.util.import_helper import *

class pol_ff_options(controller_model):
    #cc_convergence: float = 1.e-12
    do_moments: bool = True
    delta_dip: float = 1.e-3
    delta_qad: float = 1.e-4
    # e_convergence: float|None=None

    @classmethod
    def moments_possible(cls, scf_type, method):
        """ MPX and ccX(Y) have perturbative elements and hence do not allow property evaluation (lack of coefficients)! """
        return not ( 
            scf_type.lower()=='conv' and any( re.match( x,method.lower()) for x in [r'cc.*\(\d+\)',r'mp\d'])
        )

    def __init__(self, **data):
        super().__init__(**data)
        # if self.e_convergence is None:
        #     self.e_convergence=self.delta_qad**2/100 # Should be 0.01 exact
class pol_lr_options(controller_model):
    def __init__(self, **data):
        super().__init__(**data)
        # if self.e_convergence is None:
        #     self.e_convergence=self.delta_qad**2/100 # Should be 0.01 exact

class optim_options(controller_model):
    maxiter: int = 50 # the psi4 default https://psicode.org/psi4manual/master/optking.html
    check_hessian: bool = True
    # maxiter: int = 100
    # g_convergence: float = 1.e-5
    # e_convergence: float = 1.e-6

class sp_options(controller_model):
    pass