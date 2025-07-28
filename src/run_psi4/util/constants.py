
BOHR_TO_ANGSTROM=.52917721067
PERIODIC_TABLE_STR = """
H                                                                                                                           He
Li  Be                                                                                                  B   C   N   O   F   Ne
Na  Mg                                                                                                  Al  Si  P   S   Cl  Ar
K   Ca  Sc                                                          Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr
Rb  Sr  Y                                                           Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te  I   Xe
Cs  Ba  La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu  Hf  Ta  W   Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn
Fr  Ra  Ac  Th  Pa  U   Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No  Lr  Rf  Db  Sg  Bh  Hs  Mt  Ds  Rg  Cn  Nh  Fl  Mc  Lv  Ts  Og
"""
atom_types=PERIODIC_TABLE_STR.split()
Z_to_atty=dict([ (i+1, at_ty) for i, at_ty in enumerate(atom_types)])
def atomic_charge_to_atom_type(Z):
    """ Return atom type when given an atomic Number """
    assert Z in Z_to_atty.keys()
    return Z_to_atty[Z]
def atom_type_to_Z(at_ty):
    """ Return integer atomic number for given atom type"""
    answer=[ Z for Z,at in Z_to_atty.items() if at.upper()==at_ty.upper()]
    assert len(answer)==1, answer
    answer=answer[0]
    return answer