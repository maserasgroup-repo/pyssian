"""
Contains all the hardcoded basis and functionals as well as the neccesary
regex expressions to identify if a certain string is a valid basis or a valid
methods. It also contains a dictionary that maps Atomic Symbol with Atomic
Number and vice versa.
"""

import re

######################## METHODS SECTION #####################################
shellWaveFunctions = 'RO|U|R|'
mutlipleMethodRegex = 'ONIOM|IRCMAX'

methodsRegex = ['DFT','MM','Amber','Dreiding','UFF','AM1','PM3','PM3MM','PM6',
                'PDDG','HF','HFS','XAlpha','HFB','VSXC','HCTH','HCTH93',
                'HCTH147','HCTH407','tHCTH','M06L','B97D','LSDA','LC-wPBE',
                'CAM-B3LYP','wB97XD','wB97','wB97X','LC-BLYP','B3LYP','B3P86',
                'B3PW91','B1B95','mPW1PW91','mPW1LYP','mPW1PBE','mPW3PBE',
                'B98','B971','B972','PBE1PBE','B1LYP','O3LYP','TPSSh','BMK',
                'M062X','M06','M06HF','M05','M052X','X3LYP','BHandH',
                'BHandHLYP','tHCTHhyb','HSEh1PBE','HSE2PBE','HSEhPBE',
                'PBEh1PBE','CASSCF','CAS','MP2','MP3','MP4','MP5','B2PLYP',
                'B2PLYPD','mPW2PLYP','mPW2PLYPD','QCISD','CCD',r'CCSD\(T\)',
                'CCSD','CC','QCID','BD','EPT','CBS-4M','CBS-QB3','ROCBS-QB3',
                'CBS-APNO','G1','G2','G2MP2','G3','G3MP2','G3B3','G3MP2B3','G4',
                'G4MP2','W1U','W1BD','W1RO','CIS',r'CIS\(D\)','CID','CISD','TD',
                'EOMCCSD','ZINDO','DFTB','DFTBA','GVB','CNDO','INDO','MINDO',
                'MNDO','SAC-CI']
methodsRegex = [i.upper() for i in methodsRegex]
methodsRegex = '|'.join(methodsRegex)

exchangeFunctional = ['S','XA','B','PW91','MPW','G96','PBE','OPBE','O','TPSS',
                      'BRX','PKZB','WPBEH','PBEH','LC-']
exchangeFunctional = '|'.join(exchangeFunctional)

correlationFunctional = ['VWN','VWN5','LYP','PL','P86','PW91','B95','PBE',
                         'TPSS','KCIS','BRC','PKZB','VP86','V5LYP']
correlationFunctional = '|'.join(correlationFunctional)

methodRegex = f'({shellWaveFunctions})({methodsRegex})'
compoundMethodRegex = '({0})({1})[0-9]*({2})'.format(shellWaveFunctions,
                                                 exchangeFunctional,
                                                 correlationFunctional)

method_re = re.compile(methodRegex)
multiplemethod_re = re.compile(mutlipleMethodRegex)
compoundmethod_re = re.compile(compoundMethodRegex)


###################### THE PRETTY FUNCTION ######################
def is_method(candidate:str) -> bool:
    """
    Tests if a candidate string is a valid method recognized by Gaussian.

    Parameters
    ----------
    candidate : str

    Returns
    -------
    bool

    """
    candidate = candidate.upper()

    test1 = multiplemethod_re.match(candidate)
    test2 = method_re.match(candidate)
    test3 = compoundmethod_re.match(candidate)
    if test1 is not None:
        test1 = test1[0] == candidate
    if test2 is not None:
        test2 = test2[0] == candidate
    if test3 is not None:
        test3 = test3[0] == candidate
    return bool(test1 or test2 or test3)

######################## BASIS SECTION #########################################
## Pople Bases require lowercasing for design reasons.
Pople = dict(primitives = r'sto|mc|[0-9]',
             zeta = '[0-9]{1,3}',
             diffuse = r'[\+]{0,2}',
             pol1 = r'[0-9]?(d|f|d\'|f\'){1,2}',
             pol2 = r'[0-9]?(\,p|d|p\'|d\')?')
x1 = Pople['pol1']
x2 = Pople['pol2']
Pople['polarization'] = r'\(({pol1}{pol2})\)|\*|\*\*'.format(pol1=x1,pol2=x2)
Pople_Basis = r'({primitives})\-({zeta})({diffuse})g({polarization})?'
Pople_Basis = Pople_Basis.format(**Pople)
# The resulting regex is this one but I wanted the reasoning to be written down
# (sto|mc|[0-9])\-([0-9]{1,3})([\+]{0,2})g(\([0-9]{0,1}(d|f|d\'|f\'){1,2}\,{0,1}[0-9]{0,1}(p|d|p\'|d\'){0,1}\)){0,1}

## ccBases
## Assumes lowercased, it's simpler so its not divided
cc_Basis = r'((sp|d)?t?(aug|jul|jun|may|apr){1}\-)?cc\-pv(d|t|q|5|6)z'
#((sp|d){0,1}t{0,1}(aug|jul|jun|may|apr){1}\-){0,1}cc\-pv(d|t|q|5|6)z

# UGBS
UGBS_Basis = r'ugbs[1-9](p|v|o)(2?\+{1,2})?'

# SV Family
SV_Basis = r'(def2)?(sv|tzv?|qzv?)p{0,2}' #Requires the usage of match

# Family that includes indicating the number of core electrons
core_Basis = r'(sdd|shf|sdf|mhf|mdf|mwb|oldsdd|sddall)[0-9]{1,3}'

# D95 Family
D95_Basis = 'd95v?('+Pople['polarization']+')'

# CEP family
CEP_Basis = r'cep\-?[0-9]{1,3}g('+Pople['polarization']+')'

# SHC Basis
SHC_Basis = 'shc('+Pople['polarization']+')'

# Or to hardcode all the family

Hardcoded_basis = set(['sec','lanl2mb','lanl2dz','midix','epr-ii','epr-iii',
                 'mtsmall','dgdzvp','dgdzvp2','dgtzvp','cbsb7','gen','genecp'])

# Construct the full regex
all_basis = [Pople_Basis,UGBS_Basis,cc_Basis,SV_Basis,core_Basis,
             D95_Basis,CEP_Basis,SHC_Basis]

basis_expr = '|'.join([f'({i})' for i in all_basis])
#print(basis_expr)
basis_regex = re.compile(basis_expr)
######################## THE PRETTY FUNCTION ###################################
def is_basis(candidate:str) -> bool:
    """
    Tests if a candidate string is a valid basis set recognized by Gaussian.

    Parameters
    ----------
    candidate : str

    Returns
    -------
    bool

    """
    lowercandidate = candidate.lower()
    in_hardcoded = lowercandidate in Hardcoded_basis
    match = basis_regex.match(lowercandidate) # Has to be match and not search
    if match:
        in_regex = match[0] == lowercandidate
    else:
        in_regex = False
    return  bool(in_hardcoded or in_regex)
######################### PERIODIC TABLE #######################################
# It is important that X appears the first one (enumeration starts in 0)
# here X stands for a dummy atom
items = """X
     H                                                                                                  He
    Li Be  B                                                                             C   N   O   F  Ne
    Na Mg Al                                                                            Si   P   S  Cl  Ar
     K Ca Sc                                           Ti  V Cr Mn Fe Co Ni Cu  Zn  Ga  Ge  As  Se  Br  Kr
    Rb Sr  Y                                           Zr Nb Mo Tc Ru Rh Pd Ag  Cd  In  Sn  Sb  Te   I  Xe
    Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta  W Re Os Ir Pt Au  Hg  Tl  Pb  Bi  Po  At  Rn
    Fr Ra Ac Th Pa  U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo
    """
items = items.replace('\n',' ').strip().split()
PeriodicTable:dict[str|int,str|int] = dict()
for i,Sym in enumerate(items):
    PeriodicTable[i] = Sym        #  1  -> 'H'
    PeriodicTable[str(i)] = Sym   # '1' -> 'H'
    PeriodicTable[Sym] = i        # 'H' ->  1
