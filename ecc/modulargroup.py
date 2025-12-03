from ecc.utils import gcd,discfac,int_sqrt
from ecc.ringclasses import IntegerSquareMatrix
from ecc.modularpolynomials import *

###################
# Quadratic Forms #
###################


def qf_in_fundom(qf:tuple[int,int,int])->bool:
    a,b,c = qf
    if abs(b)< a and a < c:
        return True
    elif b == a and a < c:
        return True
    elif c == a and b >= 0 and b <= a:
        return True
    return False

def qf_gcd(qf:tuple[int,int,int])->int:
    a,b,c = qf
    return gcd(a,gcd(b,c))

def qf_is_prim(qf:tuple[int,int,int])->bool:
    return qf_gcd(qf)==1

def qf_make_prim(qf:tuple[int,int,int])->tuple[int,int,int]:
    a,b,c = qf
    g = gcd(a,gcd(b,c))
    if g > 1:
        a,b,c = a//g,b//g,c//g
    return (a,b,c)

def qf_disc(qf:tuple[int,int,int])->tuple[int,int,int]:
    a,b,c = qf_make_prim(qf)
    return b*b-4*a*c

def qf_to_dc(qf:tuple[int,int,int])->tuple[int,int]:
    return discfac(qf_disc(qf))

########################
# Modular group action #
########################

def qf_to_mat(qf:tuple[int,int,int])->IntegerSquareMatrix:
    a,b,c = qf
    return IntegerSquareMatrix([[2*a,b],[b,2*c]])

def mat_to_qf(m:IntegerSquareMatrix)->tuple[int,int,int]:
    if not isinstance(m,IntegerSquareMatrix) or m.dim!=2:
        raise TypeError('Input should be 2x2 integer matrix')
    arr = m.mat
    a,b1,b2,c = arr[0][0],arr[0][1],arr[1][0],arr[1][1]
    if b1!= b2:
        raise ValueError('Input should be symmetric matrix')
    if a % 2 !=0 or c %2 != 0:
        raise ValueError('Diagonal entries should be even')
    return (a//2,b1,c//2)

def act_qf(qf:tuple[int,int,int],m:IntegerSquareMatrix):
    qfm = qf_to_mat(qf)
    tm = m.trace()
    minv = IntegerSquareMatrix([[tm,0],[0,tm]])+(-1)*m
    qfm_new = minv.transpose()*qfm*minv
    return mat_to_qf(qfm_new)


def qf_to_fun_dom(qf:tuple)->tuple:
    matrix = IntegerSquareMatrix([[1,0],[0,1]])
    while not qf_in_fundom(qf):
        a,b,c = qf
        if a>c:
            m0 = IntegerSquareMatrix([[0,-1],[1,0]])
            matrix = matrix*m0
            qf = act_qf(qf,m0)
            if qf[0]>qf[2]:
                return 'Step 1 failed'
        elif a < abs(b):
            k = b//(2*a)
            if b % (2*a) >= a:
                k+=1
            m0 = IntegerSquareMatrix([[1,k],[0,1]])
            matrix = matrix*m0
            qf = act_qf(qf,m0)
            if qf[0]<abs(qf[1]):
                return 'Step 2 failed'
        elif a+b==0:
            m0 = IntegerSquareMatrix([[1,-1],[0,1]])
            matrix = matrix*m0
            qf = act_qf(qf,m0)
            if qf[0]+qf[1]==0:
                return 'Step 3 failed'
        elif a==c and b<0:
            m0 = IntegerSquareMatrix([[0,-1],[1,0]])
            matrix = matrix*m0
            qf = act_qf(qf,m0)
            if qf[1]<0:
                return 'Step 4 failed'
        else:
            return qf,matrix
    return qf, matrix

def qf_mod_gamma(qf:tuple[int,int,int])->tuple[int,int,int]:
    return qf_make_prim(qf_to_fun_dom(qf)[0])

#######################################
# Generating lists of quadratic forms #
#######################################

def get_qfs_all(d:int):
    reps_found = []
    # First we check that d is indeed a discriminant;
    # if it isn't, we simply return an empty list because there are no associated lattices
    if d % 4 > 1 or d >= 0:
        return reps_found
    b = d % 4
    while 3*b *b <= abs(d):
        num = (b*b-d)//4
        a = b
        while a *a <= num:
            if a == 0:
                a+=1
            if num % a == 0:
                c = num//a
                if qf_in_fundom((a,b,c)):
                    reps_found.append((a,b,c))
                if b!= 0 and qf_in_fundom((a,-b,c)):
                    reps_found.append((a,-b,c))
            a+=1
        b+=2
    return reps_found

def get_qfs_strict(d:int):
    return [qf for qf in get_qfs_all(d) if qf_disc(qf)==d]


#############
# Isogenies #
#############
def fricke_inv(qf:tuple[int,int,int],l:int)->tuple[int,int,int]:
    a,b,c = qf
    return (l*l * c, -l*b,a)

def gamma_0_coset_reps(p:int)->list[IntegerSquareMatrix]:
    return [IntegerSquareMatrix([[1,0],[0,1]])]+[IntegerSquareMatrix([[0,-1],[1,a]])
           for a in range(-(p//2),(p//2)+(p%2))]

def gamma_0_orb(qf:tuple[int,int,int],l:int)->list[tuple[int,int,int]]:
    return [act_qf(qf,m) for m in gamma_0_coset_reps(l)]


def qf_isogenies_all(qf:tuple[int,int,int],l:int,normalize=True):
    isoqfs = [fricke_inv(qf0,l) for qf0 in gamma_0_orb(qf,l)]
    if not normalize:
        return isoqfs
    return [qf_mod_gamma(qf0) for qf0 in isoqfs]

def qf_isogenies_hor(qf:tuple[int,int,int],l:int):
    d = qf_disc(qf)
    return [qf0 for qf0 in qf_isogenies_all(qf,l) if qf_disc(qf0)==d]

def qf_isogenies_down(qf:tuple[int,int,int],l:int):
    d = qf_disc(qf)
    return [qf0 for qf0 in qf_isogenies_all(qf,l) if qf_disc(qf0)<d]

def qf_parents(qf:tuple[int,int,int],l:int):
    d = qf_disc(qf)
    return [qf0 for qf0 in qf_isogenies_all(qf,l) if qf_disc(qf0)>d]

def qf_iso_cycle(qf0:tuple[int,int,int],l:int):
    cycle = []
    nextbatch = [qf0]
    while len(nextbatch)>0:
        qf = nextbatch[0]
        cycle.append(qf)
        nextbatch = [qf1 for qf1 in qf_isogenies_hor(qf,l) if qf1 not in cycle]
    return cycle

def d_to_rqf(d:int)->tuple:
    if d >= 0 or (d%4 >1):
        raise ValueError('Input must be a negative discriminant')
    if d % 4 == 0:
        return (1,0,(-d)//4)
    else:
        return (1,1,-(d//4))
    
def d_to_ssl_cycle_data(d:int)->dict:
    qf0 = d_to_rqf(d)
    qf0_data = {l:qf_iso_cycle(qf0,l) for l in atkin_polys_dict}
    ns = [len(qf0_data[l]) for l in qf0_data]
    nmx = max(ns)
    l0 = min([l for l in qf0_data if len(qf0_data[l])==nmx])
    data = {qf0:{l0:qf0_data[l0]}}
    for qf1 in data[qf0][l0]:
        if qf1 not in data:
            data[qf1]={}
    # look for cycles of length 2 among remaining l's
    qf2s_seen = [qf for qf in data]
    for l1 in qf0_data:
        if l1 != l0 and len(qf0_data[l1])==2:
            qf2 = qf0_data[l1][1]
            if qf2 not in qf2s_seen:
                data[qf0][l1]=qf2
                qf2s_seen.append(qf2)
                for qf1 in data:
                    if qf1!=qf0:
                        qf2s = qf_iso_cycle(qf1,l1)
                        data[qf1][l1]=qf2s[-1]
    return data
###########
# Classes #
###########

class CMLattice:
    def __init__(self,qf:tuple[int,int,int]):
        self.qf = qf
        self.disc = qf_disc(qf)
        self.lc = qf[0]
        self.height_ub = 1+((int_sqrt(-qf_disc(qf))+1)//(2*qf[0]))
    