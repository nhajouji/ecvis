from ecc.utils import gcd,gcd_list,discfac,primesBetween,divisors, mod_sfd
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
    d = qf_disc(qf)
    if d >=0:
        raise ValueError('Discriminant must be negative')
    matrix = IntegerSquareMatrix([[1,0],[0,1]])
    while not qf_in_fundom(qf):
        a,b,c = qf
        if a>c:
            m0 = IntegerSquareMatrix([[0,-1],[1,0]])
            matrix = m0*matrix
            qf = act_qf(qf,m0)
        elif a < abs(b):
            k = b//(2*a)
            if b % (2*a) >= a:
                k+=1
            m0 = IntegerSquareMatrix([[1,k],[0,1]])
            matrix =  m0*matrix
            qf = act_qf(qf,m0)
        elif a+b==0:
            m0 = IntegerSquareMatrix([[1,-1],[0,1]])
            matrix =  m0*matrix
            qf = act_qf(qf,m0)
        elif a==c and b<0:
            m0 = IntegerSquareMatrix([[0,-1],[1,0]])
            matrix =  m0*matrix
            qf = act_qf(qf,m0)
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


##########
# X_0(l) #
##########
def minv(m:IntegerSquareMatrix)->IntegerSquareMatrix:
    return -m+m.trace()

def find_rrep_g0(m:IntegerSquareMatrix,l:int)->IntegerSquareMatrix:
    reps = gamma_0_coset_reps(l)
    cands= [m0 for m0 in reps if ((m*minv(m0)).mat)[1][0]%l==0]
    if len(cands)!=1:
        raise ValueError('No unique rep')
    else:
        return cands[0]

def qf_to_gamma_0_fd(qf:tuple[int,int,int],l:int)->tuple[tuple[int,int,int],IntegerSquareMatrix]:
    qf0,m = qf_to_fun_dom(qf)
    if m.mat[1][0]%l == 0:
        return qf0, m
    ml = minv(find_rrep_g0(m,l))
    return act_qf(qf,ml),m*ml

def qf_mod_gamma_0(qf:tuple[int,int,int],l:int)->tuple[int,int,int]:
    return qf_to_gamma_0_fd(qf,l)[0]

def x0_endos(qf:tuple[int,int,int],l:int)->list[tuple[int,int,int]]:
    qf0 = qf_mod_gamma(qf)
    return [qf1 for qf1 in gamma_0_orb(qf0,l) if qf_mod_gamma(fricke_inv(qf1,l))==qf0]

def iso_taus_x0_l(qf,l):
    qf_reps = gamma_0_orb(qf_mod_gamma(qf),l)
    return [qf1 for qf1 in qf_reps if qf_disc(qf1)==qf_disc(qf)]

def isos_x0_l_all(d,l):
    qfs = get_qfs_all(d)
    iso_taus = {}
    for qf0 in qfs:
        qf1s = gamma_0_orb(qf0,l)
        for qf1 in qf1s:
            if qf_mod_gamma(fricke_inv(qf1,l)) in qfs:
                iso_taus[qf1]=fricke_inv(qf1,l)
    return iso_taus



###############
# Class group #
###############

def class_group_id(d:int):
    if d % 4 > 1:
        raise ValueError(f'{d} is not a discriminant')
    else:
        return (1,d%4,-(d//4))
    
def class_group_inv(qf:tuple[int,int,int])->tuple[int,int,int]:
    a,b,c = qf_mod_gamma(qf)
    return qf_mod_gamma((a,-b,c))

def get_qfs_with_primes(d):
    pmax = max(abs(d)+1,72)
    qfs = get_qfs_strict(d)
    qf0 = class_group_id(d)
    qf_to_prime = {}
    prime_cands = primesBetween(2,pmax)
    for p in prime_cands:
        qf1s = qf_isogenies_hor(qf0,p)
        if len(qf1s)>0 and qf1s[0] not in qf_to_prime:
            for qf1 in set(qf1s):
                qf_to_prime[qf1]=p
            if len(qf_to_prime)== len(qfs):
                return qf_to_prime
    return qf_to_prime

def cg_form_test(qf1,qf2):
    d = qf_disc(qf1)
    if qf_disc(qf2)!= d:
        raise ValueError('Incompatible forms')
    a1,b1,c1 = qf_mod_gamma(qf1)
    a2,b2,c2 = qf_mod_gamma(qf2)
    if a1 == a2 and b1+b2 == 0:
        return 0
    elif gcd_list([a1,a2,(b1+b2)//2])==1:
        return 1
    else:
        return -1
    
def cgf_get_b(qf1,qf2):
    d = qf_disc(qf1)
    if qf_disc(qf2)!=d:
        raise ValueError('Incompatible discriminants')
    a1,b1,c1 = qf1
    a2,b2,c2 = qf2
    bs = [b for b in range(2*a1*a2) if (b1-b)% (2*a1)==0 and (b2-b)% (2*a2)==0 and (b*b-d) % (4*a1*a2)==0]
    if len(bs)!= 1:
        raise ValueError('Something wrong')
    return bs[0]

def cgf_mult_table_data(d):
    qfs = get_qfs_strict(d)
    qf0 = class_group_id(d)
    if qf0 not in qfs:
        raise ValueError('Check qfs and disc')
    mult_table= {}
    for i, qf1 in enumerate(qfs):
        for qf2 in qfs[i:]:
            t = cg_form_test(qf1,qf2)
            if t == 0:
                mult_table[(qf1,qf2)]=qf0
                mult_table[(qf2,qf1)]=qf0
            elif t==1:
                b = cgf_get_b(qf1,qf2)
                a12 = qf1[0]*qf2[0]
                c12 = (b*b-d)//(4*a12)
                mult_table[(qf1,qf2)]=qf_mod_gamma((a12,b,c12))
                if qf1!= qf2:
                    mult_table[(qf2,qf1)]= mult_table[(qf1,qf2)]
    return mult_table


def completed_mat(mat):
    for row in mat:
        if '*' in row:
            return False
    return True

def extend_mat(mat):
    n = len(mat)
    ijs_unknown = [(i,j) for i in range(n) for j in range(n) if mat[i][j]=='*']
    mat0 = mat.copy()
    mat1= mat.copy()
    ch = 1
    while len(ijs_unknown)>0 and ch == 1:
        for ij in ijs_unknown:
            i, j = ij
            cands = [k for k in range(n) if k not in mat1[i]+mat1[j]]
            if len(cands)==1:
                mat1[i][j]=cands[0]
                mat1[j][i]=cands[0]
        ijs_unknown = [ij for ij in ijs_unknown if mat1[ij[0]][ij[1]]=='*']
        if mat0 == mat1:
            ch = 0
        else:
            mat0 = mat1.copy()
    return mat1


def cgf_mult_table_matdic(d):
    qf0 = class_group_id(d)
    qfs_all = get_qfs_strict(d)
    qf_to_i = {qf0:0}
    n = 1
    for qf in qfs_all:
        if qf not in qf_to_i:
            qf_to_i[qf] = n
            n+=1
    data = cgf_mult_table_data(d)
    mat = [['*' for _ in range(n)] for _ in range(n)]
    for pair in data:
        qf1, qf2 = pair
        i,j = qf_to_i[qf1],qf_to_i[qf2]
        mat[i][j] = qf_to_i[data[pair]]
    return extend_mat(mat), qf_to_i

class QFClassGroup:
    def __init__(self,d:int):
        if d % 4 >1 or d>=0:
            raise ValueError(f'{d} is not a discriminant')
        self.disc = d
        self.dc = discfac(d)
        self.fundisc = (self.dc)[0]
        self.cond = (self.dc)[1]
        self.qfs = get_qfs_strict(d)
        self.identity = class_group_id(d)
        self.qf_to_iso_degs = get_qfs_with_primes(d)
        self.mult_table_matdic = cgf_mult_table_matdic(d)
        self.mult_table_mat = (self.mult_table_matdic)[0]
        self.qf_to_index_dic = (self.mult_table_matdic)[1]
        self.mt_complete = completed_mat(self.mult_table_mat)
        self.mult_table_dic = {}
        for i1,qf1 in enumerate(self.qfs):
            for i2,qf2 in enumerate(self.qfs):
                i12= (self.mult_table_mat)[i1][i2]
                if isinstance(i12,int):
                    self.mult_table_dic[(qf1,qf2)]= self.qfs[i12]
                    self.mult_table_dic[(qf2,qf1)]= self.qfs[i12]

    def multiply_qfs_dic(self,qf1,qf2):
        qf1, qf2 = qf_mod_gamma(qf1),qf_mod_gamma(qf2)
        if qf1 not in self.qfs or qf2 not in self.qfs:
            raise ValueError('Forms not in class groups')
        elif (qf1,qf2) not in self.mult_table_dic:
            return '*'
        else:
            return self.mult_table_dic[(qf1,qf2)]
        
    def qf_power(self,qf,n):
        if n < 0:
            n = -n
            qf = class_group_inv(qf)
        qfn = class_group_id(self.disc)
        qf2n = qf
        while n > 0:
            if n % 2 == 0:
                qfn = self.multiply_qfs_dic(qfn,qf2n)
            n = n//2
            qf2n = self.multiply_qfs_dic(qf2n,qf2n)
        return qf2n

    def multiply_qfs_iso(self,qf1,qf2):
        qf1, qf2 = qf_mod_gamma(qf1),qf_mod_gamma(qf2)
        if qf1 not in self.qfs or qf2 not in self.qfs:
            raise ValueError('Forms not in class groups')
        qf0 = class_group_id(self.disc)
        # Handle trivial cases
        if qf1 == qf0:
            return qf2
        elif qf2 == qf0:
            return qf1
        elif qf1 == class_group_inv(qf2):
            return qf0
        # Make sure we have primes
        if qf1 not in self.qf_to_iso_degs or qf2 not in self.qf_to_iso_degs:
            raise ValueError('Search for primes first')
        l1 = self.qf_to_iso_degs[qf1]
        l2 = self.qf_to_iso_degs[qf2]
        # Deal with 2 torsion
        if qf1==class_group_inv(qf1):
            return list(set(qf_isogenies_hor(qf2,l1)))[0]
        elif qf2==class_group_inv(qf2):
            return list(set(qf_isogenies_hor(qf1,l2)))[0]
        # Check if one of elements is in subgroup
        # generated by other
        qf1cycle = qf_iso_cycle(qf1,l1)
        if qf2 in qf1cycle:
            i2 = [i for i in range(len(qf1cycle)) if qf1cycle[i]==qf2][0]
            return [qf1cycle[(i2+1)%len(qf1cycle)]][0]
        qf2cycle = qf_iso_cycle(qf2,l2)
        if qf1 in qf2cycle:
            i1 = [i for i in range(len(qf2cycle)) if qf2cycle[i]==qf1][0]
            return [qf2cycle[(i1+1)%len(qf2cycle)]][0]
        else:
            return [qf for qf in qf_isogenies_hor(qf1,l2) if qf in qf_isogenies_hor(qf2,l1)] [0]
