from ecc.nt import gcd,gcd_list,discfac,primesBetween,divisors, mod_sfd,quad_rec
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

def qf_2_mat(qf):
    a,b,c = qf
    return [[0,-c],[a,-b]]
def mat_2_qf(m):
    a,b,c,d = m[0][0],m[0][1],m[1][0],m[1][1]
    s = c//abs(c)
    return (s*c,s*(a-d),-s*b)

def prod_tup(t):
    p = 1
    for x in t:
        p*=x
    return p

### Computing isogeny codomains
def qf_isogs(qf0,l):
    a,b,c = qf0
    qfls = []
    if c % l == 0:
        qfls.append(qf_mod_gamma((a*l,b,c//l)))
    for t in range(l):
        qt = a+b*t+c*t*t
        if qt % l == 0:
            at = qt//l
            bt = (b+2*c*t)
            ct = c*l
            qfls.append(qf_mod_gamma((at,bt,ct)))
    return qfls

# Computing isogeny cycles
def qf_isog_cycle(qf0,l):
    cyc = qf_isogs(qf0,l)
    if len(cyc)==1:
        return [qf0,cyc[0]]
    elif len(cyc)>2:
        raise ValueError('Too many isogenies')
    cyc = [qf0,cyc[0]]
    nextbatch = [qf for qf in qf_isogs(cyc[-1],l) if qf not in cyc]
    while len(nextbatch)>0:
        cyc.append(nextbatch[0])
        nextbatch = [qf for qf in qf_isogs(cyc[-1],l) if qf not in cyc]
    return cyc

def qf_isog_cycle_power(qf0,lk):
    l,k = lk
    if k < 0:
        return qf_sibs(qf0,l)
    elif k == 0:
        return [qf0]
    cyc = qf_isog_cycle(qf0,l)
    if k == 1:
        return cyc
    n = len(cyc)
    m = gcd(n,k)
    nm = n//m
    return [cyc[(k*i) % n] for i in range(nm)]

def qf_sibs(qf0:tuple[int,int,int],l:int):
    sibs = qf_isogenies_down(qf_parents(qf0,l)[0],l)
    return [qf0]+[qf for qf in sibs if qf != qf0]

def cycs_from_ancestors(qf0):
    d, c = discfac(qf_disc(qf0))
    cycs = {}
    if c % 2 == 0:
        sibs = qf_sibs(qf0,2)
        if len(sibs)>1:
            cycs[2] = sibs
    if c % 3 == 0:
        sibs = qf_sibs(qf0,3)
        if len(sibs)<4:
            cycs[3] = sibs
    return cycs

def generate_qfs_from_ls(qf0,lset):
    qfs = [qf0]
    for l in lset:
        qfls = []
        for qf in qfs:
            qfls+=qf_isog_cycle(qf,l)
        qfs = qfls
    return qfs

def generate_qfs_from_lks(qf0,lkset):
    qfs = [qf0]
    for lk in lkset:
        qfls = []
        for qf in qfs:
            qfls+=qf_isog_cycle_power(qf,lk)
        qfs = qfls
    return qfs

def qf_binary_tree(qf0,lset):
    rows = [[qf0]]
    for l in lset:
        lastrow = rows[-1]
        newrow = []
        for qf in lastrow:
            newrow+= qf_isogs(qf,l)
        rows.append(newrow)
    return rows

def qf_binary_tree_chains(qf0,lset):
    chains = [[qf0]]
    for l in lset:
        newchains = []
        for chain in chains:
            lastqf = chain[-1]
            qfls = qf_isogs(lastqf,l)
            for qf in qfls:
                newchains.append(chain+[qf])
        chains = newchains
    return chains



def qf_cyc_data(d,lcands):
    d0, c = discfac(d)
    qf0 = class_group_id(d)
    lcands = [l for l in lcands if c % l !=0 and quad_rec(d0,l)>=0]
    cycdata = {}
    gens = []
    for l in lcands:
        cycl = qf_isog_cycle(qf0,l)
        if len(cycl)>1 and cycl[1] not in gens and cycl[-1] not in gens:
            cycdata[l] = len(cycl)
            gens.append(cycl[1])
    return cycdata


def qf_cyc_data_ext(d,lcands):
    cycdata = qf_cyc_data(d,lcands)
    ls = [l for l in cycdata]
    for l in ls:
        for k in divisors(cycdata[l])[:-1]:
            cycdata[(l,k)] = cycdata[l]//k
    anc_data= cycs_from_ancestors(class_group_id(d))
    for l in anc_data:
        cycdata[(l,-1)]=len(anc_data[l])
    return {l:cycdata[l] for l in cycdata if type(l)==tuple}

def subs_w_prod(lscores:dict,target:int):
    incomp = {():1}
    complete = []
    for l in lscores:
        newdata = {}
        for lt in incomp:
            ltx = tuple(list(lt)+[l])
            n_ltx= lscores[l]*incomp[lt]
            if n_ltx == target:
                complete.append(ltx)
            elif n_ltx<target:
                newdata[ltx] = n_ltx
        incomp.update(newdata)
    return complete

def qf_search_lgens(d,lcands):
    cycdata = qf_cyc_data_ext(d,lcands)
    cld = len(get_qfs_strict(d))
    if cld in cycdata.values():
        return [(l,) for l in cycdata if cycdata[l]==cld]
    ltups_all = subs_w_prod(cycdata,cld)
    # Get rid of combinations that involve the same prime multiple times
    ltups_all= [lks for lks in ltups_all if len(lks)==len(set([lk[0] for lk in lks]))]
    ltups_all.sort(key = len)
    lks_gen = [lk for lk in ltups_all if len(set(generate_qfs_from_lks(class_group_id(d),lk)))==cld]
    if len(lks_gen)==0:
        return []
    ml = min([len(lks) for lks in lks_gen])
    return [lks for lks in lks_gen if len(lks) == ml]

####
# TODO: Replace qf_iso_cycle below with qf_isog_cycle

def qf_iso_cycle(qf0:tuple[int,int,int],l:int):
    cycle = []
    nextbatch = [qf0]
    while len(nextbatch)>0:
        qf = nextbatch[0]
        cycle.append(qf)
        nextbatch = [qf1 for qf1 in qf_isogenies_hor(qf,l) if qf1 not in cycle]
    return cycle

def qf_iso_cycle_oriented(qf0:tuple[int,int,int],qf1:tuple[int,int,int],l):
    cycle = qf_iso_cycle(qf0,l)
    if len(cycle)<2:
        raise ValueError('No cycle in this degree')
    if cycle[1]==qf1:
        return cycle
    elif cycle[-1]==qf1:
        return [qf0]+cycle[:0:-1]
    else:
        raise ValueError('Forms are not neighbors in the graph')
    
def intersect_qf_codoms(qf0,qf1,l0,l1):
    return [qf2 for qf2 in qf_isogenies_hor(qf0,l0) if qf2 in qf_isogenies_hor(qf1,l1)]

def qf_iso_frame(qf1,qfa,la,lb,lab):
    intab = intersect_qf_codoms(qf1,qfa,lab,lb)
    if len(intab)!= 1:
        raise ValueError(f'Intersection for ab has size {len(intab)}')
    qfab = intab[0]
    intb = intersect_qf_codoms(qf1,qfab,lb,la)
    if len(intb)!= 1:
        raise ValueError(f'Intersection for b has size {len(intb)}')
    return intb[0],qfab

def qf_iso_mat_from_frame(qf1,qfa,la,lb,lab):
    qfb,qfab = qf_iso_frame(qf1,qfa,la,lb,lab)
    col0 = qf_iso_cycle_oriented(qf1,qfb,lb)
    col1 = qf_iso_cycle_oriented(qfa,qfab,lb)
    mat = []
    for i, qfbi in enumerate(col0):
        qfabi = col1[i]
        mat.append(qf_iso_cycle_oriented(qfbi,qfabi,la))
    return mat

def qf_isomat_ext(mat,l3):
    return [[qf_iso_cycle(qf,l3)[-1] for qf in row] for row in mat]


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

def qf_x0_endos(qf:tuple[int,int,int],l:int)->list[tuple[int,int,int]]:
    qf0 = qf_mod_gamma(qf)
    return [qf1 for qf1 in gamma_0_orb(qf0,l) if qf_mod_gamma(fricke_inv(qf1,l))==qf0]

def x0_endos_all(p:int)->dict:
    endos_by_trace = {}
    a = 0
    while a*a < 4*p:
        d = a*a-4*p
        qfs = get_qfs_all(d)
        endos_by_trace[a]={qf:qf_x0_endos(qf,p) for qf in qfs}
        a+=1
    return endos_by_trace

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
