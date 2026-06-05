from qfs import *
###############
# Class group #
###############

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
