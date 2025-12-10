from ecc.utils import quad_rec, gcd, hall_multiplier,axby,discfac, primefact
from ecc.modularpolynomials import *

def red_bf(abc:tuple)->bool:
    a,b,c = abc
    if abs(b)< a and a < c:
        return True
    elif b == a and a < c:
        return True
    elif c == a and b >= 0 and b <= a:
        return True
    return False

# Input: a discriminant d
# Output: The list of reduced binary quadratic forms with discriminant d

def get_cl_reps(d:int)->list:
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
                if red_bf((a,b,c)):
                    reps_found.append((a,b,c))
                if b!= 0 and red_bf((a,-b,c)):
                    reps_found.append((a,-b,c))
            a+=1
        b+=2
    return reps_found


def get_cl_reps_filtered(d:int)->dict:
    reps = get_cl_reps(d)
    ds = {}
    for abc in reps:
        a,b,c = abc
        g = gcd(a,gcd(b,c))
        a,b,c = a//g,b//g,c//g
        ds[(a,b,c)]=d//(g*g)
    grps = {ds[abc]:[] for abc in ds}
    for abc in ds:
        grps[ds[abc]].append(abc)
    return grps



## Frobenius


class Matrix2x2:
    def __init__(self, data):
        if not isinstance(data, (list, tuple)) or len(data) != 2:
            raise ValueError("Input must be a 2x2 structure (list of lists or tuple of tuples)")
        for row in data:
            if not isinstance(row, (list, tuple)) or len(row) != 2:
                raise ValueError("Each row must have 2 elements")
            for element in row:
                if not isinstance(element, int):
                    raise TypeError("All elements must be integers")
        self.data = [list(row) for row in data]  # Store as a list of lists

    def __str__(self):
        return f"[{self.data[0][0]}, {self.data[0][1]}]\n[{self.data[1][0]}, {self.data[1][1]}]"

    def __add__(self, other):
        if not isinstance(other, Matrix2x2):
            raise TypeError("Can only add with another Matrix2x2 object")
        result = [[0, 0], [0, 0]]
        for i in range(2):
            for j in range(2):
                result[i][j] = self.data[i][j] + other.data[i][j]
        return Matrix2x2(result)

    def __sub__(self, other):
        if not isinstance(other, Matrix2x2):
            raise TypeError("Can only subtract with another Matrix2x2 object")
        result = [[0, 0], [0, 0]]
        for i in range(2):
            for j in range(2):
                result[i][j] = self.data[i][j] - other.data[i][j]
        return Matrix2x2(result)

    def __mul__(self, other):
        if not isinstance(other, Matrix2x2):
            raise TypeError("Can only multiply with another Matrix2x2 object")
        result = [[0, 0], [0, 0]]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    result[i][j] += self.data[i][k] * other.data[k][j]
        return Matrix2x2(result)
    
    def mvec(self,v):
        data = self.data
        a,b,c,d = data[0][0],data[0][1],data[1][0],data[1][1]
        x,y = v[0],v[1]
        return [a*x+b*y,c*x+d*y]
    
    def det(self):
        mat = self.data
        return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]
    
    def inv(self):
        mat = self.data
        a,b,c,d = mat[0][0],mat[0][1],mat[1][0],mat[1][1]
        det = a*d-b*c
        if det**2 != 1:
            raise ValueError('Determinant must be +/- 1 to invert.')
        return Matrix2x2([[det*d,-det*b],[-det*c,det*a]])

    def gcdfac(self):
        data = self.data
        a,b,c,d = data[0][0],data[0][1],data[1][0],data[1][1]
        g = gcd(gcd(a,b),gcd(c,d))
        if g == 0:
            return self, 0
        else:
            return Matrix2x2([[a//g,b//g],[c//g,d//g]]),g

    def __pow__(self, power):
        if not isinstance(power, int):
            raise ValueError("Power must be an integer")
        if power < 0:
            power = -power
            base = self.inv()
        else:
            base = self
        result =  Matrix2x2([[1, 0], [0, 1]])  # Identity matrix
        while power > 0:
            if power % 2 == 1:
                result *=base
            power = power//2
            base *=base
        return result

    def __eq__(self, other):
        if not isinstance(other, Matrix2x2):
            return False
        return self.data == other.data

## Frobenius matrix
# The following gives a matrix that represents
# the action of Frobenius on the lattice,
# relative to the basis 1, tau.
def frobmat(ap:tuple[int],abc:tuple[int]):
    t,p = ap
    a,b,c = abc
    g = gcd(gcd(a,b),c)
    a,b,c = a//g,b//g,c//g
    df,mf = discfac(t*t-4*p)
    dt,mt = discfac(b*b-4*a*c)
    if dt!= df or mf % mt != 0:
        return 'Check disc'
    mft = mf // mt
    trdiff = t+b*mft
    if trdiff%2 != 0:
        return 'Check trace'
    t0 = trdiff//2
    return Matrix2x2([[t0,-a*c*mft],[mft,t0-b*mft]])



#############
## Kernels ##
#############

def kernel_gen_cyc(mat:list[list[int]])->dict:
    m00 = mat[0][0]
    m01 = mat[0][1]
    m10 = mat[1][0]
    m11 = mat[1][1]
    l = m00*m11 - m01*m10
    l0 = gcd(gcd(m00,l),gcd(m01,l))
    l1 = gcd(gcd(m10,l),gcd(m11,l))
    if gcd(l0,l1)>1:
        return 'Check gcds'
    n0 = hall_multiplier(l0,l)
    # Multiplying u by n0 gives us a point that generates a subgroup
    # of index l0c; the order of the point is l/l0c
    l0c = n0*l0
    # We need a point of order l0c, and index l/l0c
    # We already have a point of order l/l1
    n1 = l//(l0c*l1)
    # v already has order l/l1
    # (l/l1) *v = 0, so (l/l1)/(l0c) will have order l0c
    xg = -(n0*m01+n1*m11)%l
    yg = (n0*m00+n1*m10)%l
    # Note that gcd(xg,yg, l) should now be equal to 1.
    # If gcd(xg,yg) is not 1, we can factor it out
    gxy = gcd(xg,yg)
    if gcd(gxy,l)!= 1:
        return 'Something went wrong'
    if gxy == 0:
        return 'Something wrong'
    return {(xg//gxy,yg//gxy):l}


def divide_cyclic_gen(gen:dict,m:int)->dict:
    v = [v for v in gen][0]
    l = gen[v]
    x,y = v
    if gcd(x,y)>1:
        g =gcd(x,y)
        if gcd(g,l)>1:
            return'The generator has the wrong order'
        v = x//g, y//g
        x,y = v
    v1 = (l*x)%(l*m),(l*y)%(l*m) 
    r,s = axby(v1)
    w = ((l*s) %(l*m),(l*r)%(l*m))
    return {v:m*l,w:m}

# MW generators
def mw_gens(ap:tuple[int],abc:tuple[int],n:int)->dict:
    fmat = frobmat(ap,abc)
    cmat,m = (fmat**n - (fmat**0)).gcdfac()
    if abs(cmat.det()) == 1:
        return {(1,0):m,(0,1):m}
    elif m == 1:
        return kernel_gen_cyc(cmat.data)
    else:
        return divide_cyclic_gen(kernel_gen_cyc(cmat.data),m)

def divide_cyclic_gen_v1(gen:dict,m:int)->dict:
    v = [v for v in gen][0]
    l = gen[v]
    x,y = v
    if gcd(x,y)>1:
        g =gcd(x,y)
        if gcd(g,l)>1:
            return'The generator has the wrong order'
        x = x//g
        y = y//g
        v = (x,y)
    r,s = axby(v)
    w = (-l*s %(l*m),l*r%(l*m))
    return {v:m*l,w:m}


def divide_cyclic_gen_old(gen:dict,m:int)->dict:
    v = [v for v in gen][0]
    l = gen[v]
    x,y = v
    if gcd(x,y)>1:
        g =gcd(x,y)
        if gcd(g,l)>1:
            return'The generator has the wrong order'
        v = x//g, y//g
        x,y = v
    v1 = (l*x)%(l*m),(l*y)%(l*m) 
    r,s = axby(v1)
    w = ((-l*s) %(l*m),(l*r)%(l*m))
    return {v:m*l,w:m}

def divide_cyclic_gen(gen:dict,m:int)->dict:
    v = [v for v in gen][0]
    l = gen[v]
    x,y = v
    if gcd(x,y)>1:
        g =gcd(x,y)
        if gcd(g,l)>1:
            return'The generator has the wrong order'
        v = x//g, y//g
        x,y = v
    v1 = (l*x)%(l*m),(l*y)%(l*m) 
    r,s = axby(v1)
    w = ((l*s) %(l*m),(l*r)%(l*m))
    return {v:m*l,w:m}

# MW generators
def mw_gens(ap:tuple[int],abc:tuple[int],n:int)->dict:
    fmat = frobmat(ap,abc)
    cmat,m = (fmat**n - (fmat**0)).gcdfac()
    if abs(cmat.det()) == 1:
        return {(1,0):m,(0,1):m}
    elif m == 1:
        return kernel_gen_cyc(cmat.data)
    else:
        return divide_cyclic_gen(kernel_gen_cyc(cmat.data),m)
    

def gens_to_list(gens:dict)->list:
    den = max(gens.values())
    pts = [(0,0)]
    for g in gens:
        ge = gens[g]
        gx, gy = g
        pts_new = []
        for pt in pts:
            ptx,pty = pt
            pts_new+=[((ptx+e*gx)%den,(pty+e*gy)%den)
                      for e in range(ge)]
        pts = pts_new
    return pts


import numpy as np

def abc_to_tau(abc:tuple[int,int,int])->np.array:
    tau_np= np.roots(abc)[0]
    return np.array([tau_np.real,tau_np.imag])

def mw_arr_from_gens(gens:dict,abc:tuple):
    den = max(gens.values())
    one = np.array([1,0])
    tau = abc_to_tau(abc)
    pts = [np.array([0,0])]
    for gen in gens:
        x,y = gen
        pts = [(pt0+m*np.array([x,y]))%den 
               for pt0 in pts for m in range(gens[gen])]
    return np.array([pt[0]*one+pt[1]*tau for pt in pts])/den

def mw_arr(ap:tuple[int,int],abc:tuple[int,int,int],n:int):
    gens_dict = mw_gens(ap,abc,n)
    den = max(gens_dict.values())
    one = np.array([1,0])
    tau = abc_to_tau(abc)
    pts = [np.array([0,0])]
    for gen in gens_dict:
        x,y = gen
        pts = [(pt0+m*np.array([x,y]))%den 
               for pt0 in pts for m in range(gens_dict[gen])]
    return np.array([pt[0]*one+pt[1]*tau for pt in pts])/den
    

### Exporting

def export_points(grp:list,filename:str):
    file = open(filename,'a')
    for i in range((len(grp)//3)+1):
        stri = ''
        for vec in grp[3*i:3*(i+1)]:
            stri += '['+str(vec[0])[:10]+','+str(vec[1])[:10]+'],'
        stri+='\n'
        file.write(stri)
    file.close()


### Testing

def mwgens_test(ap:tuple[int],abc:tuple[int],n:int):
    fmat = frobmat(ap,abc)
    bmat = fmat**n - fmat**0
    gens = mw_gens(ap,abc,n)
    den = max(gens.values())
    results = {}
    for v in gens:
        fv = bmat.mvec(v)
        results[v] = ((fv[0])%den,fv[1]%den)
    return results

def mwsize_test(gens:dict):
    expected_size = 1
    for e in gens.values():
        expected_size*=e
    pts = gens_to_list(gens)
    return {expected_size,len(set(pts))}



## Moduli ##

### Quadratic forms

def clear_common_fac_bqf(abc:tuple[int,int,int]):
    a,b,c  = abc
    g = gcd(gcd(a,b),c)
    if g == 0:
        return abc
    else:
        return a//g,b//g,c//g

def bqf_disc(abc:tuple[int,int,int])->int:
    a,b,c = clear_common_fac_bqf(abc)
    return b**2-4*a*c


def bqf_fr(abc:tuple,p:int)->tuple:
    a,b,c = abc
    return (p*p*c,-p*b,a)

def sym2(m2:np.matrix)->np.matrix:
    a = m2[0,0]
    b = m2[0,1]
    c = m2[1,0]
    d = m2[1,1]
    return np.array([
        [a*a,a*c,c*c],
         [2*a*b,b*c+a*d,2*c*d],
         [b*b,b*d,d*d]
         ])

def sym2D(m2:np.array)->np.array:
    a = m2[1,1]
    b = -m2[0,1]
    c = -m2[1,0]
    d = m2[0,0]
    return np.array([
        [a*a,a*c,c*c],
         [2*a*b,b*c+a*d,2*c*d],
         [b*b,b*d,d*d]
         ])

def act_bqf(bqf:tuple,mat:np.array)->tuple:
    a,b,c = tuple(np.matmul(sym2D(mat),bqf))
    return int(a),int(b),int(c)

def gamma0orb(bqf:tuple,p:int):
    mats = [np.array([[1,0],[0,1]])]+[np.array([[0,-1],[1,a]]) for a in range(-(p//2),(p//2)+(p%2))]
    return [act_bqf(bqf,mat) for mat in mats]

def reduce(bqf:tuple)->tuple:
    matrix = np.array([[1,0],[0,1]])
    while True:
        if bqf[0] > bqf[2]:
            m0 = np.array([[0,-1],[1,0]])
            matrix = np.matmul(m0,matrix)
            bqf = act_bqf(bqf,m0)
            if bqf[0]>bqf[2]:
                return 'something wrong case 1'
        elif bqf[0] < abs(bqf[1]):
            a = bqf[0]
            b = bqf[1]
            k = b//(2*a)
            if b % (2*a) >= a:
                k+=1
            m0 = np.array([[1,k],[0,1]])
            matrix = np.matmul(m0,matrix)
            bqf = act_bqf(bqf,m0)
            if bqf[0]<abs(bqf[1]):
                return 'something wrong case 2'
        elif bqf[0]+bqf[1]==0:
            m0 = np.array([[1,-1],[0,1]])
            matrix = np.matmul(m0,matrix)
            bqf = act_bqf(bqf,m0)
            if bqf[0]+bqf[1]==0:
                return 'something wrong case 3'
        elif bqf[0]==bqf[2] and bqf[1]<0:
            m0 = np.array([[0,-1],[1,0]])
            matrix = np.matmul(m0,matrix)
            bqf = act_bqf(bqf,m0)
            if bqf[1]<0:
                return 'something wrong case 4'
        else:
            return tuple(bqf),matrix
        

def bqf_iso_verts(abc:tuple[int,int,int],l:int):
    return [clear_common_fac_bqf(reduce(bqf_fr(tup,l))[0]) for tup in gamma0orb(abc,l)]

def bqf_parents(abc:tuple[int,int,int],l:int):
    d = bqf_disc(abc)
    return list({qf for qf in bqf_iso_verts(abc,l) if bqf_disc(qf)*(l**2)==d})
    
def bqf_iso_nextrow(abc:tuple[int,int,int],l:int):
    d0 = bqf_disc(abc)
    return list(set([tup for tup in bqf_iso_verts(abc,l) if bqf_disc(tup)<d0]))

def bqf_horiz_isogs(abc:tuple[int,int,int],l:int):
    return [bqf for bqf in bqf_iso_verts(abc,l) if bqf_disc(abc)==bqf_disc(bqf)]

def get_cycles_qf_fromd(d:int,l:int):
    cycles = []
    qfs = [abc for abc in get_cl_reps(d) if bqf_disc(abc)==d]
    while len(qfs)>0:
        cycle = []
        nextbatch = [qfs[0]]
        while len(nextbatch)>0:
            q0 = nextbatch[0]
            cycle.append(q0)
            nextbatch = [q for q in bqf_horiz_isogs(q0,l) if q not in cycle]
        cycles.append(cycle)
        qfs = [q for q in qfs if q not in cycle]
    return cycles

def check_ssl_gen_qf(d:int):
    d0, c = discfac(d)
    cycs = {'Solved':False}
    for l in atkin_polys_dict:
        if c % l != 0 and quad_rec(d0,l)==1:
            cycs[l] = get_cycles_qf_fromd(d,l)
            if len(cycs[l])==1:
                return {l:cycs[l][0],'Solved':True}
    return cycs

def iso_tree_dm_labels(d:int,m:int):
    v0s = get_cl_reps(d)
    labelsdic = {}
    for i,abc in enumerate(v0s):
        labelsdic[chr(i+ord('a'))] = abc
    pfm = primefact(m)
    for p in pfm:
        lastbatch = [s for s in labelsdic]
        for e in range(pfm[p]):
            newbatch = []
            for s in lastbatch:
                abc = labelsdic[s]
                nextrow = bqf_iso_nextrow(abc,p)
                for j, abc1 in enumerate(nextrow):
                    labelsdic[s+str(p)+chr(j+ord('a'))] = abc1
                    newbatch.append(s+str(p)+chr(j+ord('a')))
            lastbatch = newbatch
    return labelsdic


### Export all point groups

def ptgroupdata_all(ap:tuple[int,int],n:int):
    a,p = ap
    d,m = discfac(a*a-4*p)
    labels = iso_tree_dm_labels(d,m)
    data = {}
    for s in labels:
        abc = labels[s]
        gens = mw_gens(ap,abc,n)
        den = max(gens.values())
        data[s] = {'tau_minpoly':abc,
                   'generators':gens,
                   'denominator':den}
    return data
