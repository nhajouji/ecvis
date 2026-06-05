## Frobenius matrix
# The following gives a matrix that represents
# the action of Frobenius on the lattice,
# relative to the basis 1, tau.
def frobmat(ap:tuple[int],abc:tuple[int])->IntegerSquareMatrix:
    t,p = ap
    a,b,c = qf_make_prim(abc)
    df,cf = discfac(t*t-4*p)
    dt,ct = discfac(b*b-4*a*c)
    if dt!= df or cf % ct != 0:
        raise ValueError('incompatible discriminants') 
    cft = cf // ct
    trdiff = t+b*cft
    if trdiff%2 != 0:
        raise ValueError('Check trace')
    t0 = trdiff//2
    return IntegerSquareMatrix([[t0,-a*c*cft],[cft,t0-b*cft]])

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
    mat = IntegerSquareMatrix(cmat.mat).mat
    cdet = mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1]
    if abs(cdet) == 1:
        return {(1,0):m,(0,1):m}
    elif m == 1:
        return kernel_gen_cyc(cmat.mat)
    else:
        return divide_cyclic_gen(kernel_gen_cyc(cmat.mat),m)

    