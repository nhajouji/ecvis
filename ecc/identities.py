from ecc.nt import gcd,quad_rec,gen_quad_symb,primefact,divisors,discfac,no_odd_prime_facs

def twisted_phi(d:int,m:int)->int:
    if m < 0:
        return 0
    pfm = primefact(m)
    phim = 1
    for p in pfm:
        phim*=(p-quad_rec(d,p))*(p**(pfm[p]-1))
    if d == -3 and m > 1:
        return phim//3
    elif d == - 4 and m > 1:
        return phim//2
    else:
        return phim

def twisted_phi_sum(d:int,m:int)->int:
    divsm = divisors(m)
    return sum([twisted_phi(d,m) for m in divsm])


def clgr_size_fd(d):
    if abs(d)<13:
        return 1
    den = 2-quad_rec(d,2)
    return sum([gen_quad_symb(d,l) for l in range(1,(abs(d)+1)//2) if gcd(l,d)==1])//den

def clgr_size_gen(dc):
    if abs(dc)<13:
        return 1
    d,c = discfac(dc)
    return clgr_size_fd(d)*twisted_phi(d,c)


def clgr_2len(d):
    if d % 4 == 1:
        return no_odd_prime_facs(d)
    n = d//(-4)
    r = no_odd_prime_facs(d)
    if n % 4 == 3:
        return r
    elif n % 4 in [1,2] or n % 8 == 4:
        return r+1
    else:
        return r+2
    
def clgr2_size(d):
    return 2**(clgr_2len(d)-1)