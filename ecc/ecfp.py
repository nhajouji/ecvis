from ecc.utils import quad_rec,discfac,divisors
from ecc.ringclasses import *
from ecc.lattices import get_cl_reps
from ecc.modularpolynomials import *

def j_to_fg(j:int,char = 0):
    if j == 0:
        return (0,1)
    elif j == 1728 or char>0 and (j-1728)%char == 0:
        return (1,0)
    else:
        f = -3*j*(j-1728)
        g = 2*j*((j-1728)**2)
        if char == 0:
            return (f,g)
        else:
            return (f % char, g% char)
        
def fg_to_j(fg:tuple[int,int],char =0):
    f,g = fg
    if f == 0 or char>0 and f%char ==0:
        return 0
    elif g == 0 or char>0 and g % char ==0:
        return 1728
    else:
        f3 = 4*(f**3)
        jnum = 1728 *f3
        jden = f3+27*(g**2)
        if char == 0:
            if jden == 0:
                raise ZeroDivisionError('Singular curve')
            elif jnum % jden == 0:
                return jnum//jden
            else:
                return jnum/jden
        else:
            jnum = jnum % char
            jden = jden % char
            if jden == 0:
                raise ZeroDivisionError('Singular curve')
            jdeninv = pow(jden,-1,char)
            return (jnum*jdeninv)%char


def cubic_qrs(fg:tuple[int,int],p:int)->list[int]:
    f,g = fg
    return [quad_rec((x**3+f*x+g)%p,p) for x in range(p)]

def trace_and_2tor(fg:tuple[int,int],p:int)->tuple[int,int]:
    qrs = cubic_qrs(fg,p)
    return -sum(qrs),len([s for s in qrs if s == 0])

def trace_frob(fg:tuple[int,int],p:int)->int:
    f,g = fg
    return - sum([quad_rec(x**3+f*x+g,p) for x in range(p)])

def supp(p:int)->list[int]:
    ds = []
    a = 0
    while a*a - 4*p<0:
        ds.append(a*a-4*p)
        a+=1
    return ds

def support_closed(p:int)->list[int]:
    suppbd = supp(p)
    supp_clos = []
    for d in suppbd:
        d0,c = discfac(d)
        supp_clos+=[d0*(c0**2) for c0 in divisors(c)]
    ds = list(set(supp_clos))
    ds.sort(key=abs)
    return ds

def supp_in_hilbdb(p:int)->dict:
    ds = support_closed(p)
    return {d:hilb_polys_dict[d] for d in hilb_polys_dict if d in ds}

def endo_db_check(j:int,p:int)->list[int]:
    cands = supp_in_hilbdb(p)
    return [d for d in cands if poly_eval_mod(cands[d][::-1],j,p)==0]


def fp_isog_codomains(j:int,l:int,p:int):
    y0s = atk_at_j(j,l,p).mod(p).find_roots_BrFo()
    atkin_linear = atk_poly_a(l,p).mod(p)
    return [(atkin_linear.eval(y0)-j)%p for y0 in y0s]




def get_endo_disc_cands(fg:tuple[int,int],p)->list[int]:
    j0 = fg_to_j(fg,p)
    hds = supp_in_hilbdb(p)
    db_check = [d for d in hds if
                poly_eval_mod(hds[d][::-1],j0,p)==0]
    if len(db_check)>0:
        return db_check
    f,g = fg
    if (f*g) %p != 0:
        a, t2 = trace_and_2tor(fg,p)
        d = a*a-4*p
        d0, c = discfac(d)
        if c == 1:
            return [d]
        elif c == 2:
            if t2==1:
                return [d]
            else:
                return [d0]
        cands = [d0*(c0**2) for c0 in divisors(c)]
        # If the conductor of Z[phi] is even, we can use 2torsion
        # to gain information about the order of 2 in the conductor of O_E
        if c % 2 == 0:
            cands1 = []
            cands2 = []
            for d0 in cands:
                if (d//d0)%2 == 1:
                    cands1.append(d0)
                else:
                    cands2.append(d0)
            if t2 == 1:
                cands = cands1
            else:
                cands= cands2
        #We've already ruled out some discriminants -
        #if they happen to appear in the set of candidates
        #we can discard them immediately
        cands = [d1 for d1 in cands if d1 not in hds]
        return cands
    elif f == 0:
        if p % 3 == 1:
            return [-3]
        else:
            return [-4*p]
    else:
        if p % 4 == 1:
            return [-4]
        qr = quad_rec(-f,p)
        if qr == 1:
            return [p]
        else:
            return [4*p]
    

def j_to_disc(j:int,p:int)->int:
    if j == 0:
        if p % 3 == 1:
            return -3
        else:
            return -4*p
    elif (j-1728)%p == 0:
        if p % 4 == 1:
            return - 4
        else:
            return -4*p
    else:
        f = -3*j*(j-1728)
        g = 2*j*(j-1728)**2
        a = trace_frob((f,g),p)
        return a**2 - 4*p
    
def supsingtrace(p,l):
    trace = 0
    d3seen = False
    d4seen = False
    a = 0
    while a*a < 4*l:
        d,c = discfac(a*a-4*l)
        qr = quad_rec(d,p)
        if qr < 1:
            while c % p == 0:
                c = c // p
            d0 = d*c*c
            h = len(get_cl_reps(d0))
            if d % p == 0 or d % l == 0:
                if d == -3: 
                    if not d3seen:
                        trace += h
                        d3seen = True
                    else:
                        trace+=h-1
                elif d == -4:
                    if not d4seen:
                        trace += h
                        d4seen = True
                    else:
                        trace+=h-1
                else:  
                    trace += h
            else:
                if d == -3: 
                    if not d3seen:
                        trace += 2*h
                        d3seen = True
                    else:
                        trace+=2*(h-1)
                elif d == -4:
                    if not d4seen:
                        trace += 2*h
                        d4seen = True
                    else:
                        trace+=2*(h-1)
                else:  
                    trace += 2*h
        a+=1                
    return trace

def x0l_fp_card(p,l):
    card = 0
    if quad_rec(-p,l)==1:
        card += (len(get_cl_reps(-4*p)))
    a = 1
    cond0 = []
    cond1728 = []
    while a*a < 4*p:
        d = a*a-4*p
        d0,c = discfac(d)
        if c % l == 0:
            d1 = d0*((c//l)**2)
            card1 = len(get_cl_reps(d))
            card2 = len(get_cl_reps(d1))
            if d0 == -3:
                card+=(card1-1)+l*(card2-1)
                cond0.append(c)
            elif d0 == -4:
                card+=(card1-1)+l*(card2-1)
                cond1728.append(c)
            else:
                card+=card1+l*card2
        else:
            qr = quad_rec(d0,l)
            if qr >=0:
                card+=(1+qr)*len(get_cl_reps(d))
        a+=1
    # We ignored ordinary curves j = 0, 1728
    if len(cond0) > 0:
        card+=1+quad_rec(-3,l)
        if min([c% l for c in cond0]) == 0:
            if l % 3 == 1:
                card += (l-1)//3
            else:
                card += (l+1)//3
        else:
            card+=1+quad_rec(-3,l)
    if len(cond1728) > 0:
        card+=1+quad_rec(-4,l)
        if min([c% l for c in cond1728]) == 0:
            if l % 4 == 1:
                card += (l-1)//2
            else:
                card += (l+1)//2
    return card