import numpy as np
import matplotlib.pyplot as plt

## General stuff
def nearest_int(x:float)->int:
    cands = [int(x),int(x)+1,int(x)-1]
    cands.sort(key = lambda t:abs(t-x))
    return cands[0]

def tau_to_verts(tau:np.array)->np.array:
    one = np.array([1,0])
    zero = np.array([0,0])
    return np.array([zero,one,one+tau,tau,zero])

def lattice_pts(tau0:np.array,tau1:np.array,r:float)->np.array:
    l0 = np.linalg.norm(tau0)
    l1 = np.linalg.norm(tau1)
    m = int(r//min(l0,l1))
    cands = [a*tau0+b*tau1 for a in range(-m,m+1) for b in range(-m,m+1)]
    return np.array([z for z in cands if np.linalg.norm(z)<r])

def lattice_pt_grs(ap:tuple,tau_a:np.array,n):
    a,p = ap
    i = np.roots([1,0,1])[0]
    tau = np.dot(tau_a,[1,i])
    xi = np.roots([1,-a,p])[0]
    r = max([1,abs(tau),abs(tau+1)])*1.3
    lat0 = lattice_pts(np.array([1,0]),np.array([tau.real,tau.imag]),r)
    latgroups = [lat0]
    for k in range(1,n+1):
        xi_k = (xi**k)-1
        t1 = 1/xi_k
        t2 = tau/xi_k
        latgroups.append(
            lattice_pts(np.array([t1.real,t1.imag]),np.array([t2.real,t2.imag]),r)
            )
    return latgroups

def check_fd(pt,tau,ep):
    ych = (pt[1]+ep)*(tau[1]+ep-pt[1])
    if ych < 0:
        return False
    slope = tau[0]/tau[1]
    x_left_bd = pt[1]*slope
    x_right_bd = pt[1]*slope +1
    return (pt[0]-x_left_bd+ep)*(x_right_bd-pt[0]+ep)>=0

def trim_groups(pt_groups,tau,ep):
    return [np.array([pt for pt in group if check_fd(pt,tau,ep)])
                   for group in pt_groups]

def lattice_pt_grs_slim(ap,tau,n,ep):
    groups = lattice_pt_grs(ap,tau,n)
    return trim_groups(groups,tau,ep)

def to_torus(R,r,t1,t2):
    pi = np.pi
    ct1 = np.cos(2*pi*t1)
    st1 = np.sin(2*pi*t1)
    ct2 = np.cos(2*pi*t2)
    st2 = np.sin(2*pi*t2)
    x = (R+r*ct1)*ct2
    y = (R+r*ct1)*st2
    z = r*st1
    return np.array([x,y,z])

def lattice_pt_grs_3D_reconly(ap,y,n,ep):
    groups = lattice_pt_grs_slim(ap,np.array([0,y]),n,ep)
    groups3d = []
    for group in groups:
        group3d = []
        for pt in group:
            group3d.append(to_torus(1+y,1,pt[0],pt[1]/y))
        groups3d.append(np.array(group3d))
    return groups3d

def make_3d_plot(ap,y,sizes,colors):
    n = len(sizes)-1
    groups = lattice_pt_grs_3D_reconly(ap,y,n,0.001)
    ax = plt.axes(projection ="3d")
    for i,group in enumerate(groups):
        ax.scatter3D(xs = group[::,0],ys= group[::,1],zs = group[::,2],
                     s = sizes[i], c = colors[i%len(colors)])
    

## Class group reps

def red_bf(abc:tuple)->bool:
    a,b,c = abc
    if abs(b)< a and a < c:
        return True
    elif b == a and a < c:
        return True
    elif c == a and b >= 0 and b <= a:
        return True
    return False
    
def get_cl_reps(d:int)->list:
    reps_found = []
    if d % 4 > 1:
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

def get_taus(d:int)->np.array:
    cl_reps = get_cl_reps(d)
    return [np.roots(bqf)[0] for bqf in cl_reps]


class EllCurveCPVis:
    def __init__(self,frob_trace,prime):
        self.tr = frob_trace
        self.p = prime
        self.disc = frob_trace**2 - 4*prime
        self.card = prime+1-frob_trace
        self.taus = get_taus(frob_trace**2 - 4*prime)
    def __repr__(self):
        p = self.p
        c = self.card
        head = 'Visualizations of elliptic curves mod '
        mid = ' with cardinality '
        return head+str(p)+mid+str(c)
    def tau_arr(self):
        return np.array([np.array([z.real,abs(z.imag)]) for z in self.taus])
    def make_plots(self,ep=0.0,k = 0,sizes=[10,3],colors = ['black','red']):
        taus = self.tau_arr()
        tau = taus[k]
        ap = self.tr,self.p
        n = len(sizes)-1
        ptsets = lattice_pt_grs_slim(ap,tau,n,ep)
        xlist = []
        ylist = []
        slist = []
        clist = []
        n = 0
        for arr in ptsets:
            xarr = arr[::,0]
            yarr = arr[::,1]
            xlist = np.concatenate((xlist,xarr))
            ylist = np.concatenate((ylist,yarr))
            slist+=len(xarr)*[sizes[n%len(sizes)]]
            clist+=len(xarr)*[colors[n%len(colors)]]
            n+=1
        plt.figure()
        plt.scatter(x=xlist,y=ylist,s=slist,c=clist)
        plt.plot(tau_to_verts(tau)[::,0],
                tau_to_verts(tau)[::,1],
                linestyle = '--')
        plt.gca().set_aspect('equal')

    def make_plots_exp(self,ep=0.0,k = 0,sizes=[10,3],colors = ['black','red'],normalize=False):
        taus = self.tau_arr()
        tau = taus[k]
        ap = self.tr,self.p
        n = len(sizes)-1
        ptsets = lattice_pt_grs_slim(ap,tau,n,ep)
        i = np.roots([1,0,1])[0]
        xlist = []
        ylist = []
        slist = []
        clist = []
        n = 0
        for arr in ptsets:
            exparr = np.exp(2*np.pi*i*np.dot(arr,[1,i]))
            expxys = [np.array([z.real,z.imag]) for z in exparr]
            #add normalization
            expxys = np.array(expxys)
            xarr = expxys[::,0]
            yarr = expxys[::,1]
            xlist = np.concatenate((xlist,xarr))
            ylist = np.concatenate((ylist,yarr))
            slist+=len(xarr)*[sizes[n%len(sizes)]]
            clist+=len(xarr)*[colors[n%len(colors)]]
            n+=1
        plt.figure()
        plt.scatter(x=xlist,y=ylist,s=slist,c=clist)
        plt.gca().set_aspect('equal')

    def make_plots_exp_rn(self,ep=0.0,k = 0,sizes=[10,3],colors = ['black','red']):
        taus = self.tau_arr()
        tau = taus[k]
        ap = self.tr,self.p
        n = len(sizes)-1
        ptsets = lattice_pt_grs_slim(ap,tau,n,ep)
        i = np.roots([1,0,1])[0]
        xlist = []
        ylist = []
        slist = []
        clist = []
        n = 0
        for arr in ptsets:
            expl = []
            for pt in arr:
                t = 2*np.pi*pt[0]
                ct = np.cos(t)
                st = np.sin(t)
                l = (1+pt[1])
                expl.append(l*np.array([ct,st]))
            expxys = np.array(expl)
            xarr = expxys[::,0]
            yarr = expxys[::,1]
            xlist = np.concatenate((xlist,xarr))
            ylist = np.concatenate((ylist,yarr))
            slist+=len(xarr)*[sizes[n%len(sizes)]]
            clist+=len(xarr)*[colors[n%len(colors)]]
            n+=1
        plt.figure()
        plt.scatter(x=xlist,y=ylist,s=slist,c=clist)
        plt.gca().set_aspect('equal')

    

        

####
# Lattice to primes
###

from primes import primesBetween

def lattice_primes(tau:np.array,r:float,c=1):
    m = int(r)
    cands = [a*np.array([1,0])+b*tau for a in range(-m,m+1) for b in range(-m,m+1)]
    latt_small= np.array([z for z in cands if np.linalg.norm(z)<r and z[1]>=0])
    norms = [nearest_int(np.linalg.norm(x)**2) for x in latt_small]
    primes = primesBetween(0,m*m+1)
    colors = [[0,0,1,0.2] for a in latt_small]
    sizes = [3 for a in latt_small]
    primes_c = []
    for i in range(len(latt_small)):
        if norms[i] in primes and latt_small[i][1]==c*tau[1]:
            colors[i] = 'red'
            sizes[i] = 8
            primes_c.append(norms[i])
    plt.figure()
    plt.scatter(x=latt_small[::,0],y=latt_small[::,1],c=colors,s=sizes)
    ax = plt.gca()
    for p in primes_c:
        ax.add_patch(plt.Circle((0, 0), np.sqrt(p),fill = False,alpha=0.4))
    ax.set_aspect('equal')
    plt.ylim((0,r+1))
    return ax

def lattice_primesX(tau_list,clist,r:float):
    n = len(tau_list)
    m = int(r)
    fig, axs = plt.subplots(n,1)
    for j,tau in enumerate(tau_list):
        cands = [a*np.array([1,0])+b*tau for a in range(-m,m+1) for b in range(-m,m+1)]
        latt_small= np.array([z for z in cands if np.linalg.norm(z)<r and z[1]>=0])
        norms = [nearest_int(np.linalg.norm(x)**2) for x in latt_small]
        primes = primesBetween(0,m*m+1)
        colors = [[0,0,1,0.2] for a in latt_small]
        sizes = [3 for a in latt_small]
        primes_c = []
        for i in range(len(latt_small)):
            if norms[i] in primes:
                colors[i] = 'green'
                if nearest_int(latt_small[i][1]/tau[1]) in clist[j]:
                    sizes[i] = 8
                    primes_c.append(norms[i])
                    colors[i] = 'green'
                else:
                    colors[i] = [1,0,0,0.3]
        axs[j].scatter(x=latt_small[::,0],y=latt_small[::,1],c=colors,s=sizes)
        for p in primes_c:
            axs[j].add_patch(plt.Circle((0, 0), np.sqrt(p),fill = False,alpha=0.4))
        axs[j].set_aspect('equal')
        axs[j].set_ylim([0,r+1])


