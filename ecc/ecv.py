import matplotlib as mpl
from matplotlib.patches import Polygon,RegularPolygon,Rectangle,Circle
from matplotlib import colormaps
import matplotlib.pyplot as plt
import numpy as np
from ecc.utils import discfac
from ecc.qfs import *
from ecc.ecfp import mw_gens, get_j_to_qfs_dict,trfr_to_models,frobmat
from ecc.modularpolynomials import *


def n_to_color(n):
    n = n %17
    if n < 9:
        return colormaps['Set1'](n)
    else:
        return colormaps['Dark2'](n-9)

#################
# MW Group Pics #
#################

def abc_to_tau(abc:tuple[int,int,int])->np.array:
    tau_np= np.roots(abc)[0]
    return np.array([tau_np.real,tau_np.imag])

def abc_to_tau_str(abc):
    a,b,c = qf_make_prim(abc)
    d,cond = discfac(b*b-4*a*c)
    den = 2*a
    g = gcd_list([den,b,cond])
    if g>1:
        den,b,cond = den//g, b//g, cond//g
    if d % 4 == 0 and b % 2 == 0 and den % 2 == 0:
        d = d//4
        b = b//2
        den = den//2
    num_str = f'sqrt({d})'
    if cond > 1:
        num_str = f"{cond}"+num_str
    if b!=0:
        num_str = f"{-b}+"+num_str
    if den == 1:
        return num_str
    if '+' in num_str:
        return '('+num_str+f')/{den}'
    else:
        return num_str+'/'+str(den)


def abc_to_parallelogram_verts(qf):
    tau_arr = abc_to_tau(qf)
    one_arr = np.array([1,0])
    return [0*one_arr,one_arr,one_arr+tau_arr,tau_arr]

def mw_arr_from_gens(abc:tuple,gens:dict):
    den = max(gens.values())
    one = np.array([1,0])
    tau = abc_to_tau(abc)
    pts = [np.array([0,0])]
    for gen in gens:
        x,y = gen
        pts = [(pt0+m*np.array([x,y]))%den 
               for pt0 in pts for m in range(gens[gen])]
    return np.array([pt[0]*one+pt[1]*tau for pt in pts])/den

def lattice_pic_data(ap:tuple[int,int],abc:tuple[int,int,int],k=1):
    gens = mw_gens(ap,abc,k)
    pts = mw_arr_from_gens(abc,gens)
    verts = abc_to_parallelogram_verts(abc)
    return pts, verts

def make_lattice_pic(ap:tuple[int,int],abc:tuple[int,int,int],c='black',k=1,s=10):
    pts, verts = lattice_pic_data(ap,abc,k)
    xs = pts[::,0]
    ys = pts[::,1]
    fundom = Polygon(verts)
    fig, ax = plt.subplots()
    ax.add_patch(fundom)
    ax.scatter(x=xs,y=ys,c='black',s=s)
    plt.gca().set_aspect('equal')
    plt.show()

def get_fp_points(fg,p):
    f,g = fg
    pts = []
    rtsmodp = {y:[] for y in range(p)}
    for y in range((p//2)+1):
        rtsmodp[pow(y,2,p)]=list({y,-y})
    for x in range(-(p//2),(p//2)+1):
        rhs = (x**3+f*x+g)%p
        for y in rtsmodp[rhs]:
            pts.append([x,y])
    return np.array(pts)

def make_both_pics(fg:tuple[int,int],ap:tuple[int,int],abc:tuple[int,int,int],
                   c='black',k=1,sz=10):
    pts, verts = lattice_pic_data(ap,abc,k)
    xs = pts[::,0]
    ys = pts[::,1]
    fundom = Polygon(verts,alpha=0.3,facecolor='gray')
    fig, axes = plt.subplots(nrows=1,ncols=2)
    axes[0].add_patch(fundom)
    axes[0].scatter(x=xs,y=ys,c=[c for x in pts],s=sz)
    axes[0].axis('off')
    f,g = fg
    p = ap[1]
    pts1 = get_fp_points(fg,p)
    amb = np.array([[x,y] for x in range(-(p//2),(p+1)//2)
                    for y in range(-(p//2),(p+1)//2)])
    axes[1].scatter(x = amb[::,0],y=amb[::,1],c='gray',alpha=0.3,s=0.7*sz)
    axes[1].scatter(x=pts1[::,0],y=pts1[::,1],c=[c for x in pts1],s=sz)
    axes[1].axis('off')
    plt.gca().set_aspect('equal')
    plt.title(f'y^2 = x^3 + {f}x +{g} mod {p}')
    plt.show()



##########
# X_0(l) #
##########

def modular_curve_patch(p:int,ymx:float):
    fig, ax = plt.subplots()
    cmap = mpl.colormaps['Pastel1']
    c0= np.array(cmap(1/(p+2)))
    c0[-1]=0.3
    c1 = np.array(cmap((p+1)/(p+2)))
    c1[-1]=0.3
    ax.add_patch(mpl.patches.Rectangle((-0.5,0),width=1,height=ymx,facecolor=c0))
    ax.add_patch(mpl.patches.Circle((0,0),1,facecolor=c1,edgecolor='black'))
    for m in range(1,(p//2)+1):
        if m % 2== 0:
            cm = np.array(cmap((2*m+p)/(2*p+2)))
        else:
            cm = np.array(cmap(1-(2*m+p)/(2*p+2)))
        cm[-1]=0.5
        cm2 = np.matmul([[0,1,0,0],[0,0,1,0],[1,0,0,0],[0,0,0,1]],cm)
        ax.add_patch(mpl.patches.Circle((1/(2*m-1),0),1/abs(2*m-1),facecolor=cm,edgecolor='black'))
        ax.add_patch(mpl.patches.Circle((-1/(2*m-1),0),1/abs(2*m-1),facecolor=cm2,edgecolor='black'))
    ax.add_patch(mpl.patches.Circle((1/(p),0),1/abs(p),facecolor='white',edgecolor='magenta'))
    ax.add_patch(mpl.patches.Circle((-1/(p),0),1/abs(p),facecolor='white',edgecolor='magenta'))
    for m in range(2,(p//2)+1):
        x = m/(m*m-1)
        cr = 1/(m*m-1)
        ax.add_patch(mpl.patches.Circle((x,0),cr,facecolor='white',edgecolor='black'))
        ax.add_patch(mpl.patches.Circle((-x,0),cr,facecolor='white',edgecolor='black'))
    ax.set_aspect('equal')
    ax.set_xlim([-0.5,0.5])
    ax.set_ylim([0,ymx])
    ax.axvline(x = -0.5, color = 'black',linestyle='dashed')
    ax.axvline(x = 0.5, color = 'black',linestyle='dashed')
    ax.axis('off')
    return fig, ax

def x0_isogs_data(qflist,l):
    q_to_c = {}
    for i,qf in enumerate(qflist):
        q_to_c[qf_mod_gamma(qf)]=n_to_color(i)
    iso_qs = {}
    for q0 in q_to_c:
        q0_orb = gamma_0_orb(q0,l)
        for q1 in q0_orb:
            q2 = qf_mod_gamma(fricke_inv(q1,l))
            if q2 in q_to_c:
                iso_qs[q1] = [q_to_c[q0],q_to_c[q2]]
    return iso_qs

def x0_iso_pic(qflist,l):
    colordata = x0_isogs_data(qflist,l)
    taus = np.array([abc_to_tau(qf) for qf in colordata])
    ymx = max(taus[::,1])
    fig,ax = modular_curve_patch(l,ymx+0.3)
    for qf in colordata:
        c1,c2 = colordata[qf][0],colordata[qf][1]
        tau = abc_to_tau(qf)
        x,y = tau[0],tau[1]
        ax.scatter(x = [x],y=[y],c = [c2],s=[40])
        ax.scatter(x = [x],y=[y],c = [c1],s=[8])
    return fig, ax
        


class IsogenyClassFp:
    def __init__(self,a:int,p:int):
        self.char = p
        self.trace_frob = a
        self.disc = a**2-4*p
        d,c = discfac(self.disc)
        self.fundisc = d
        self.cond = c
        self.fgs = trfr_to_models(a,p)
        self.qfs = get_qfs_all(self.disc)
        self.j_to_qf_dict = get_j_to_qfs_dict(a,p)
        self.j_list = [j for j in self.j_to_qf_dict]
        self.is_supersingular = (a%p==0)
        self.card = p-a+1
        self.qfdata = {}
        self.hgt = 1
        for j in self.j_to_qf_dict:
            qf = (self.j_to_qf_dict)[j]
            self.qfdata[qf]= {'j':j,'frob_mat':frobmat((a,p),qf),'fg':self.fgs[j],'tau_arr':abc_to_tau(qf),
            'tau_str':abc_to_tau_str(qf)}
            self.hgt = max(self.hgt,abc_to_tau(qf)[1])
        self.ht_ub = int((self.hgt+1)//1)


    def j_to_qf(self,j0:int):
        if j0 not in self.j_to_qf_dict:
            raise ValueError(f'{j0} not found')
        return (self.j_to_qf_dict)[j0]
    
    def uhfd_pic(self):
        ymax = self.ht_ub
        ymin = 0
        xmin = -0.7
        xmax = 0.7
        ud = Circle((0,0),1,alpha=0.3,facecolor='gray')
        fig, ax = plt.subplots(figsize = (4.8,self.ht_ub))
        for i,j in enumerate(self.j_list):
            qf = self.j_to_qf(j)
            cj = n_to_color(i)
            xj = self.qfdata[qf]['tau_arr'][0]
            yj = self.qfdata[qf]['tau_arr'][1]
            strj = self.qfdata[qf]['tau_str']
            ax.scatter(x = [xj],y=[yj],c=[cj],label=strj)
        ax.add_patch(ud)
        ax.vlines(x=-0.5,ymin=np.sqrt(3)/2,ymax=ymax)
        ax.vlines(x=0.5,ymin=np.sqrt(3)/2,ymax=ymax)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.legend(loc='lower right', bbox_to_anchor=(-0.6,0))
        return fig, ax
    
    def x0_pic_pair(self,l):
        ymax = self.ht_ub
        ymin = 0
        xmin = -0.7
        xmax = 0.7
        ud = Circle((0,0),1,alpha=0.3,facecolor='gray')
        fig, axs = plt.subplots(ncols=2,nrows=1,figsize = (4.8,self.ht_ub))
        q_to_c = {}
        for i,j in enumerate(self.j_list):
            qf = self.j_to_qf(j)
            cj = n_to_color(i)
            q_to_c[qf]=cj
            xj = self.qfdata[qf]['tau_arr'][0]
            yj = self.qfdata[qf]['tau_arr'][1]
            strj = self.qfdata[qf]['tau_str']
            axs[0].scatter(x = [xj],y=[yj],c=[cj],label=strj)
        axs[0].add_patch(ud)
        axs[0].vlines(x=-0.5,ymin=np.sqrt(3)/2,ymax=ymax)
        axs[0].vlines(x=0.5,ymin=np.sqrt(3)/2,ymax=ymax)
        axs[0].set_xlim(xmin,xmax)
        axs[0].set_ylim(ymin,ymax)
        axs[0].set_aspect('equal')
        axs[0].axis('off')
        axs[0].legend(loc='lower right', bbox_to_anchor=(-0.6,0.5))
        ##########
        # X0 pic #
        ##########
        isoqs = x0_isogs_data(self.qfs,l)
        isotaus = [abc_to_tau(qf) for qf in isoqs]
        if len(isotaus)==0:
            ymx2 = ymax
        else:
            ymx2 = max([xy[1] for xy in isotaus])+0.3
        cmap = mpl.colormaps['Pastel1']
        c0= np.array(cmap(1/(l+2)))
        c0[-1]=0.3
        c1 = np.array(cmap((l+1)/(l+2)))
        c1[-1]=0.3
        axs[1].add_patch(mpl.patches.Rectangle((-0.5,0),width=1,height=ymx2,facecolor=c0))
        axs[1].add_patch(mpl.patches.Circle((0,0),1,facecolor=c1,edgecolor='black'))
        for m in range(1,(l//2)+1):
            if m % 2== 0:
                cm = np.array(cmap((2*m+l)/(2*l+2)))
            else:
                cm = np.array(cmap(1-(2*m+l)/(2*l+2)))
            cm[-1]=0.5
            cm2 = np.matmul([[0,1,0,0],[0,0,1,0],[1,0,0,0],[0,0,0,1]],cm)
            axs[1].add_patch(mpl.patches.Circle((1/(2*m-1),0),1/abs(2*m-1),facecolor=cm,edgecolor='black'))
            axs[1].add_patch(mpl.patches.Circle((-1/(2*m-1),0),1/abs(2*m-1),facecolor=cm2,edgecolor='black'))
        axs[1].add_patch(mpl.patches.Circle((1/(l),0),1/abs(l),facecolor='white',edgecolor='magenta'))
        axs[1].add_patch(mpl.patches.Circle((-1/(l),0),1/abs(l),facecolor='white',edgecolor='magenta'))
        for m in range(2,(l//2)+1):
            x = m/(m*m-1)
            cr = 1/(m*m-1)
            axs[1].add_patch(mpl.patches.Circle((x,0),cr,facecolor='white',edgecolor='black'))
            axs[1].add_patch(mpl.patches.Circle((-x,0),cr,facecolor='white',edgecolor='black'))
        axs[1].set_aspect('equal')
        axs[1].set_xlim([-0.6,0.6])
        axs[1].set_ylim([0,ymx2])
        axs[1].axvline(x = -0.5, color = 'black',linestyle='dashed')
        axs[1].axvline(x = 0.5, color = 'black',linestyle='dashed')
        axs[1].axis('off')
        if len(isoqs)>0:
            for qf in isoqs:
                c1 = q_to_c[qf_mod_gamma(qf)]
                c2 = q_to_c[qf_mod_gamma(fricke_inv(qf,l))]
                tau = abc_to_tau(qf)
                x,y = tau[0],tau[1]
                axs[1].scatter(x = [x],y=[y],c = [c2],s=[40])
                axs[1].scatter(x = [x],y=[y],c = [c1],s=[8])
        return fig, axs

    def qfs_with_mw_groups(self,k=1):
        if k < 1:
            raise ValueError('Field exponent should be at least 1')
        qfdata = (self.qfdata).copy()
        ap = self.trace_frob,self.char
        for qf in qfdata:
            qfdata[qf]['mw_gens']=mw_gens(ap,qf,k)
        return qfdata
    
    def make_pics_single(self,j0:int,k=1,sz=20,clr='black'):
        qf = self.j_to_qf(j0)
        fg = (self.fgs)[j0]
        ap = self.trace_frob, self.char
        pts, verts = lattice_pic_data(ap,qf,k)
        xs = pts[::,0]
        ys = pts[::,1]
        fundom = Polygon(verts,alpha=0.3,facecolor='gray')
        fig, axes = plt.subplots(nrows=1,ncols=2)
        axes[0].add_patch(fundom)
        axes[0].scatter(x=xs,y=ys,c=[clr for x in pts],s=sz)
        axes[0].axis('off')
        axes[0].set_title('tau = '+abc_to_tau_str(qf))
        axes[1].set_aspect('equal')
        f,g = fg
        p = ap[1]
        pts1 = get_fp_points(fg,p)
        amb = np.array([[x,y] for x in range(-(p//2),(p+1)//2)
                        for y in range(-(p//2),(p+1)//2)])
        axes[1].scatter(x = amb[::,0],y=amb[::,1],c='gray',alpha=0.3,s=0.7*sz)
        axes[1].scatter(x=pts1[::,0],y=pts1[::,1],c=[clr for x in pts1],s=sz)
        axes[1].axis('off')
        axes[1].set_aspect('equal')
        axes[1].set_title(f'y^2 = x^3 + {f}x +{g} mod {p}')
        return fig, axes
        
        
    def make_pics_all(self,k=1,sz=20):
        a,p = self.trace_frob, self.char
        ap = a,p
        data = self.qfs_with_mw_groups(k)
        amb = np.array([[x,y] for x in range(-(p//2),(p+1)//2)
                    for y in range(-(p//2),(p+1)//2)])
        qfs = [qf for qf in data]
        npics = len(qfs)
        if npics == 1:
            j0 = data[qfs[0]]['j']
            d = a**2-4*p
            hd = hilb_polys_dict[d]
            jz = -hd[0]
            fig, axes = self.make_pics_single(data[qfs[0]]['j'])
            fig.suptitle(f'The elliptic curve over F_{p} \n with j = {jz} = {j0} mod {p} \n and trace of Frobenius {a}')
            return fig, axes
        else:
            fig, axs = plt.subplots(2,npics,figsize=(4*npics,4*2))
            for i,qf in enumerate(qfs):
                ci = n_to_color(i)
                fg = data[qf]['fg']
                f,g = fg
                pts0 = get_fp_points(fg,p)
                axs[0,i].scatter(x = amb[::,0],y=amb[::,1],c='gray',alpha=0.3,s=0.7*sz)
                axs[0,i].scatter(x=pts0[::,0],y=pts0[::,1],
                                c=[ci for _ in pts0],s=sz)
                axs[0,i].set_title(f'y^2 = x^3 + {f}x +{g} mod {p}')
                axs[0,i].axis('off')
                axs[0,i].set_aspect('equal')
                pts, verts = lattice_pic_data(ap,qf,k)
                fundom = Polygon(verts,alpha=0.3,facecolor='gray')
                axs[1,i].add_patch(fundom)
                axs[1,i].scatter(x=pts[::,0],y=pts[::,1],c=[ci for _ in pts],s=sz)
                axs[1,i].axis('off')
                axs[1,i].set_title('tau = '+abc_to_tau_str(qf))
                axs[1,i].set_aspect('equal')
            fig.suptitle(f'Elliptic curves mod {p} with trace of Frobenius {a}')
            return fig, axs
