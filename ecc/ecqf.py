import itertools
from ecc.qfs import *
from ecc.nt import discfac,quad_rec,divisors,primesBetween
from ecc.identities import *
from ecc.ecfp import fp_isog_codomains, trfr_to_js
from ecc.modularpolynomials import *

ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
heeg_dics = {}
for d in heeg_js:
    d0 = discfac(d)[0]
    if d0 == d:
        heeg_dics[d]={heeg_js[d]:class_group_id(d)}
    else:
        heeg_dics[d] = {heeg_js[d0]:class_group_id(d0),heeg_js[d]:class_group_id(d)}

def small_bij_check(d):
    if d in heeg_dics:
        return heeg_dics[d]
    else:
        return {}
#########################
# Load precomputed data #
#########################
import json 

def strtup_to_tup(s):
    return tuple([int(s0) for s0 in s[1:-1].split(',')])
with open('ecc/data/ecqf_ord_pcbij_4to256.json', 'r') as f:
    ecqf_ord_pcbij_4to256_loaded = json.load(f)
ecqf_ord_256_pc = {strtup_to_tup(aps):{int(ns):tuple(ecqf_ord_pcbij_4to256_loaded[aps][ns] )
                                       for ns in ecqf_ord_pcbij_4to256_loaded[aps]} 
                                       for aps in ecqf_ord_pcbij_4to256_loaded}

###########################
# Obtaining generating ls #
###########################

def scoredic_to_tups(scores:dict):
    ls = [l for l in scores]
    ls.sort(reverse=True)
    return tuple(ls),tuple([scores[l] for l in ls])

def qf_ev(qf,m):
    a,b,c = qf
    return list({a*x*x+b*x*y+c*y*y for x in range(-m,m+1) for y in range(-m,m+1)})


def qf_reps_pm(d):
    return [qf for qf in get_qfs_strict(d) if qf[1]>=0]


def disc_to_ssls(d):
    ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
    qfsfound = {l:[] for l in ssprimes}
    qf_reps = qf_reps_pm(d)
    ls = {}
    qf0 = class_group_id(d)
    for qf in qf_reps:
        qfls = [l for l in qf_ev(qf,10) if l in ssprimes]
        if len(qfls)>0:
            l0 = min(qfls)
            ls[l0]=len(qf_isog_cycle(qf0,l0))
    return ls

def disc_to_ssl_qfs(d):
    ssprimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
    qfsfound = {l:[] for l in ssprimes}
    qf_reps = qf_reps_pm(d)
    ls = {}
    qf0 = class_group_id(d)
    for qf in qf_reps:
        qfls = [l for l in qf_ev(qf,10) if l in ssprimes]
        if len(qfls)>0:
            for l in qfls:
                ls[qf] = min(qfls)
    return ls

def gen_qfs_d_ls(d,lset):
    qf0 = class_group_id(d)
    qfs = [qf0]
    for l in lset:
        qfs_new = []
        for qf in qfs:
            qfs_new+=qf_isog_cycle(qf,l)
        qfs = list(set(qfs_new))
    return qfs

def score_tuple(ls,scores):
    prod = 1
    for l in ls:
        prod*=scores[l]
    return prod

def check_lset(d,ls):
    cld = clgr_size_gen(d)
    lgr_size = len(gen_qfs_d_ls(d,ls))
    prd = 1
    nlist = []
    for l in ls:
        nl = len(gen_qfs_d_ls(d,[l]))
        prd*=nl
    return [cld,lgr_size,prd]

def minimize_lset(d,ls):
    ls = [l for l in ls]
    n1 = len(gen_qfs_d_ls(d,ls))
    l0s= [ls[:i]+ls[i+1:] for i in range(len(ls))]
    l0s_g =[l0 for l0 in l0s if len(gen_qfs_d_ls(d,l0))==n1]
    while len(l0s_g)>0:
        l0s_g.sort(key = max)
        ls = l0s_g[0]
        l0s= [ls[:i]+ls[i+1:] for i in range(len(ls))]
        l0s_g =[l0 for l0 in l0s if len(gen_qfs_d_ls(d,l0))==n1]
    return tuple(ls)

def get_spanl2s(d):
    l2s = [ln[0] for ln in disc_to_ssls(d).items() if ln[1]==2]
    mx = min(len(l2s),clgr_2len(d)-1)
    if mx < 2:
        return [(l,) for l in l2s]
    all_ls = {}
    for k in range(2,mx+1):
        subs = []
        for lks in itertools.combinations(l2s,k):
            if len(gen_qfs_d_ls(d,list(lks)))==2**k:
                subs.append(lks)
        all_ls[k]=subs
    ks_nz = [k for k in all_ls if len(all_ls[k])>0]
    k = max(ks_nz)
    return all_ls[k]

def clgr_n2k_ssls(d):
    cld = clgr_size_gen(d)
    if cld == 1:
        return [()]
    ldata = disc_to_ssls(d)
    if len(ldata)==0:
        raise ValueError('No data found')
    ml = max(ldata.values())
    lms = [l for l in ldata if ldata[l]==ml]
    if ml == cld:
        return [(l,) for l in lms]
    ls2s = get_spanl2s(d)
    if ml == 2:
        return ls2s
    ln2s = []
    for lm in lms:
        for ls2 in ls2s:
            ln2 = [lm]+list(ls2)
            if len(gen_qfs_d_ls(d,ln2))==cld:
                return [minimize_lset(d,ln2)]
    return []


def disc_search_ssls_by_size(d,k=None):
    cld = clgr_size_gen(d)
    ssls = disc_to_ssls(d)
    ls_all = [l for l in ssls]
    if k == None:
        k = len(minimize_lset(d,ls_all))
    ls_gen = {}
    for ls in itertools.combinations(ls_all,k):
        if len(gen_qfs_d_ls(d,list(ls))) == cld:
            ls_gen[ls]=score_tuple(ls,ssls)
    if len(ls_gen)== 0:
        return []
    m = min(ls_gen.values())
    return [{l:ssls[l] for l in ls} for ls in ls_gen if ls_gen[ls]==m]

def disc_search_ssl_gen(d):
    n2k_ls = clgr_n2k_ssls(d)
    if len(n2k_ls)>0:
        return n2k_ls
    else:
        scores = disc_search_ssls_by_size(d)
        return disc_search_ssls_by_size(d)

def qf_trios_from_ltrio(d,l123):
    qf0 = class_group_id(d)
    l1,l2,l3 = l123
    qfs_l1 = qf_isogs_hor(qf0,l1)
    qfs_l2 = qf_isogs_hor(qf0,l2)
    qfs_l3 = qf_isogs_hor(qf0,l3)
    qftrios = []
    for qf1 in qfs_l1:
        qf1_l2_qfs = qf_isogs_hor(qf1,l2)
        for qf2 in qfs_l2:
            qfs_l123 = [qf for qf in qf_isogs_hor(qf2,l1) if qf in qf1_l2_qfs and qf in qfs_l3]
            if len(qfs_l123)>0:
                qftrios.append([qf1,qf2,qfs_l123[0]])
    return qftrios

def qf_nn_rig_basis(d,ldata=None):
    if ldata == None:
        ldata = disc_rig_ssl_search(d)
    if len(ldata['ls'])!=2 or ldata['l_sum'] == None:
        raise ValueError('Use diff algorithm')
    l1,l2 = ldata['ls']
    l3 = ldata['l_sum']
    trios = qf_trios_from_ltrio(d,(l1,l2,l3))
    if len(trios)== 0:
        raise ValueError('No basis found')
    qf1,qf2,qf3 = tuple(trios[0])
    return {():class_group_id(d),(l1,):qf1,(l2,):qf2,(l1,l2):qf3,l3:qf3}


def qf_rig_ssl_data(d,lscores):
    ls_big= [l for l in lscores if lscores[l]>2]
    ls,ns = scoredic_to_tups(lscores)
    data = {'ls':ls,'ns':ns,'needs_sum':False,'l_sum':None}
    if len(ls_big)<2:
        return data
    data['needs_sum']=True
    qf0 = class_group_id(d)
    qfs = [qf0]
    for l in ls_big:
        qfs_ext = []
        for qf in qfs:
            qfs_ext+=qf_isogs_hor(qf,l)
        qfs = qfs_ext
    qf_to_l_data = disc_to_ssl_qfs(d)
    for qf in qfs:
        if qf in qf_to_l_data:
            data['l_sum'] = qf_to_l_data[qf]
            if len(ls_big)==2:
                rigbasis = qf_nn_rig_basis(d,data)
                data['nn_frame'] = rigbasis
            return data
    return data



def disc_rig_ssl_search(d):
    cld = clgr_size_gen(d)
    if cld == 1:
        return {'ls':(),'ns':(),'needs_sum':False,'l_sum':None}
    ssls = disc_to_ssls(d)
    ls_n2k = clgr_n2k_ssls(d)
    if len(ls_n2k)>0:
        ls = ls_n2k[0]
        n0 = max({ssls[l] for l in ls})
        ns = tuple([n0]+[2 for _ in range(len(ls)-1)])
        return {'ls':ls,'ns':ns,'needs_sum':False,'l_sum':None}
    ls_all = [l for l in ssls]
    k = len(minimize_lset(d,ls_all))
    ls_gen = {}
    for ls in itertools.combinations(ls_all,k):
        if len(gen_qfs_d_ls(d,list(ls))) == cld:
            ls_gen[ls]=score_tuple(ls,ssls)
    if len(ls_gen)== 0:
        return []
    m = min(ls_gen.values())
    if m != cld:
        return [ls for ls in ls_gen if ls_gen[ls]==m]
    ls_gen_l = [l for l in ls_gen]
    for ls in ls_gen_l[:-1]:
        rigdata = qf_rig_ssl_data(d,{l:ssls[l] for l in ls})
        if (rigdata['l_sum']!=None) or (not rigdata['needs_sum']):
            return rigdata
    return qf_rig_ssl_data(d,{l:ssls[l] for l in ls_gen_l[-1]})



###################################
# Reconstructing the Cayley graph #
###################################

def qf_isopwr_ngbrs(qf0,lk):
    cyc = qf_isog_cycle_power(qf0,lk)
    if len(cyc)<2:
        return []
    elif len(cyc)==2:
        return [cyc[1]]
    else:
        return [cyc[k] for k in [1,-1]]

def qf_isopwr_intrs(qfl0,qfl1):
    qf0,l0 = qfl0
    qf1,l1 = qfl1
    return [qf2 for qf2 in qf_isopwr_ngbrs(qf0,l0) if qf2 in qf_isopwr_ngbrs(qf1,l1)]


def fp_isog_cycle(jp,l):
    j,p = jp
    cyc = [j]
    nbrs = fp_isog_codomains(j,l,p)
    if len(nbrs)==0:
        return cyc
    cyc.append(nbrs[0])
    if len(nbrs)==1:
        return cyc
    nextbatch = [j2 for j2 in fp_isog_codomains(cyc[-1],l,p) if j2 not in cyc]
    while len(nextbatch)>0:
        cyc.append(nextbatch[0])
        nextbatch = [j2 for j2 in fp_isog_codomains(cyc[-1],l,p) if j2 not in cyc]
    return cyc

def fp_isog_cycle_power(jp,lk):
    l,k = lk
    cycl = fp_isog_cycle(jp,l)
    n = len(cycl)
    kord = len(set((k*a)%n for a in range(n)))
    return [cycl[(k*a)%n] for a in range(kord)]

def fp_isopwr_ngbrs(jp,lk):
    cyc = fp_isog_cycle_power(jp,lk)
    if len(cyc)<2:
        return []
    elif len(cyc)==2:
        return [cyc[1]]
    else:
        return [cyc[k] for k in [1,-1]]

def fp_isopwr_intrs(jpl0,jpl1):
    jp0,l0 = jpl0
    jp1,l1 = jpl1
    return [j for j in fp_isopwr_ngbrs(jp0,l0) if j in fp_isopwr_ngbrs(jp1,l1)]



#######################
# End to end bijection #
########################

# Currently only works for n2k groups #

def tree_edges_to_ancestors(nbrdata):
    leaves = [v for v in nbrdata if len(nbrdata[v])==1 and v != 0]
    anc_data = {}
    nextbatch = []
    for v0 in leaves:
        v1 = nbrdata[v0][0]
        anc_data[v0]=v1
        nextbatch.append(v1)
    while len(nextbatch)>0:
        currentbatch = nextbatch
        nextbatch = []
        for v in currentbatch:
            v_ancs = [v1 for v1 in nbrdata[v] if v1 not in anc_data]
            if len(v_ancs)==1:
                anc_data[v]=v_ancs[0]
                nextbatch+=v_ancs
    return anc_data



def js_to_rabs(js,p):
    ab_to_js = {}
    for i,j1 in enumerate(js):
        for j2 in js[i:]:
            ab_to_js[((-j1-j2)%p,(j1*j2)%p)] = (j1,j2)
    return ab_to_js

def get_nbr_data(ap,l):
    a,p = ap
    js = trfr_to_js(a,p)
    rabs = js_to_rabs(js,p)
    nbrdata = {j:[] for j in js}
    for x in range(p):
        evx = eval_atk(x,l,p)
        if evx in rabs:
            j1,j2 = rabs[evx]
            nbrdata[j1].append(j2)
            nbrdata[j2].append(j1)
    return nbrdata


def get_ancestor_data_ord(ap):
    a,p = ap
    d,c = discfac(a*a-4*p)
    js = trfr_to_js(a,p)
    if c == 1:
        return {'ancestor_data':{},'leaves':js,'js_all':js}
    elif c == 2:
        if d == -4:
            return {'ancestor_data':{287496%p:1728%p},'leaves':[287496%p],'js_all':js}
        leaves = [j for j in js if quad_rec(j-1728,p)==-1]
        ancs = {}
        nbrs = get_nbr_data(ap,2)
        for j in leaves:
            ancs[j]=nbrs[j][0]
        return {'ancestor_data':{2:ancs},'leaves':leaves,'js_all':js}
    ls = [l for l in primesBetween(1,c+1) if c%l ==0]
    anc_data = {}
    leaf_cands = [j for j in js if j!= 0 and (j-1728)%p !=0]
    for l in ls:
        nbrs_l = get_nbr_data(ap,l)
        anc_data[l]=tree_edges_to_ancestors(nbrs_l)
        leaf_cands = [j for j in leaf_cands if len(nbrs_l[j])==1]
    return {'ancestor_data':anc_data,'leaves':leaf_cands,'js_all':js}

def compute_ecqf_bij_ord_n2_leaves(ap,leaves=None,lcands = ssprimes):
    a,p = ap
    if a == 0 or a**2 > 4*p:
        raise ValueError('use different trace')
    d = a*a-4*p
    if leaves == None:
        vertical_iso_data = get_ancestor_data_ord(ap)
        leaves = [j for j in vertical_iso_data['leaves'] if j*(j-1728)%p!=0]
    if len(leaves)<4:
        # Any bijection works 
        # Will return a random one to avoid weird issues with j = 0,1728
        qfs = get_qfs_strict(d)
        assert len(qfs)== len(leaves)
        return {ji:qfs[i] for i,ji in enumerate(leaves)}
    j0 = leaves[0]
    # Look up degrees that are guaranteed to give connected graph
    ldata = disc_rig_ssl_search(d)
    if ldata['needs_sum']:
        raise ValueError('Use different algorithm')
    j_to_qf = {j0:class_group_id(d)}
    qf_to_j = {class_group_id(d):j0}
    ls = ldata['ls']
    for lk in ls:
        newassignments = {}
        for qf in qf_to_j:
            j = qf_to_j[qf]
            jcycle = fp_isog_cycle((j,p),lk)
            qfcycle = qf_isog_cycle(qf,lk)
            assert len(jcycle) == len(qfcycle)
            for j1,qf1 in zip(jcycle,qfcycle):
                if j1 in j_to_qf:
                    assert j_to_qf[j1] == qf1
                else:
                    newassignments[j1] = qf1
        for j1,qf1 in newassignments.items():
            j_to_qf[j1] = qf1
            qf_to_j[qf1] = j1
    return j_to_qf

    
def vert_isog_ext(j_to_qf:dict,vertical_iso_data:dict)->dict:
    for l in vertical_iso_data['ancestor_data']:
        ancsl = vertical_iso_data['ancestor_data'][l]
        nextbatch = [j for j in j_to_qf if j in ancsl and ancsl[j] not in j_to_qf]
        while len(nextbatch)>0:
            currentbatch = nextbatch.copy()
            nextbatch = []
            for j0 in currentbatch:
                j1 = ancsl[j0]
                if j1 not in j_to_qf:
                    qf0 = j_to_qf[j0]
                    qf1 = qf_parents(qf0,l)[0]
                    j_to_qf[j1] = qf1
                    if j1 in ancsl:
                        nextbatch.append(j1)
    return j_to_qf

def compute_ecqf_bij_n2_ord(ap):
    a,p = ap
    assert a**2 < 4*p
    if a == 0:
        raise ValueError('Use diff algo')
    d = a**2 - 4*p
    #check if d is Heegner number
    hd = small_bij_check(d)
    if len(hd)>0:
        data = {}
        return {j%p:hd[j] for j in hd}
    #If we have conductor 1, we don't need 2 steps - just return the output of leaf alg
    if discfac(a**2-4*p)[1] == 1:
        return compute_ecqf_bij_ord_n2_leaves(ap)
    #We're now in the general situation. We compute the vertical isogenies
    #to obtain the leaves.
    vert_data = get_ancestor_data_ord(ap)
    leaves = vert_data['leaves']
    # Compute the bijection between leaves
    j_to_qf_leaves = compute_ecqf_bij_ord_n2_leaves(ap,leaves)
    # Extend to the rest of the category using the data again
    return vert_isog_ext(j_to_qf_leaves,vert_data)
