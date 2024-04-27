def get_sqrt_dic(p:int)->dict:
    sqrts = {a:[] for a in range(p)}
    sqrts[0] = [0]
    for r in range(1,(p+1)//2):
        s = (r*r)%p
        sqrts[s]=[-r,r]
    return sqrts

def get_points_wtwst(coefs:tuple,sqrts:dict)->list:
    p = len(sqrts)
    f,g = coefs
    d = p-1
    while pow(d,p//2,p)==1 and d>1:
        d-=1
    pts = []
    for x in range(p):
        y2 = (pow(x,3,p)+f*x+g)%p
        if y2 == 0:
            pts.append([x,0,0])
        elif len(sqrts[y2])>0:
            pts+=[[x,r,1] for r in sqrts[y2]]
        else:
            pts+=[[x,r,-1] for r in sqrts[(d*y2)%p]]
    return pts