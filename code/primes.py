def primesBetween(a,b):
    a = max(2,a)
    if b < a:
        return []
    elif b == 2:
        return [2]
    elif b == 3:
        if a == 2:
            return [2,3]
        else:
            return [3]
    m = b
    # The following code is going to compute the set of primes < m
    # First, we create a list that will keep track of the primes-
    # the set contains m+1 elements, and for integers k <= m,
    # cands[k] = 0 records the fact that k is composite.
    # To begin, every even element other than 2 is set to 0,
    # to record the fact that 2 is the only even prime.
    cands = [0,0,1]+[(i % 2) for i in range(3,m+1)]
    # Next, we record that multiples of 3 other than 3 are composite.
    # Note that 6 is already known to be composite, 
    # so we can start at 9.
    for i in range(9,m+1,3):
        cands[i] = 0
    # Now we start going through the list of candidates to find primes p,
    # and updating the list by recording that all multiples of the prime,
    # which lie between p^2 and m inclusive, are composite.
    # We can stop once p^2 > m, since all multiples of p^2 will exceed m.
    # We start at p = 5, and we will only check odd numbers
    # that are not multiples of 3.
    # To avoid the multiples of 3, we will jump ahead by either 2 or 4
    # at each step. The number e records whether we need to jump by 2 or 4
    # at each step.
    e = -1
    p = 5
    while p**2 <= m:
        # This assumes p is prime. This is the case when p = 5,
        # and will be the case at the end of the loop when p is redefined.
        # As mentioned above, we start by recording that all multiples of
        # p between p^2 and m are composite.
        for pm in range(p**2,m+1,p):
            cands[pm]=0
        # Now we look for the next prime.
        # We alternate adding 2 and 4 to p and then checking whether
        # p is prime by checking if cands[p] is 0 or not.
        p+=3+e
        e*=-1
        while cands[p] == 0:
            p+=3+e
            e*=-1
    # When the loop ends, cands[p] = 0 if and only if p is 0,
    # so we can obtain the set of primes:
    primes = [p for p in range(a,m+1) if cands[p] == 1]
    return primes

pfs = {}

primes = [2]

def prime_fac(n):
    if n in pfs:
        return prime_fac[n]
    pf = {}
    if n < 2:
        return {}
    if n % 2 == 0:
        e = 0
        while n % 2 == 0:
            e+=1
            n = n//2
        pf[2]=e
    if n % 3 == 0:
        e = 0
        while n % 3 == 0:
            e+=1
            n = n//3
        pf[3] = e
    d = 5
    s = -1
    while n > 1 and d*d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                e+=1
                n = n //d
            pf[d] = e
        d+=3+s
        s*=-1
    if n > 1:
        pf[n] = 1
    pfs[n] = pf
    return pf

def sqrtFind(n):
    r = 1
    while r**1 < n:
        b = 2
        while (r+b)**2 < n:
            b*=2
        r+= (b//2)
    return r


def primes_from_d(d:int,m:int):
    primes_m = primesBetween(max(4,abs(d)//4),m+1)
    primes_found = {}
    for p in primes_m:
        a2 = 4*p-abs(d)
        a = sqrtFind(a2)
        if a**2 == a2:
            primes_found[p]=a
    return primes_found