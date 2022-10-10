# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:43:16 2019

@author: Michele Boreale

################### USAGE ##################################
1) itersdm(p,F,rho,vars,n,dmon=transpi)
INPUT:
    p = polynomial expression in the variables in vars
    F = list of polynomials expressions defining the rhs of variables in vars
    rho = list of reals, initial conditions
    vars = list of variables
    n = nonnegative integer, number of stream elements to be generated
    dmon = transition function for product pi. Can be user defined, predefined functions are
        * dmLie: shuffle   
        * dmCon: convolution   
        * dmHad: Hadamard
        * dmInf: infiltration
        * dmCon2: convolution, alternate definition
RETURNS first n elements of the stream defined by p and F

2) post(pt,F,P,vars0,vars,pl,dmon=transpi)
INPUT:
    F,vars,dmon: as above
    pt: polynomial template (see below)
    pl: list of template parameters in pt (see below)
    P: list of polynomials in the variables in vars0, define initial conditions
    vars0: variables representing the (generic) initial value of functions in vars
RETURNS a pair (qt,G) where
    * qt is a template, representing all instantiations of pt that lie in the morphism kernel (are mapped to the 0 stream)
    * G is a Groebner basis for the smallest ideal that: a) is an invariant (closed w.r.t. transition), 
      b) is contained in the  morphism kernel, and c) contains those instances.
    
3) pt,pl=genpt(vars,k)
INPUT:
    vars: list of variables
    k: nonegative integer
RETURNS a pair (pt,pl) where
    * pt is a template, that is a formal sum of all monomials over vars, up to degree k, with undetermined coefficients (parameters)
    * pl is a list of the parameters occurring in pt
    
    
#################### EXAMPLES ##############################
var('pos,f,g,f0,g0,x0,pos0')  # define a few symbols

# Fibonacci
xfib=[x,f,g]
Fib=[1,g,f+g]
rhoFib=[0,0,1]
itersdm(f,Fib,rhoFib,xfib,10,dmon=dmCon)  # generates first 11 elements of Fibonacci sequence f; can use any product here
pt,pl=genpt([x,f],3)
resFib=post(pt,Fib,[x0,f0,g0-1],[x0,f0,g0],xfib,pl,dmCon) # valid polynomial identities with Convolution product, up to degree 3 involving x and f (=Fibonacci)
# a generating function can be obtained from those equations: f=-x/(x**2 + x - 1)


# Catalan numbers
xcat=[x,f]
Cat=[1,f**2]
rhoCat=[0,1]
itersdm(f,Cat,rhoCat,xcat,10,dmon=dmCon)  # generates first 11 elements of Catalan sequence f
pt,pl=genpt([x,f],3)
resCat=post(pt,Cat,[x0,f0-1],[x0,f0],xcat,pl,dmCon) # valid polynomial identities up to degree 3 involving x and f (=Catalan numbers)
# a generating function can be obtained from those equations: f= (1-sqrt(-4*x + 1))/(2*x)

# Factorial numbers
# Same equations, but with Shuffle product (Lie derivative)  
itersdm(f,Cat,rhoCat,xcat,10,dmon=dmLie)  # generates first 11 elements of Factorial sequence 
resFact=post(pt,Cat,[x0,f0-1],[x0,f0],xcat,pl,dmLie) # valid polynomial identities up to degree 3 involving x and f (=Factorial numbers)
# an exponential generating function can be obtained from those equations: f=1/(1-x) 

# Same equations, but with Infiltration product
itersdm(f,Cat,rhoCat,xcat,5,dmon=dmInf) 
# [1, 1, 3, 29, 1667, 3254781]
# A255597, upper bound on the number of different Euler diagrams for n classes (OIES). 


# Factorial, with Hadamard
xfactH=[x,pos,f]
FactH=[1,pos+1,f*pos]
rhoFactH=[0,1,1]
itersdm(f,FactH,rhoFactH,xfactH,10,dmon=dmHad)
pt,pl=genpt([x,pos,f],3)
resFactH=post(pt,FactH,[x0,pos0-1,f0-1],[x0,pos0,f0],xfactH,pl,dmHad)

# Same equations as above, but with different products. We get, according to OEIS, the following sequences:
itersdm(f,FactH,rhoFactH,xfactH,10,dmon=dmCon) # A001333, numerators of continued fraction convergents to sqrt(2)
itersdm(f,FactH,rhoFactH,xfactH,10,dmon=dmLie) # A217924, row sums of triangle A217537
itersdm(f,FactH,rhoFactH,xfactH,7, dmon=dmInf) # A322661, number of graphs with loops spanning n labeled vertices. 
pt,pl=genpt([x,f],3)
resFactC=post(pt,FactH,[x0,pos0-1,f0-1],[x0,pos0,f0],xfactH,pl,dmCon) 
# gives the equation f*x**2 + f*x - f + x = 0, hence the g.f. f=(x - 1)/(x**2 + 2*x - 1)
# no interesting equation for Shuffle product.

# Double factorial of odd numbers
xdb=[x,f]
Fdb=[1,f**3]
itersdm(f,Fdb,[0,1],xdb,10,dmon=dmLie) # A001147, double factorial of odd numbers
pt,pl=genpt([x,f],3)
resDF=post(pt,Fdb,[x0,f0-1],[x0,f0],xdb,pl,dmLie)
# gives the equation f**2*x - f**2/2 + 1/2 = 0, solving which we obtain the exponential g.f. f=sqrt(-1/(2*x - 1)) for A001147

"""

from sympy import *
import time
 
x, y, z, w, t , u, v, r= symbols('x y z w t u v r')
v1,v2, k, x1, x2, vx, vy= symbols('v1 v2 k x1 x2 vx vy')
a, b, c, d, e, f, g, h = symbols('a b c d e f g h')



#####################



maxiter = 200

def genpt(Xlist,deg):
    monlist=list(itermonomials(Xlist,deg))
    l=len(monlist)
    parlist = list(var('a%d' % j) for j in range(l))
    prod=((Matrix(monlist)).T*Matrix(parlist))[0,0]
    return prod, parlist

def linind(list,p):     # stub: should check linear independence of p from res
    return not(p==0)

def instantiate(ptlist,zeropar,var):
    res = []
    for pt in ptlist:
        if zeropar.keys()==[]:
            if linind(res,pt):
                res.append(pt)
        else:
            for a in zeropar.keys():
                tmp=zeropar.copy()
                tmp.update({a:1})
                p=Poly(pt.subs(tmp),var,domain='QQ')
                if linind(res,p):
                    res.append(p)
    return res # returns list of polynomials in R

def check(base,newinst,var):
    G = groebner(base,var,order='grevlex',domain='QQ')
    for p in newinst:
        if p!=0:
            _,r=reduced(p,G)
            if r!=0:
                return False,G
    return True,G


step = 10

expandflag = true

def seqsub(sigma,pt):   # sequentially applies a list og substitutions in sigma to template pt
    for s in sigma:
        pt=pt.subs(s)
    return pt

def expansion(p):
    if expandflag:
        return expand(p)
    else:
        return p



def post(pt,F,P,rho,Xlist,Plist,dmon):
    start_time = time.time()
    sigma0={Xlist[i]:rho[i] for i in range(len(Xlist))}
    if P==[]:
        dummy=var('dummy')
        extra=[dummy]
        P=[dummy]       # avoids empty set of generators; Python 3 bug
    else:
        extra=[]
    G0=groebner(P,Plist+rho+Xlist+extra,order='lex',domain='QQ')
    print('G0=',G0)
    #pt=reduced(pt,G0)[1]/1
    A,rt=reduced(Poly(pt.subs(sigma0),Plist+rho+Xlist+extra,domain='QQ'),G0)
    print("pi(0)= ",pt)
    print("r_0  = ",rt)
    coeffs = Poly(rt,rho).coeffs() # list of linear expressions (constraints) extracted as coefficients of rt    
    sigma=[]    # list of individual substitutions for parameters 
    print('coeffs=',coeffs)    
    for lin in coeffs:      # solve the list of linear constraints, one by one
        lin=seqsub(sigma,lin)  # apply previously accumulated subst. to current constraint lin
        if Add(lin,0).free_symbols==set({}):
            if lin==0:
                C=True
            else:
                sigma=sigma+[{lin:0}]
                print("Linear constraint for V_"+str(j+1)+": "+str(sigma))
                print("Cannot find instance of given template that is an invariant.")
                print('m='+str(j+1))
                print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                return None,False
        else:
            C=solve(lin,dict=True)      # solve one linear constraint lin
            if C==[]:
                print('lin=',lin)
                print('No solution found')
                print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                return None,False
            s = C[0]                          # obtain subst. s from solution
            sigma=sigma+[s]                   # append s to list of substitutions  
    
    print("Linear constraint for V_0: "+str(sigma))        

    pt= seqsub(sigma,pt)                                     # apply sigma to pt
    freepar =  Add(pt/1,0).free_symbols -set(rho+Xlist) 
    zeropar = {a:0 for a in freepar}
    derlist=[pt]
    base = []
    updbase = False
        
    for j in range(maxiter):
        print(j)
        newpt=sdm(pt,F,rho,Xlist,dmon)/1#reduced(sdm(pt,F,rho,Xlist,dmon),G0)[1]/1#sdm(pt,F,rho,Xlist,dmon)##expansion((Matrix(pt.gradient(Xlist))*VF)[0,0])   # compute Lie-derivative of pt
        #newpt=reduced(newpt,G0)[1]/1        
        print("pi("+str(j+1)+")= ",newpt)
        rr,rt=reduced(Poly(pt.subs(sigma0),Plist+rho+Xlist+extra,domain='QQ'),G0)#RU(JR.reduce(RU(newpt)))            # reduce newpt modulo JR
        print("r_"+str(j+1)+"= ",rt)
        coeffs = Poly(rt,rho).coeffs() # list of linear expressions extracted as coefficients of rt    
        sigma=[]                                                
        for lin in coeffs:      # solve the list of linear constraints, one by one
            lin=seqsub(sigma,lin)  # apply previously accumulated subst. to current constraint lin
            if Add(lin,0).free_symbols==set({}):
                if lin==0:
                    C=True
                else:
                    sigma=sigma+[{lin:0}]
                    print("Linear constraint for V_"+str(j+1)+": "+str(sigma))
                    print("Cannot find instance of given template that is an invariant.")
                    print('m='+str(j+1))
                    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                    return None,False
            else:
                C=solve(lin,dict=True)      # solve one linear constraint lin
                if C==[]:
                    print('lin=',lin)
                    print('No solution found')
                    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))
                    return None,False
                s = C[0]                          # obtain subst. s from solution
                sigma=sigma+[s]                   # append s to list of substitutions  
        
        print("Linear constraint for V_"+str(j+1)+": "+str(sigma))

        if sigma==[]:
            print("Vector space equality detected: V_"+str(j)+" = V_"+str(j+1))
            print("Checking ideal equality J_"+str(j)+" = J_"+str(j+1)+"...")
            if not(updbase):
                zeropar = { a:0 for a in freepar }
                base = instantiate(derlist,zeropar,rho+Xlist)   
                updbase = True
            newinst=instantiate([newpt],zeropar,rho+Xlist)
            #print('newinst=',newinst)
            #print('base=',base)
            flag,Jm=check(base,newinst,rho+Xlist)
            if flag:
                print("Equality holds, both chains stabilized")
                print('m='+str(j))
                print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))                
                return (derlist[0],Jm)
            else:
                print("Equality does not hold, chains not yet stabilized")
                derlist.append(newpt)
                base = base + newinst
                pt=newpt
        else:
            updbase=False
            derlist = [  seqsub(sigma,ptt)  for ptt in derlist]
            pt = seqsub(sigma,newpt)
            if not(pt== (0)):
                derlist.append(pt)
            print('s=',s)
            freepar = freepar-set(s.keys()) 
    print('m='+str(j))
    print("--- Elapsed time: %s seconds ---" % (time.time() - start_time))    
    return (derlist[0],base,'WARNING: maximum number of iterations reached, result is likely to be not correct.')
 
    
    
    
############### stream derivatives
def conv(f,g,n):
    return sum([f**i*g**(n-i) for i in range(n+1)])


def sdm0(m,F,rho,Xlist):
    if m==():
        return 0
    if m[0]==0:
        sdm0(m[1:],F[1:],rho[1:],Xlist[1:])
    #print(m1,Xlist[1:])
    if Xlist[1:]==[]:
        tau=1
    else:
        tau=Poly({m[1:]:1},Xlist[1:])/1
    #print('tau=',tau)
    A=F[0]*tau*conv(Xlist[0],rho[0],m[0]-1)
    B=(rho[0]**m[0])*sdm0(m[1:],F[1:],rho[1:],Xlist[1:])
    return Poly(A+B,Xlist)/1
    
def sdm0H(m,F,rho,Xlist):
    return Poly(prod([F[i]**m[i] for i in range(len(Xlist))]),Xlist)/1
    
    
def sdm(p,F,rho,Xlist,dmon):
    if dmon=='lie':
        return lie(p,F,Xlist)
    if type(p)!=type(Poly(0,Xlist)):
        p=Poly(p,Xlist)
    d=p.as_dict()
    s=Poly(sum([dmon(m,F,rho,Xlist)*d[m] for m in d.keys()]),Xlist)
    return s
    
def itersdm(p,F,rho,Xlist,k,dmon=sdm0):
    L=[p]
    for j in range(k):
        p=sdm(p,F,rho,Xlist,dmon)
        L.append(Poly(p,Xlist)/1)
    sigma={Xlist[i]:rho[i] for i in range(len(Xlist)) }
    return L, [Poly(e,Xlist).subs(sigma)/1 for e in L]
    
def checkeq(p,F,rho,Xlist,dmon=sdm0):
    sigma={Xlist[i]:rho[i] for i in range(len(Xlist)) }    
    if p==0:
        return True
    if p.subs(sigma)!=0:
        return False
    L=[p]
    G=[0]    
    while(true):
        p=sdm(p,F,rho,Xlist,dmon)
        if p.subs(sigma)!=0:
            return G,False
        G=groebner(L,Xlist,domain='QQ')
        cf,r=reduced(p,G)
        if r==0:
            return L,p/1,G,cf,True
        L.append(p/1)

def grad(f, Xlist):
    fe=Add(f,0)
    g=[]
    for x in Xlist:
        g.append(fe.diff(x))
    return(g)
        
def lie(p,F,Xlist):
    i=0
    Nabla = grad(p,Xlist)
    e=0
    for d in Nabla:        
        e=e+d*F[i]
        i=i+1
    return(expand(e))

def logist(r,x0,n=4):
    L=[x0]
    x=x0
    for i in range(n):
        x=r*(x-x**2)
        L.append(x)
    return L
 
zero, unit,left,right,dleft,dright,left0,right0,x=var('zero, unit,left,right,dleft,dright,left0,right0,x')   

pc=dleft*right + dright*left0
pc2= left*dright+right*dleft-x*dleft*dright
psh=dleft*right + dright*left
ph=dleft*dright
pinf=left*dright+right*dleft+dleft*dright

def dermon(m,F,rho,q,p,Xlist):  
    if m==():
        return 0
    if m==(0,):
        return Add(q,0).subs({unit:1,zero:0})  
    if m==(1,):
        return F[0]
    if len(m)==1:
        d=dermon((m[0]-1,),F,rho,q,p,Xlist)
        return p.subs({left:Xlist[0], dleft:F[0], right:Poly({(m[0]-1,):1},Xlist)/1, dright:d, left0:rho[0], right0:rho[0]})
    m0=m[:1]
    m1=m[1:]
    d0=dermon(m0,F[:1],rho[:1],q,p,Xlist[:1])
    d1=dermon(m1,F[1:],rho[1:],q,p,Xlist[1:])
    sigma0={Xlist[i]:rho[i] for i in range(len(Xlist))}
    p1=Poly({m1:1},Xlist[1:])/1
    return Poly(p.subs({left:Xlist[0]**m[0], dleft:d0, right:p1, dright:d1, left0:rho[0]**m[0], right0:p1.subs(sigma0)  }),Xlist)/1
    
dmLie=lambda m,F,rho,Xlist: dermon(m,F,rho,zero,psh,Xlist)    
dmCon=lambda m,F,rho,Xlist: dermon(m,F,rho,zero,pc,Xlist)    
dmHad=lambda m,F,rho,Xlist: dermon(m,F,rho,unit,ph,Xlist)

dmCon2=lambda m,F,rho,Xlist: dermon(m,F,rho,zero,pc2,Xlist)
dmInf= lambda m,F,rho,Xlist: dermon(m,F,rho,zero,pinf,Xlist)

