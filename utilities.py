
import math as m
import matplotlib.pyplot as plt

#read  and write
def read_write(A,z,name):
    if z=="r":
        with open(A,z) as fhand:
            M=[]
            N=[]
            for line in fhand:
              line=line.rstrip()
              li=list(line.split(","))
              c=len(li)
              M.append(li)
            r=len(M)
            c=len(M[0])
            A=[[0 for y in range(c)]for x in range(r)]

            for i in range(r):
              for j in range(c):
                  A[i][j]=int(M[i][j])
            return(A,c)
    if z=="w" or z=="a":
        file1=open(name,z)
        file1.write(A)
        file1.close()


#partial pivoting
def partial_pivot(a,b,col):
    """This function does partial pivoting of passed matrices"""
    r=len(a)
    for i in range(r):
        if a[i][i]==0:
            for k in range(i,col):
                if k==i or a[i][i]!=0:
                    continue
                else:
                    if abs(a[k][i])>abs(a[i][i]):
                        # c=b[i][col-1]
                        # b[i][col-1]=b[k][col-1]
                        # b[k][col-1]=c
                        for j in range(r):
                            pivot=a[i][j]
                            a[i][j]=a[k][j]
                            a[k][j]=pivot
                        for z in range(col):
                            c=b[i][z]
                            b[i][z]=b[k][z]
                            b[k][z]=c
    return a,b


#Gauss_Jordan elemination
def Gauss_Jordan(a,b,col):
    """Gauss Jordan method of decomposition"""
    for q in range(r):
        pivot=a[q][q]
        for l in range(q,r):
            a[q][l]= a[q][l]/pivot
            b[q][col]=b[q][col]/pivot
        for w in range(r):
            if a[w][q]==0 or q==w:
                continue
            else:
                factor=a[w][q]
                b[w][col]=b[w][col]-factor*b[q][col]
                for c in range(q,r):
                    a[w][c]=a[w][c]-factor*a[q][c]

    return a,b


#multiplication of matrice
def mult(M,N,col_N):
    r=len(M)
    E=[[0 for y in range(col_N)]for x in range(r)]
    I=[[1,0,0],[0,1,0],[0,0,1]]
    for i in range(r):
        for j in range(col_M):
            for k in range(r):
                E[i][j]+=M[i][k]*N[k][j]
    if E==I:
        print("result is identity matrix")
    for r in E:
        print(r)


    def L_Udec(A):
        for j in range(c_A):
            for i in range(len(A)):

                #diagonal
                if i==j:
                    sum=0
                    for u in range(i):
                        sum=sum+A[i][u]*A[u][i]
                    A[i][i]=A[i][i]-sum

                    #elements of upper triangle
                if i<j:
                    sum=0
                    for k in range(i):
                        sum=sum+A[i][k]*A[k][j]
                    A[i][j]=A[i][j]-sum

                    #elements of lower triangle
                if i>j:
                    sum=0
                    for z in range(j):
                        sum=sum+A[i][z]*A[z][j]
                    A[i][j]=(A[i][j]-sum)/A[j][j]
        return(A)


def forw_backw(A,B,col):
    r=len(A)
    Y=[[0 for x in range(col)] for y in range(r)]
    X=[[0 for t in range(col)] for w in range(r)]
    for g in range(col):
        #forwardd substitution
        for i in range(r):
            sum=0
            for k in range(i):
                sum=sum+A[i][k]*Y[k][g]
            Y[i][g]=B[i][g]-sum

        #backward substitution
        for l in range(r-1,-1,-1):
            sum=0
            for m in range(l+1,r):
                sum=sum+A[l][m]*X[m][g]
            X[l][g]=(Y[l][g]-sum)/A[l][l]
            X[l][g]=round(X[l][g],4)
    print("matrix Y",Y)
    print("inverse matrix is",X)


# bracketing of roots
def bracketing(a,b,equation):
    f_a=equation(a)
    f_b=equation(b)
    if equation(a)*equation(b)<0:
        return(a,b)

    if f_a*f_b>0:
        i=0
        while equation(a)*equation(b)>0 and  i<15:
            if abs(f_a)<abs(f_b):
                    a=a-0.5*(b-a)
                    i=i+1
            else:
                b=b+0.5*(b-a)
                i=i+1

        if equation(a)*equation(b)<0:
            return(a,b)
        if i>14:
            return("please provide another range")


# bisection method
def bisection(a,b,func,A):
    num=[]
    error=[]
    s="{:<16s}{:s}\n".format("iteration","absolute error")
    read_write(s,'w',A)
    i=0;
    while abs(a-b)>=0.000001 and i<200 :
        c=(a+b)/2
        if func(a)*func(c)<0:
            b=c
            i=i+1
            num.append(i)
            error.append(abs(a-b))
            s="{:<16d}{}\n".format(i,abs(a-b))
            read_write(s,'a',A)
        if func(b)*func(c)<0:
            a=c
            i=i+1
            num.append(i)
            error.append(abs(a-b))
            s="{:<16d}{}\n".format(i,abs(a-b))
            read_write(s,'a',A)
    if i<200:
        return(a,b,num,error)
    else:
        return("Please provide another range")


#regula_falsi method
def Regula_falsi(a,b,func,A):
    n=[]
    e=[]
    a,b=bracketing(a,b,func)
    i=0;k=0;error=1
    s="{:<16s}{:s}\n".format("iteration","absolute error")
    read_write(s,'w',A)
    while abs(a-b)>=0.000001 and i<200 and abs(error)>=0.000001:
        c=b-(((b-a)*func(b))/(func(b)-func(a)))
        error=k-c
        if func(a)*func(c)<0:
            b=c
            k=c
            i=i+1
            if i!=1:
                n.append(i-1)
                e.append(abs(error))
                s="{:<16d}{}\n".format(i-1,abs(error))
                read_write(s,'a',A)
        if func(a)*func(c)>0:
            a=c
            k=c
            i=i+1
            if i!=1:
                n.append(i-1)
                e.append(abs(error))
                s="{:<16d}{}\n".format(i-1,abs(error))
                read_write(s,'a',A)
        if func(c)==0:
            return(c)
    if abs(a-b)<0.000001:
        return(a,n,e)
    if i<200:
        return(c,n,e)
    else:
        print("Please provide another range")


# derivations
def first_derivative(func,x0,h):
    f_=(func(x0+h)-func(x0-h))/(2*h)
    return(f_)

def derivative1(a,func,x0,n,h):
    f_=(func(a,x0+h,n)-func(a,x0-h,n))/(2*h)
    return(f_)

# sec_derivative function
def sec_derivative(func,x0,h):
    f2=(func(x0+h)-2*func(x0)+func(x0-h))/(h*h)
    return(f2)
def derivative2(a,func,x0,n,h):
    f2=(func(a,x0+h,n)-2*func(a,x0,n)+func(a,x0-h,n))/(h*h)
    return(f2)

# Newton_Raphson method
def Newton_Raphson(func,x0,h,A):
    i=0;error=abs(x0);
    N=[]
    ab_er=[]
    s="{:<16s}{:s}\n".format("iteration","absolute error")
    read_write(s,'w',A)
    while error>=0.000001 and i<=200 and func(x0)!=0:
        k=x0
        f_=first_derivative(func,x0,h)
        x0=x0-func(x0)/f_
        error=abs(x0-k)
        i=i+1
        N.append(i)
        ab_er.append(error)
        s="{:<16d}{}\n".format(i,error)
        read_write(s,'a',A)
    if i<200:
        return(x0,N,ab_er)
    else:
        print("Please provide another range")


# potting function
def plotting(n,e,x,y,u):
    plt.plot(n,e,marker=".")
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(u)
    plt.show()

#polynomial calculator
def polynomial(a,x,n):
    j=0;f=0
    while j<=(n):
        t=a[j]*x**(n-j)
        f=f+t
        j=j+1
    return(f)


# deflection using synthetic division
def deflection(poly,x0,n):
    k=0; d=[0 for y in range(n+1)]
    while k<=n:
        if k==0:
            d[k]=poly[k]
            k=k+1
        else:
            d[k]=(x0*d[k-1])+poly[k]
            k=k+1
    return(d,n-1)

#Laguerres_method function
def Laguerres_method(c,x0,h,n):
    X=[0 for m in range(n)]; z=1;g=n
    while z<=4:
        func=polynomial
        f=func(c,x0,n)
        error=1;i=0
        while error>=0.000001 and i<=200 and func(c,x0,n)!=0:
            G=derivative1(c,func,x0,n,h)/func(c,x0,n)
            H=G**2-(derivative2(c,func,x0,n,h)/func(c,x0,n))
            if G<0:
                x_=x0
                a=n/(G-m.sqrt((n-1)*(n*H-G**2)))
                x0=x_-a
                error=abs(x0-x_)
                i=i+1
            else:
                x_=x0
                a=n/(G+m.sqrt((n-1)*(n*H-G**2)))
                x0=x_-a
                error=abs(x0-x_)
                i=i+1
        x1=x0
        X[z-1]=x0
        c,n=deflection(c,x0,n)
        z=z+1
    return(X)
