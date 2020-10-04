import utilities as li
import math as m
import matplotlib.pyplot as plt
def equation(x):
    f=-x-m.cos(x)
    return(f)
func=equation
a,b=li.bracketing(-1,-1.5,func)
print("a for bisection",a)
print("b for bisection",b)
r1,r2,iter,err=li.bisection(a,b,func,"Q1b_bisection")
print("root is from bisection",r1,"+/-", 10**(-6))
li.plotting(iter,err,"iteration","absolute error","Q1b_bisection")

c,num,error=li.Regula_falsi(-1,-1.5,func,"Q1b_Regula.f")
print("root from regula_falsi a and b ",c,"+/-", 10**(-6))
li.plotting(num,error,"iteration","absolute error","Q1b_Regula.f")

#Newton_Raphson method
x0,N,abse=li.Newton_Raphson(func,-1.5,0.01,"Q1b_Newton.R")
print("root from Newton_Raphson ",x0,"+/-", 10**(-6))
li.plotting(N,abse,"iterations","absolute error","Q1b_Newton.R")


# Outputa for bisection -0.375
# b for bisection -1.5
# root is from bisection -0.7390846610069275 +/- 1e-06
# root from regula_falsi a and b  -0.7390850477639415 +/- 1e-06
#root from Newton_Raphson  -0.7390851332149027 +/- 1e-06
