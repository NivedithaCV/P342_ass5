import utilities as li
import math as m
import matplotlib.pyplot as plt

def equation(x):
    f=m.log(x)-m.sin(x)
    return(f)
func=equation

a,b=li.bracketing(1.5,2.5,func)
print("a bracket is",a)
print("b bracket is",b)

#Bisection method
r1,r2,iter,err=li.bisection(a,b,func,"Q1_bisection")
print("root from bisection ",r1,"+/-", 10**(-6))
li.plotting(iter,err,"iterations","absolute error","Q1_bisection")

#Regula_falsi method
c,num,error=li.Regula_falsi(1.5,2.5,func,"Q1_Regula.f")
print("root from regula_falsi ",c,"+/-", 10**(-6))
li.plotting(num,error,"iterations","absolute error","Q1_Regula.f")

#Newton_Raphson method
x0,N,abse=li.Newton_Raphson(func,1.5,0.01,"Q1_Newton.R")
print("root from Newton_Raphson ",x0,"+/-", 10**(-6))
li.plotting(N,abse,"iterations","absolute error","Q1_Newton.R")

# Output
# a bracket is 1.5
# b bracket is 2.5
# root from bisection  2.219106674194336 +/- 1e-06
# root from regula_falsi  2.2191071418525734 +/-1.1e-06
# root from Newton_Raphson  2.2191071418525734 +/- 1.1e-06
