import utilities as li

# coefficients of olynomial
a=[1,-3,-7,27,-18]

i=li.Laguerres_method(a,-3.5,10**(-4),4)
print("root are")
for j in i:
    print(j)
    j=j+1
print("error value +/-",10**(-16))

# output:
#         root are
#         -3.0
#         1.0
#         2.0000000000000004
#         2.9999999999999996
#         error value +/- 1e-16
