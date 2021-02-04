######################## Examination ##############################
import library
import math

#Question no. 1- Wein's displacement law

#Define the function given
def f1(x):
    return ((x-5) * (math.exp(x))) + 5
    
#Calculation of 'b'
def Wein(x):
    y = ((6.626*10**(-34))*(3*10**8))/((1.381*10**(-23))*x)
    return y  
    

    
    
#Main function

A,B=library.newton_raphson(f1,2,100)
print('The value of x for the given equation is coming out to be: ')
#print(A)
#print(B[len(B)-1])
for n, v in zip(A, B):
        print("{} = {}".format(n, v))
        
#Call function to calculate b 
ans=Wein(B[len(B)-1])
print('The Weins cosntant b came out to be: ')
print(ans)
        

        


