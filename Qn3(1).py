#Question 3-Least square fitting
import library
import math 






#main function

#Part 1

X=[0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3,3.3]
Y=[2.2,1.96,1.72,1.53,1.36,1.22,1.1,1,0.86,0.75,0.65,0.6]

b1,a1,P1=library.linear_fit(X,Y)
print('The value of w_c for straight line fit came out to be: ')
print(b1)
print('The value of w_0 for straight line fit came out to be: ')
print(a1)
print('The pearson factor sor the staright line fit came out to be: ')
print(P1)

'''b2,a2,P2=library.linear_fit(X,Y)
print('The value of -w_c in case of exponential fit came out to be: ')
print(a2)
print('The value of w_0 for the exponential fit came out to be: ')
print(math.exp(b2))'''

#Pearson formula for exponential fit 
Z=[0 for i in range(len(Y))]
i=0 
while i<len(Y):
    Z[i]=math.log(Y[i],2.718)
    i=i+1 

b3,a3,P3=library.linear_fit(X,Z)
print('The value of w_c for the exponential plot came out to be: ')
print(b3)
print('The value of log(w_0) for the exponential fir came out to be: ')
print(a3)
print('The pearson factor for the exponential fit is: ')
print(P3)



################################################################################

'''
The value of w_0 for straight line fit came out to be: 
2.029102564102564
The value of w_c for straight line fit came out to be: 
-0.4747086247086248
-0.7655606603764613
The value of w_0 for the exponential plot came out to be: 
0.7903594771921305
The value of log(w_c) for the exponential fir came out to be: 
-0.39563719590281327
The pearson factor for the exponential fit is: 
-0.26146415804622253





'''








