#Question 2-integration in angular variables

import library
import math

#Define the given function to integrate
def Timeperiod(x):
    z=(4*0.3193)/((1 - 0.1464*(math.sin(x))**2)**(0.5))
    return z 
    
#Main function

#Call simpson function from the library
ans=library.simpson(0,1.57,Timeperiod,10)
print('the Time period of the given Pendulum is: ')
print(ans,'sec inverse')



