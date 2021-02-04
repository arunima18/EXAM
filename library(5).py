import math
import random


# First Derivative of a Function
def d1(f, x):
    h = 10 ** (-5)
    return (f(x + h) - f(x - h)) / (2 * h)


#partial pivoting
def partial_pivot(Mat1,n,Mat2):
    
    i=0
    j=0
    while i<=(n-1):
        if Mat1[i][i]==0:
            j=i
            
            while j <= n:
                
                if Mat1[j+1][i]!=0:
                    Mat2[j]=Mat1[j+1]
                    Mat1[j+1]=Mat1[i]
                    Mat1[i]=Mat2[j]
                    j=(n+2)
                else:
                    j=j+1
        i=i+1
    return Mat1 
   
   
   
#crout's method
def crout(MatA):
    i=0
    j=0
    k=0
    sum1=0
    sum2=0
    U=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    L=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    
    while j<=3:
        
        i=0
        while i<=3:
            sum1=0
            sum2=0
            k=0
            if i ==0:
                U[i][j]=MatA[i][j]
                i=i+1
                
            elif i>j:
                #if k<=(j-1):
                while k<=(j-1):
                    sum1=sum1-(L[i][k]*U[k][j])/U[j][j]
                        #L[i][j]=(MatA[i][j]-L[i][k]*U[k][j])/U[j][j]
                    k=k+1
                L[i][j]=(MatA[i][j]/U[j][j])+sum1
                
                   
                #else:
                 #   L[i][j]=MatA[i][j]/U[j][j]
                i=i+1
            elif i<=j and i !=0:
               
                while k<=(i-1):
                    sum2=sum2-(L[i][k]*U[k][j])
                   
                        #L[i][j]=(MatA[i][j]-L[i][k]*U[k][j])/U[j][j]
                    k=k+1
                U[i][j]=(MatA[i][j])+sum2
                
                i=i+1
                   
               
        j=j+1
    return L ,U 
    
    

# gauss-jordan method
def gj(am):
    n = len(am)
    # diagonal = 1
    for r in range(n):
        partpiv(am)

        pivot = am[r][r]
        for j in range(r, 2 * n):
            am[r][j] /= pivot

        # non-diagonals = 0
        for i in range(n):
            if i == r or am[i][r] == 0:
                continue
            else:
                temp = am[i][r]
                for j in range(r, 2 * n):
                    am[i][j] -= temp * am[r][j]

    return am
    
    
    
    
#Solution for Ly=b: forward substituion 
def forward_sub(L,y,b):
    i=0
    j=0
    
    sumf=0
    while i<=3:
        j=0
        sumf=0
        while j<=(i-1):
            sumf=sumf-L[i][j]*y[j]
            j=j+1 
        y[i]=b[i]+sumf 
        #print('***')
        #print(y[i])
        #print('***')
        i=i+1 
        
    return y 
            
        
    
#Solution for Ux=y:backward substitution 
def backward_sub(U,x,y,n):
    i=2 
    j=0
    sumb=0
    while i>=(-1):
        sumb=0
        j=i+1 
        while j<=(n-1):
            sumb=sumb-(U[i+1][j+1]*x[j+1])/U[i+1][i+1]
            j=j+1 
        x[i+1]=(y[i+1]/U[i+1][i+1]) + sumb
        
        i=i-1 
    return x 
    
    
#LU decoposition    
def LU_dec(A):
    n = len(A)
    L = [[0 for e in range(n)] for f in range(n)]
    U = [[0 for g in range(n)] for h in range(n)]
    
    # solve column-wise.
    for j in range(n):
        # main diagonal elements of L are 1
        L[j][j] = 1
        # elements of you
        # i <= j
        for i in range(j + 1):
            temp = 0
            for k in range(i):
                temp += U[k][j] * L[i][k]
            U[i][j] = A[i][j] - temp
            del temp
            # replace elements of A
            A[i][j] = U[i][j]
    
        # elements of L
        # i > j
        for i in range(j + 1, n):
            temp = 0
            for k in range(i):
                temp += U[k][j] * L[i][k]
            L[i][j] = (A[i][j] - temp) / U[j][j]
            del temp
            # replace elements of A
            A[i][j] = L[i][j]
    
    return L, U, A
    
    
def for_backward(A, B):
    n = len(A)
    # forward substitution
    # solve row-wise
    for i in range(n):
        # j < i, lower triangular
        for j in range(i):
            B[i] -= A[i][j] * B[j]

    # backward substitution
    # run loop in decreasing value of index
    for i in range(n - 1, -1, -1):
        # j > i, upper triangular
        for j in range(i + 1, n):
            B[i] -= A[i][j] * B[j]
        B[i] /= A[i][i]

    return B

    

#matrix multiplication    
def multi(Mat1,Mat2,arr,n):
    i=0
    k=0
    j=0
    sum=0
    #arr=[[0,0,0],[0,0,0],[0,0,0]]
    while i<=n:
        j=0
        while j<=n:
            k=0
            sum=0
            while k<=n:
                sum=sum+ Mat1[i][k]*Mat2[k][j]
                k=k+1
                arr[i][j]=sum
            j=j+1
        i=i+1
    return arr 
    
    
    
################################################################################
################### Linear equations ###########################################

# bisection method

def bisec(f, a, b, max_iter):
    root_in_each_iter = []

    for i in range(max_iter):
        root = (a + b) * .5

        if abs(a - b) < 10 ** (-6):
            continue

        if f(a) * f(b) > 0:
            error_statement = 'Root not bracketed properly. Try again. '
            return error_statement

        if f(a) * f(root) > 0:
            a = root
        elif f(a) * f(root) < 0:
            b = root
        elif f(root) == 0:
            return root

        root_in_each_iter.append(root)

    return root_in_each_iter
    
    
# Regula Falsi (False-Position) Method


def regula_falsi(f, a, b, max_iter):
    root_in_each_iter = []

    for i in range(max_iter):
        root = b - (((b - a) * f(b)) / (f(b) - f(a)))

        if abs(a - b) < 10 ** (-6):
            continue

        if f(a) * f(b) > 0:
            error_statement = 'Root not bracketed properly. Try again. '
            return error_statement

        if f(a) * f(root) > 0:
            a = root
        elif f(a) * f(root) <= 0:
            b = root
        elif f(root) == 0:
            return root

        root_in_each_iter.append(root)

    return root_in_each_iter



# Newton-Raphson Method


def newton_raphson(f, guess, max_iter):
    x_list = []
    x = guess
    A = []
    for i in range(max_iter):
        h = f(x) / d1(f, x)

        if abs(h) < 10 ** (-6):
            continue

        x -= h
        x_list.append(x)
        A.append(i)

    return A,x_list

    
    
################################################################################
######################## Integration ###########################################


#Rectangle method for integration
def Rect_(a,b,f,N):
    i=1
    x=0
    y=0
    sum1=0
    h=(b-a)/N
    while i<=N:
        
        x=a+((2*i-1)*(h/2))
        y=f(x)
        sum1=sum1+(h*y) 
        #print(sum1)
        i=i+1 
    return sum1 
    
    
    
def trapezium(a,b,f,N):
    # f = integrand, lowr_lim = lower limit, uppr_lim = upper limit, N = no. of partitions

    
    sum2 = 0

    h = (b - a) / N

    for i in range(1, N):
        
        sum2 = sum2 + h * f(a + i * h)
        #print(h/2 * (f(a) + f(b)) + sum2)

    return h/2 * (f(a) + f(b)) + sum2
    
    
#Simpson's Method 
def simpson(a,b,f,N):
    h = (b - a) / N
    sum3 = 0
    x = a + h
    for i in range(1, int(N / 2) + 1):
        
        
        sum3 += 4 * f(x)
        #print(h/3 * (f(a) + f(b) + sum3))
        x += 2 * h
    x = a + 2 * h
    for i in range(1, int(N / 2)):
        
        
        sum3 += 2 * f(x)
        #print(h/3 * (f(a) + f(b) + sum3))
        x += 2 * h

    return h/3 * (f(a) + f(b) + sum3)
    

#Monte Carlo method
def montecarlo(f,a,b,N):
    x=[]
    sum=0.0

    for i in range(N):
        n=random.random()       #random number in the range (0,1)
        r=a+((b-a)*n)           #random number in range (a,b)
        x.append(r)

    for i in range(N):
        sum+=f(x[i])

    value=((b-a)*sum)/N

    return value
    
    
    
################################################################################
############################## D.E #############################################
    

#Euler's method for solving D.E    
def Euler(f,h,x_0,y_0,n):
    y=y_0
    x=x_0
    i=1 
    arr1=[x_0]
    arr2=[y_0]
    while i<=n:
        y=y+(h*f(x,y))
        x=x+h
        arr1.append(x)
        arr2.append(y)
        i=i+1 
    return arr1,arr2
    

def rungeKutta(f1,f2,x0, y0,z0, x_lim, h): #input- initial value problem x0, y(x0), y'(x0), end point, interval width
    # Blank lists for storing the coordinate values
    X=[]
    Y=[]
    n = (int)((x_lim - x0) / h) #defining number of intervals
    #assigning values to variables and storing them in Blank lists
    x = x0
    X.append(x0)
    y = y0
    Y.append(y0)
    z = z0
    for i in range(1, n + 1):
        x = x + h #updating the value of x per each iteration
        X.append(x) #storing the x value

        #defining the parameters for 4th order Runge_Kutta method
        ky1 = f1(x,y,z)
        kz1 = f2(x,y,z)

        ky2 = f1(x + (0.5 * h), y + (0.5 * ky1 * h), z + (0.5 * kz1 * h))
        kz2=  f2(x + (0.5 * h), y + (0.5 * ky1 * h), z + (0.5 * kz1 * h))

        ky3 = f1(x + (0.5 * h), y + (0.5 * ky2 * h), z + (0.5 * kz2 * h))
        kz3 = f2(x + 0.5 * h, y + 0.5 * ky2 * h, z + 0.5 * kz2 * h)

        ky4 = f1(x + h, y + ky3 * h, z + kz3 * h)
        kz4 = f2(x + h, y + ky3 * h, z + kz3 * h)

        # Update next value of x,y,z
        y = y + (h / 6.0) * (ky1 + 2 * ky2 + 2 * ky3 + ky4)
        z= z+ (h / 6.0) * (kz1 + 2 * kz2 + 2 * kz3 + kz4)
        Y.append(y) #updating the value of y per each iteration

    return X,Y
def shooting(f1,f2,x0,y0,xf,yf,z_g,h):
    tol=10**-6
    X,Y=rungeKutta(f1,f2,x0,y0,z_g,xf,h)
    yn=Y.pop() 
    if abs(yn-yf)>tol: 
        if yn>yf: #if the obtained y(x_f) overshoots the analytical boundary value
            z_high=z_g 
            print('Overshooting Guess:',z_high)
            y_high=yn

            while yn>=yf: 
                z_g=z_g-1
                X,Y=rungeKutta(f1,f2,x0,y0,z_g,xf,h)
                yn = Y.pop()

            z_l=z_g #storing the undershoot guess
            print('Obtained Undershooting Guess:', z_l)
            yl=yn
        elif yn<=yf: #if the obtained y(x_f) undershoots the analytical boundary value
            z_l=z_g
            print('Undershooting Guess:', z_l)
            yl=yn
            while yn<=yf: #running loop by increasing the guess value by 1 each time, till a overshooting guess is obtained
                z_g=z_g+1
                X,Y=rungeKutta(f1,f2,x0,y0,z_g,xf,h)
                yn = Y.pop()

            z_high=z_g #storing the undershoot guess
            print('Obtained Overshooting Guess:', z_high)
            y_high=yn
        #after having the guesses bracketing the analytical solution we find the actual value of slope at the initial point
        z= z_l + (z_high - z_l) * (yf-yl)/(y_high-yl) #lagrange interpolation

        X,Y=rungeKutta(f1,f2,x0,y0,z,xf,h)
        y_i=Y.pop() 
        while (y_i-yf)>tol: 
            z = z_l + (z_high - z_l) * (yf - yl) / (y_high - yl)
            X, Y = rungeKutta(f1, f2, x0, y0, z, xf, h)
            y_i = Y.pop()
        return rungeKutta(f1, f2, x0, y0, z, xf, h)
    else: #the guess bang on lands at the solution
        print("The guess z=",z_g," is perfect!")
        return rungeKutta(f1, f2, x0, y0, z_g, xf, h)
        
        
        
        
################################################################################
######################### Square fitting #######################################

#Average value of an array 
def Avg(X):
    return sum(X)/len(X)
        

def linear_fit(X,Y):

    if len(X)!=len(Y):
        print('Entries Not valid!')
    else:
        n=len(X)
        xbar=Avg(X)
        ybar=Avg(Y)
        S_xx=0
        S_xy = 0
        std_X=0
        std_Y=0
        for i in range(n):
            S_xx += X[i]**2-xbar**2
            S_xy += (X[i]*Y[i])-(xbar*ybar)
            std_X += ((X[i]-xbar)**2)/n
            std_Y += ((Y[i]-xbar)**2)/n
        b=S_xy/S_xx
        a=ybar-b*xbar
        sigma_x= math.sqrt(std_X)
        sigma_y =math.sqrt(std_Y)
        PearsonR=S_xy/((n)*sigma_x*sigma_y)
        return b,a,PearsonR
        
        





