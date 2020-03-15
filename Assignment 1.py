print("**** ISA calculator ****")
print()

def menu () :
    print("1. Calculate ISA for altitude in meter \n2. Calculate ISA for altitude in feet \n3. Calculate ISA for altitude in FL")
    print()
menu()
choice=input("Enter your choice: ")
print()
if choice=="1":
    h=float(input("Enter altitude: "))
elif choice=="2":
    h=float(input("Enter altitude: "))*0.3048
elif choice=="3":
    h=float(input("Enter altitude: "))*30.48

from math import exp

h0=[0,11000,20000,32000,47000,51000,71000,86000]
T=[0 for i in range(8)]
a=[-0.0065,0,0.001,0.0028,0,-0.0028,-0.002]
p=[0 for i in range(8)]
i=0
T[0]=288.15
p[0]=101325
while i<7:
    T[i+1]=T[i]+a[i]*(h0[i+1]-h0[i])
    if i!=1 and i!=4:
        p[i+1]=p[i]*((T[i+1]/T[i])**(-9.80665/(a[i]*287)))
    else:
        p[i+1]=p[i]*exp(-9.80665*(h0[i+1]-h0[i])/(T[i]*287))
    i=i+1    

if h<=h0[7]:
    if h<=h0[1]:
        T1=T[0]+a[0]*(h-h0[0])
        p1=p[0]*((T1/T[0])**(-9.80665/(a[0]*287)))

    elif h<=h0[2]:
        T1=T[1]+a[1]*(h-h0[1])
        p1=p[1]*exp(-9.80665*(h-h0[1])/(T1*287))

    elif h<=h0[3]:
        T1=T[2]+a[2]*(h-h0[2])
        p1=p[2]*((T1/T[2])**(-9.80665/(a[2]*287)))
    
    elif h<=h0[4]:
        T1=T[3]+a[3]*(h-h0[3])
        p1=p[3]*((T1/T[3])**(-9.80665/(a[3]*287)))
    
    elif h<=h0[5]:
        T1=T[4]+a[4]*(h-h0[4])
        p1=p[4]*exp(-9.80665*(h-h0[4])/(T1*287))
    
    elif h<=h0[6]:
        T1=T[5]+a[5]*(h-h0[5])
        p1=p[5]*((T1/T[5])**(-9.80665/(a[5]*287)))
    
    elif h<=h0[7]:
        T1=T[6]+a[6]*(h-h0[6])
        p1=p[6]*((T1/T[6])**(-9.80665/(a[6]*287)))

    rho=p1/(287*T1)

    print()
    print("Temperature: ",round(T1,2),"K (",round(T1-228.15,1),"'C)")
    print("Pressure: ",round(p1,2),"Pa (",round((p1/101325)*100,1),"% SL)")
    print("Density: ",round(rho,4),"kg/m3 (",round((rho/1.225)*100,1),"% SL)")
    
else :
    print()
    print("I'm sorry, you are out of reach of the atmosphere.")

print()
print("Ready.")

