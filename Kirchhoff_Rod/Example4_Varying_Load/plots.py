import numpy as np
import matplotlib.pyplot as plt
import sys        

filename=sys.argv[1]
f=open(filename,"r")

il=1
XX=[]
YY=[]
ZZ=[]
Folds=[]
Tip_x=[]
Tip_y=[]
Tip_z=[]

S=[]
X1=[]
Y1=[]
Z1=[]
X1_end=[];
Y1_end=[];
Z1_end=[];
LP=[]
ii=0
th=2*np.pi
Nl=14/7+1
rib=0.1
arm_x=[]
arm_y=[]
arm_z=[]
Arm=[]
bif=[]
lam=[]
while(1):
    s=[]
    x1=[]
    y1=[]
    z1=[]
    D1=[]
    D2=[]
    D3=[]
    A=f.readline()
    if (A==''):
        break
    K=A.split()
    ii=ii+1
    if int(K[2])==5:
        LP.append(ii)
    for i in range(0,int(K[8])-2):
        B=f.readline()
        Bb=B.split()
        if (i%Nl==0 and i< (int(K[6])*Nl)):    
            s.append(float(Bb[0]))
            x1.append(float(Bb[1]))
            y1.append(float(Bb[2]))
            z1.append(float(Bb[3]))
            q1=(float(Bb[4]))
            q2=(float(Bb[5]))
            q3=(float(Bb[6]))
            r1=float(Bb[1])
            r2=float(Bb[2])
            r3=float(Bb[3])
            
            if (i==((int(K[6])-1)*3) + 0):
                q11=(float(Bb[4]))
                q12=(float(Bb[5]))
                q13=(float(Bb[6]))
                        
        
                
        if (i%Nl==1 and i< (int(K[6])*Nl)):    
            q4=(float(Bb[0]))
            d1=[q1*q1 - q2*q2 - q3*q3 + q4*q4,2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4)]
            d2=[2*(q1*q2 - q3*q4),-q1*q1 + q2*q2 - q3*q3 + q4*q4,2*(q2*q3 + q1*q4)]
            d3=[2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4),-q1*q1 - q2*q2 + q3*q3 + q4*q4]
            D1.append(d1);
            D2.append(d2);
            D3.append(d3);
            
            
            if (i==((int(K[6])-1)*3) + 1):
                q14=(float(Bb[0]))
                mu11=(float(Bb[1]))
                mu12=(float(Bb[2]))
                mu13=(float(Bb[3]))
                mu14=(float(Bb[4]))
                n1=(float(Bb[5]))
                n2=(float(Bb[6]))
                
        if (i==((int(K[6])-1)*3) + 0):
            #print(Bb[0])
            r3=(float(Bb[3]))    
                
        if (i%Nl==2 and i< (int(K[6])*Nl)):
            if (i==((int(K[6])-1)*3) + 2):
                n3=(float(Bb[0]))
            
    B=f.readline();
    Bb=B.split()
    load=float(Bb[3])
    arm=float(Bb[6])
    B=f.readline();
    Bb=B.split()
    Tht=float(Bb[0])
    lam.append(load)
    alf=float(Bb[1])
    

    Arm.append([r1 + arm*d2[0]*np.sin(Tht)+ arm*d3[0]*np.cos(Tht),r2 +arm*d2[1]*np.sin(Tht)  + arm*d3[1]*np.cos(Tht),r3 + arm*d2[2]*np.sin(Tht)  + arm*d3[2]*np.cos(Tht)]);

    N_norm=np.sqrt( (n1*n1) + (n2*n2) + (n3*n3))
    bif.append((0.0*(r1 + arm*0.0*d1[0] + arm*np.sin(Tht)*d2[0] + arm*np.cos(Tht)*d3[0])) + (np.sin(-alf)*(r2 + arm*0.0*d1[1] + arm*np.sin(Tht)*d2[1] + arm*np.cos(Tht)*d3[1])) + (np.cos(-alf)*(r3 + arm*0.0*d1[2] + arm*np.sin(Tht)*d2[2] + arm*np.cos(Tht)*d3[2])) )
    X1.append(x1);
    Y1.append(y1);
    Z1.append(z1);
    S.append(s);

for k in range(0,len(X1)):
    X1_end.append(X1[k][-1]);
    Y1_end.append(Y1[k][-1]);
    Z1_end.append(Z1[k][-1]);
Tip_x.append(X1_end)
Tip_y.append(Y1_end)
Tip_z.append(Z1_end)

print("!!!Successfully Loaded solutions");    


plt.plot(lam,bif,color='k')
plt.scatter([lam[0]],bif[0],color='b')
plt.xlim([-30,30])
plt.xlabel('Load P')
plt.ylabel('$\\gamma$')
plt.grid()
plt.show()
