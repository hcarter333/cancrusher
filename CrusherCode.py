
{{{id=34|
import numpy
from numpy.linalg import inv
rstvty(temperature)=(4.5+1.133*10^-2*temperature)*10^-8
#Returns the specific heat in J/Kg/K
specheat(temperature) = 4.18*(0.16+1.5e-4*temperature)*1000.0    
#setup the constants
perm0 = 1.25663706*1e-6
nmov = 6
nfix = 3

#need to change this for Pb
density = 2824.0

thickness = .00016

#Setup an array for coil dimensions
r = numpy.array([3.43e-2, 3.43e-2, 3.43e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2], dtype=float)
#Here's the array for coil z coordinates
z = numpy.array([-2.67e-03, 0.0e-00, 2.67e-3,-6.67e-3, -4.0E-3, -1.33E-3, 1.33E-3, 4.0E-3, 6.67E-3], dtype=float)
#Here's the array for coil-can separtion
dr = numpy.array([1.33E-3, 1.33E-3, 1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3], dtype=float)

delta = abs(z[nmov+1] - z[nmov+2])

#Here are two more arrays that are used, but I'm not sure what for yet
vr = numpy.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], dtype=float)
vz = numpy.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], dtype=float)

#NOTE:  r, z, and vr will all be modified by move can

#array stores the temerature of the can under each coil
temp = numpy.array([293.0, 293.0, 293.0, 293.0, 293.0, 293.0, 293.0], dtype=float)

#in the following, we'll need the resistances and capacitances that were pulled 
#from the input data file
#The following are the resistance, capacitance, inductance, admittance, 
#current, and charge arrays with one entry per can modelling coil
#They have been filled in here with the values from the model data file
#In the model data file, there are 7 values as opposed to the six I expected
#I need todo check on this
res = numpy.array([8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)

cap = numpy.array([0.00005, -1, -1, -1, -1, -1, -1], dtype=float)

lckt = numpy.array([0.05, 0, 0, 0, 0, 0, 0], dtype=float)

cur = numpy.array([0.0, 0, 0, 0, 0, 0, 0], dtype=float)

V0 = numpy.array([5000.0, 0, 0, 0, 0, 0, 0], dtype=float)

#setup the admittance matrix.  I've made this a single statement here, keep in mind that 
#there's a there is a loop treatment in the original code
adc = numpy.array([1/cap[0], 0, 0, 0, 0, 0, 0], dtype=float)
qq = numpy.array([V0[0]*cap[0], 0, 0, 0, 0, 0, 0], dtype=float)

#Here's a clue to how the radius positions of coil and can should be worked out
#rcan = r[nfix+1]
#the indexing above seems odd
rcan = r[nfix+1]
mass = 2*pi*rcan*density*thickness*delta

#the next step is to fill in the res array with calls to rstvty
#This begins to make more sense after a fahion.  The one bigger is intentional, the 
#first value in the data file was real, all the others are filled in
for i in range(1, nmov+1):
    res[i] = 1e3*2*pi*rcan*rstvty(temp[i])/(thickness*delta)


#arrays to hold the reduced mutual inductance matrix
mm = numpy.array(range((nmov+1)^2), dtype=float)
mm.shape=(nmov+1,nmov+1)
mmold = mm.copy()

#setting up the mutual inductance raw flux array
mfull = numpy.array(range((nmov+nfix)*(nmov+nfix)), dtype=float)
mfull.shape=(nmov+nfix,nmov+nfix)

#now setup the variables for stepping through the simulation
ntim = 600
ddt = 0.00002
nchange = 20

coilI = numpy.array(range(ntim+1), dtype=float)
bzero = numpy.array(range(ntim+1), dtype=float)
zeroline = numpy.array(range(ntim+1), dtype=float)
heatenrg = numpy.array(range(ntim+1), dtype=float)
work = numpy.array(range(ntim+1), dtype=float)
enrgtot = numpy.array(range(ntim+1), dtype=float)
ptime = numpy.array(range(ntim+1), dtype=float)
denrg = 0.0
dheat = 0.0

dt = ddt
time = 0.0
aaa = 1.0
ntim1 = ntim-1
#This concludes all the setup of variables that is done by the initialize function


#setup the flux calculation that mimics what goes in in dbcoil
#call it the same name as the original code for readability
#leave the unneccasary parts out though
t(rcoil, zc, rc) = (rc+rcoil)^2+zc^2
argm(rcoil, zc, rc)= 4*rcoil*rc/t(rcoil,zc,rc)
dbcoilflux(rcoil, zc, rc, curren) = (((2*perm0*curren*rcoil*rc)/argm(rcoil, zc, rc))/(t(rcoil, zc, rc)^0.5))*((2-argm(rcoil, zc, rc))*elliptic_kc(argm(rcoil, zc, rc))-(2*elliptic_ec(argm(rcoil, zc, rc))))

#here's the z direction magnetic field calculation from dbcoil
#There are two versions, one for r = 0 and the other for r <> 0
dbcoilbzrzero(rcoil, zc, rc, curren) = perm0*rcoil^2*curren/(2*(rcoil^2+zc^2)^(3/2))
bzfac(rcoil, zc, rc) = ((rcoil^2-rc^2)-zc^2)/(((rcoil-rc)^2)+zc^2)
dbcoilbz(rcoil, zc, rc, curren) = (bzfac(rcoil, zc, rc)*elliptic_ec(argm(rcoil, zc, rc))+elliptic_kc(argm(rcoil, zc, rc)))+((2*perm0/(2*pi))/(t(rcoil, zc, rc)^0.5))

#now take mfull created above and use it to generate the reduced mutual induction matrix mm
#mm is dimensioned for only the moving coils plus one extra row and column
#ultimately we'll need to track the immediately previous version of mm
#in the original code this is called mmold, and here too!

#The mutual inductance array has to be calculated multiple times, so I'm putting 
#it into a Python funcion, (I hope)
def find_mutual_inductance(mfullarray):
    #first for the off-diagonal elements
    for i in range(0, nmov+nfix-1):
       for j in range(i+1, nmov+nfix):
            mfullarray[i,j] = dbcoilflux(r[i], z[j]-z[i], r[j]-dr[j], 1)
    
    #next for the diagonal elements
    for i in range(0, nmov+nfix):
        mfullarray[i,i] = dbcoilflux(r[i], z[i]-z[i], r[i]-dr[i], 1)
   
    #finally copy the off diagonal elements to the other side of the array
    for i in range(1, nmov+nfix):
        for j in range(0, i):
            mfullarray[i,j] = mfull[j,i]
            
    return mfullarray


#also create a function here that encapsulates the reduced array creation
#since this will need to be called multiple times
def make_reduced_matrix(mreduced, mfullwork):
    global mmold
    mmold = mm.copy()
    jmwork = numpy.array(range(nmov), dtype=float)
    #nt contains the fixed coil portion of the mutual inductance array
    #it's actually just a 
    ntwork = mfull[0:(nfix), 0:(nfix)].copy()
    nttwork = sum(ntwork)
    
    #jt is used to total the mutual inductance associated with each moveable coil
    #the result is placed in the jmwork array
    for i in range(0, nmov):
        jtwork=mfull[i+nfix,0:nfix].copy()
        jmwork[i]=sum(jtwork)
    
    #Here's why mm was dimensioned one bigger than it should have been
    #I don't have a good reason for this yet other than it came out of the old code
    #notice, however that we do overwrite the 0 row and column next
    mreduced = mfullwork[nfix-1:nmov+nfix, nfix-1:nmov+nfix].copy()
    mreduced = 1e6*mreduced
    #print "mfullwork"
    #print mfullwork
    #print "mreduced"
    #print mreduced
    #mm = mm*1e6

    for i in range(1, nmov+1):
        mreduced[i,0] = 1e6*jmwork[i-1]
        mreduced[0,i] = 1e6*jmwork[i-1]
    
    mreduced[0,0] = 1e6*nttwork

    return mreduced
///
}}}

{{{id=30|
global mfull
mfull = find_mutual_inductance(mfull)

mm = make_reduced_matrix(mm, mfull)
///
}}}

{{{id=32|
mm
///
array([[ 1.00491907,  0.21781373,  0.27722786,  0.32387533,  0.32387533,
         0.27722786,  0.21781373],
       [ 0.21781373,  0.13039726,  0.09826708,  0.07426194,  0.05946457,
         0.04906901,  0.04127011],
       [ 0.27722786,  0.09826708,  0.13039726,  0.09826708,  0.07433023,
         0.05946457,  0.04906901],
       [ 0.32387533,  0.07426194,  0.09826708,  0.13039726,  0.09838592,
         0.07433023,  0.05946457],
       [ 0.32387533,  0.05946457,  0.07433023,  0.09838592,  0.13039726,
         0.09826708,  0.07426194],
       [ 0.27722786,  0.04906901,  0.05946457,  0.07433023,  0.09826708,
         0.13039726,  0.09826708],
       [ 0.21781373,  0.04127011,  0.04906901,  0.05946457,  0.07426194,
         0.09826708,  0.13039726]])
}}}

{{{id=31|
inv(mm)
///
array([[ 19.46621276,  -0.98195254,  -9.96888971, -19.46257346,
        -19.46257346,  -9.96888971,  -0.98195254],
       [ -0.98195254,  17.88356992, -12.8407891 ,   1.69472885,
          0.32025774,   0.26410608,  -0.34204355],
       [ -9.96888971, -12.8407891 ,  32.91489114,  -3.92050343,
         11.16755566,   4.6360251 ,   0.26410608],
       [-19.46257346,   1.69472885,  -3.92050343,  47.35792694,
          5.4720753 ,  11.16755566,   0.32025774],
       [-19.46257346,   0.32025774,  11.16755566,   5.4720753 ,
         47.35792694,  -3.92050343,   1.69472885],
       [ -9.96888971,   0.26410608,   4.6360251 ,  11.16755566,
         -3.92050343,  32.91489114, -12.8407891 ],
       [ -0.98195254,  -0.34204355,   0.26410608,   0.32025774,
          1.69472885, -12.8407891 ,  17.88356992]])
}}}

{{{id=11|
plot(dbcoilbzrzero(2, zc, 0.1, 1), 0, 5)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=2|
#here's what the resistivity function looks like
plot(rstvty(temperature), 293, 297)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=3|

///
}}}

{{{id=4|
#Now graph the flux along the z coordinate at the edge of a coil and see if it makes sense
plot(dbcoilflux(2, zc, 2, 1), 0, 5)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=16|
#Test along the radial direction in the plane of the coil
plot(dbcoilflux(2, 0, rc, 1), 2, 5)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=15|
#and inside the coil, being careful to avoid the point r = 0
plot(dbcoilflux(2, 0, rc, 1), 0.01, 2)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=18|
#The above was to provide the magnetic flux at one coil due to another for finding 
#the mutual inductance of a set of two coils.  The following models the coils and 
#calls the above functions to fill in the mutual induction matrix
#Utilize a python for loop here and I guess go ahead and use sage matrices, although 
#I'm not sure there's not a reason to use numpy matrices or plain old Python arrays
#The documentation for numpy arrays is the easiest for me to read, so I'm going that 
#direction.  I went ahead and put the import at the top of the worksheet
///
}}}

{{{id=17|
#Checking that mfull is filled up appropriately.  There were issues with the manner in which
mfull
///
array([[  3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05,   1.09795070e-05,   9.67626152e-06,
          8.64261251e-06,   7.79619058e-06,   7.08735215e-06],
       [  1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05,
          9.67626152e-06,   8.64261251e-06,   7.79619058e-06],
       [  1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05,
          1.09795070e-05,   9.67626152e-06,   8.64261251e-06],
       [  1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05,   1.09795070e-05,   9.67626152e-06],
       [  1.09795070e-05,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05],
       [  9.67626152e-06,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05],
       [  8.64261251e-06,   9.67626152e-06,   1.09795070e-05,
          1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05],
       [  7.79619058e-06,   8.64261251e-06,   9.67626152e-06,
          1.09795070e-05,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05],
       [  7.08735215e-06,   7.79619058e-06,   8.64261251e-06,
          9.67626152e-06,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05]])
}}}

{{{id=9|
#next, check that I get the same results when calling the function
mfullcheck = numpy.array(range((nmov+nfix)*(nmov+nfix)), dtype=float)
mfullcheck.shape=(nmov+nfix,nmov+nfix)
find_mutual_inductance(mfullcheck)
mfullcheck
///
array([[  3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05,   1.09795070e-05,   9.67626152e-06,
          8.64261251e-06,   7.79619058e-06,   7.08735215e-06],
       [  1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05,
          9.67626152e-06,   8.64261251e-06,   7.79619058e-06],
       [  1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05,
          1.09795070e-05,   9.67626152e-06,   8.64261251e-06],
       [  1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05,   1.09795070e-05,   9.67626152e-06],
       [  1.09795070e-05,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05],
       [  9.67626152e-06,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05],
       [  8.64261251e-06,   9.67626152e-06,   1.09795070e-05,
          1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05],
       [  7.79619058e-06,   8.64261251e-06,   9.67626152e-06,
          1.09795070e-05,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05],
       [  7.08735215e-06,   7.79619058e-06,   8.64261251e-06,
          9.67626152e-06,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05]])
}}}

{{{id=10|
#mm is now the mutual inductance matrix in microHenrys
#The calculation portions of the initialize procedure are now completely implemented
#correctly???





#Here are the arrays that will hold the data from the simulation
rstor = numpy.array(range(nmov*(ntim+1)), dtype=float)
rstor.shape = (nmov, ntim+1)

zstor = numpy.array(range(nmov*(ntim+1)), dtype=float)
zstor.shape = (nmov, ntim+1)

#At this point, I believe everything from the initialize function has been implemented

#make a copy of the reduced mutual inductance array
minv = mm.copy()

cntr = 0

ptime[0] = 0.0
coilI[0] = 0.0
bzero[0] = 0.0
zeroline[0] = 0.0
pctr = 0
heatenrg[0] = 0.0
enrgtot[0] = 0.0
work[0] = 0.0

#load the rstor and zstor arrays with the initial positions of the moving coils
for jj in range(0, nfix):
    rstor[jj,0] = r[jj+nfix]
    zstor[jj,0] = z[jj+nfix]
///
}}}

{{{id=26|
import numpy as np
#getting started on the compute_current function
def compute_current():
    global qq
    global adc
    global cur
    global coilI
    global cntr
    global denrg
    global zeroline
    global z
    global r
    global dheat
    global temp
    #first, we'll need some working arrays
    ccur = numpy.array(range(nmov+nfix), dtype=float)
    for i in range(0,nfix):
        ccur[i] = cur[0]
    for i in range(nfix,nfix+nmov):
        ccur[i] = cur[i-nfix+1]
    
    #variable names are copied from the original code
    nc1 = nmov+1
    nc11 = nmov
    idot = numpy.array(range(nc1), dtype=float)
 
    vv = numpy.array(range(nc1), dtype=float)
    
    #and here's the rather large array
    #that in fact never gets used for anything in the code, so I'm going to 
    #comment it out.  It does look like it would make a good debug tool though
    #pcur = numpy.array(range(nc1*ntim), dtype=float)
    #pcur.shape = (nc1, ntim)
    
    yp = numpy.array(range(ntim), dtype=float)
    tim = numpy.array(range(ntim), dtype=float)
    tim[0] = 0
    
    
    #add the external inductance into the reduced mutual inductance array
    for i in range(0, nmov+1):
        mm[i,i] = mm[i,i]+lckt[i]
    
    
    #NOTE:  Here's why we left a copy of the original mm hanging around in mmold
    dmmdt = (mm-mmold)/dt
    
    #make a copy of mm
    mident = mm.copy
    
    #now, invert mm
    minv = inv(mm)
    
    #make the inverse symmetric
    for i in range(0, nmov):
        for j in range(i+1, nmov+1):
            minv[i,j] = (minv[i,j]+minv[j,i])/2.0
            minv[j,i] = minv[i,j]
            
    #now compute the currents in the coils
    #the following is an operation on matrices, but appears not to be a matrix multiply
    #which would be denoted in IDL by #  I'm guessing it just multiplies element by element
    #eyer looks like a voltage
    eyer = cur*res
    print "Here's cur"
    print cur
    print "Here's res"
    print res
    #c = q/v so quec looks like a voltage
    quec = qq*adc
    print "Here's qq"
    print qq
    print "Here's adc"
    print adc
    #The following actually is a matrix multiply
    #it's cute because it's not l*di/dt
    idldt = np.dot(dmmdt, cur)
    
    #Make the arrays symmetric as usual
    for i in range(1, (nmov/2)+1):
        iii=nmov+1-i
        eyer[i]=(eyer[i]+eyer[iii])/2.0  
        eyer[iii]=eyer[i]
        quec[i]=(quec[i]+quec[iii])/2.0  
        quec[iii]=quec[i]
        idldt[i]=(idldt[i]+idldt[iii])/2.0  
        idldt[iii]=idldt[i]
        
    #voltage drop due to the admittance, (although, seems to be just cancelling the 
    #capacitance back out
    #minus the resistive voltage drop minus the inductive voltage drop
    #keep in mind that each argument is an array
    vv=quec-eyer-idldt
    print "Here's quec"
    print quec
    print "Here's eyer"
    print eyer
    print "Here's idldt"
    print idldt
    print "Here's vv"
    print vv
    print "Here's minv"
    print minv
    
    #what is this for?
    idot=np.dot(minv, vv)
    
    #Looks like we're making idot symmetric as well
    for i in range(1, (nmov/2)+1):
        iii=nmov+1-i
        idot[i]=(idot[i]+idot[iii])/2 
        idot[iii]=idot[i]
        
    cur=cur+idot*dt
    
    qq=qq-adc*cap*cur*dt
    coilI[cntr]=cur[cntr] 
    
    denrg=(cur[0]*qq[0]*dt)/cap[0]
    
    zeroline[cntr]=0.0
    
    #this variable is actually local!
    sbz = 0.0
    
    for j in range(0, nfix+nmov):
        zz = z[j]
        rr = 0.0
        rso = r[j]
        curr = 1e3*ccur[j]
        #implement the portion of dbcoil that returns bz
        sbz = sbz + dbcoilbzrzero(rso, zz, rr, curr)
        #sbz = sbz + bz
        
    #This is the field in kGauss
    bzero[cntr] = 10*sbz
    tms = 1e3*time
    
    #Calculate the heat
    dheat = 0.0
    for i in range(1, nmov+1):
        tt=temp[i]
        sig = specheat(tt)
        enrg = cur[i]^2*res[i]*dt
        dheat = dheat + enrg
        dtemp = dheat/(sig*mass)
        temp[i] = temp[i] + dtemp
        tt = temp[i]
        rho = rstvty(tt)
        #Cool!  Updating the resistance based on the temperature change
        #in the can
        res[i] = 1e3*rho*2*pi*rcan/(thickness*delta)
        
        
    
    
    #testing
    return cur
///
}}}

{{{id=25|
#now, let's test the above code by calling it at least
#now, find the mutual inductance
global mfull
mfull = find_mutual_inductance(mfull)
#then, reduce the mutual inductance array again
global mm
mm = make_reduced_matrix(mm, mfull)
my_current = compute_current()
///
}}}

{{{id=28|
#Interestingly, the current is already progression on each move through
cur
///
array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.])
}}}

{{{id=19|
delta
///
0.0026699999999999996
}}}

{{{id=20|
#make sure that the reduced matrix function works properly
mmtest = numpy.array(range(49), dtype=float)
mmtest.shape=(7,7)

mmtest = make_reduced_matrix(mm, mfull)
///
}}}

{{{id=12|
mmtest
///
array([[  1.07314394e+02,   2.78989073e+01,   2.36848067e+01,
          2.06557685e+01,   1.83188740e+01,   1.64388031e+01,
          1.48835427e+01],
       [  2.78989073e+01,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05,
          9.67626152e-06],
       [  2.36848067e+01,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05,
          1.09795070e-05],
       [  2.06557685e+01,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05],
       [  1.83188740e+01,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05],
       [  1.64388031e+01,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05],
       [  1.48835427e+01,   9.67626152e-06,   1.09795070e-05,
          1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05]])
}}}

{{{id=13|
#The simulation code lives in this cell
#currently, the full time range takes a while with ntim
#for kk in range(0,ntim):
for kk in range(0,9):
    #if the counter has advanced beyond nchange, then make the time step larger
    if cntr >= nchange:
        dt = ddt*10
    time = time + dt
    #store the current time in microseconds
    ptime[cntr] = time*1e3
    
    #Even those these funcitons have been called in initialize, it's important to call them even on 
    #the first loop through here.  Otherwise, mmold winds up with junk in it.
    #now, find the mutual inductance
    global mfull
    mfull = find_mutual_inductance(mfull)
    #then, reduce the mutual inductance array again
    global mm
    mm = make_reduced_matrix(mm, mfull)
    #now, finally, the first new simulation step, compute the currents
    compute_current()
///
WARNING: Output truncated!  
<html><a target='_new' href='/home/admin/4/cells/13/full_output.txt' class='file_link'>full_output.txt</a></html>



Here's cur
[ 0.  0.  0.  0.  0.  0.  0.]
Here's res
[  8.         35.8833313  35.8833313  35.8833313  35.8833313  35.8833313
  35.8833313]
Here's qq
[ 0.25  0.    0.    0.    0.    0.    0.  ]
Here's adc
[ 20000.      0.      0.      0.      0.      0.      0.]
Here's quec
[ 5000.     0.     0.     0.     0.     0.     0.]
Here's eyer
[ 0.  0.  0.  0.  0.  0.  0.]
Here's idldt
[ 0.  0.  0.  0.  0.  0.  0.]
Here's vv
[ 5000.     0.     0.     0.     0.     0.     0.]
Here's minv
[[  9.86474830e+00  -4.97616806e-01  -5.05186032e+00  -9.86290404e+00
   -9.86290404e+00  -5.05186032e+00  -4.97616806e-01]
 [ -4.97616806e-01   1.78591381e+01  -1.30888235e+01   1.21048367e+00
   -1.63987444e-01   1.60717265e-02  -3.66475355e-01]
 [ -5.05186032e+00  -1.30888235e+01   3.03968192e+01  -8.83661356e+00
    6.25144554e+00   2.11795318e+00   1.60717265e-02]
 [ -9.86290404e+00   1.21048367e+00  -8.83661356e+00   3.77600522e+01
   -4.12579942e+00   6.25144554e+00  -1.63987444e-01]
 [ -9.86290404e+00  -1.63987444e-01   6.25144554e+00  -4.12579942e+00
    3.77600522e+01  -8.83661356e+00   1.21048367e+00]
 [ -5.05186032e+00   1.60717265e-02   2.11795318e+00   6.25144554e+00
   -8.83661356e+00   3.03968192e+01  -1.30888235e+01]
 [ -4.97616806e-01  -3.66475355e-01   1.60717265e-02  -1.63987444e-01
    1.21048367e+00  -1.30888235e+01   1.78591381e+01]]
array([ 0.98647483, -0.04976168, -0.50518603, -0.9862904 , -0.9862904 ,
       -0.50518603, -0.04976168])
Here's cur
[ 0.98647483 -0.04976168 -0.50518603 -0.9862904  -0.9862904  -0.50518603
 -0.04976168]
Here's res
[  8.          35.88333176  35.88337899  35.88355901  35.88373904
  35.88378627  35.88378673]
Here's qq
[ 0.24998027  0.          0.          0.          0.          0.          0.        ]
Here's adc
[ 20000.      0.      0.      0.      0.      0.      0.]
Here's quec
[ 4999.60541007     0.             0.             0.             0.             0.
     0.        ]
Here's eyer
[  7.89179864  -1.78562621 -18.12788471 -35.39169869 -35.39169869
 -18.12788471  -1.78562621]
Here's idldt
[ 0.  0.  0.  0.  0.  0.  0.]
Here's vv
[  4.99171361e+03   1.78562621e+00   1.81278847e+01   3.53916987e+01
   3.53916987e+01   1.81278847e+01   1.78562621e+00]
Here's minv
[[  9.86474830e+00  -4.97616806e-01  -5.05186032e+00  -9.86290404e+00
   -9.86290404e+00  -5.05186032e+00  -4.97616806e-01]
 [ -4.97616806e-01   1.78591381e+01  -1.30888235e+01   1.21048367e+00

...

Here's idldt
[ 0.  0.  0.  0.  0.  0.  0.]
Here's vv
[ 4937.27110799    14.84272923   119.46341169   230.14030276   230.14030276
   119.46341169    14.84272923]
Here's minv
[[  9.86474830e+00  -4.97616806e-01  -5.05186032e+00  -9.86290404e+00
   -9.86290404e+00  -5.05186032e+00  -4.97616806e-01]
 [ -4.97616806e-01   1.78591381e+01  -1.30888235e+01   1.21048367e+00
   -1.63987444e-01   1.60717265e-02  -3.66475355e-01]
 [ -5.05186032e+00  -1.30888235e+01   3.03968192e+01  -8.83661356e+00
    6.25144554e+00   2.11795318e+00   1.60717265e-02]
 [ -9.86290404e+00   1.21048367e+00  -8.83661356e+00   3.77600522e+01
   -4.12579942e+00   6.25144554e+00  -1.63987444e-01]
 [ -9.86290404e+00  -1.63987444e-01   6.25144554e+00  -4.12579942e+00
    3.77600522e+01  -8.83661356e+00   1.21048367e+00]
 [ -5.05186032e+00   1.60717265e-02   2.11795318e+00   6.25144554e+00
   -8.83661356e+00   3.03968192e+01  -1.30888235e+01]
 [ -4.97616806e-01  -3.66475355e-01   1.60717265e-02  -1.63987444e-01
    1.21048367e+00  -1.30888235e+01   1.78591381e+01]]
array([ 7.37166794, -0.48367082, -3.76322907, -7.23142955, -7.23142955,
       -3.76322907, -0.48367082])
Here's cur
[ 7.37166794 -0.48367082 -3.76322907 -7.23142955 -7.23142955 -3.76322907
 -0.48367082]
Here's res
[  8.          35.88345878  35.89210175  35.92427002  35.95643812
  35.96508101  35.9652085 ]
Here's qq
[ 0.24932124  0.          0.          0.          0.          0.          0.        ]
Here's adc
[ 20000.      0.      0.      0.      0.      0.      0.]
Here's quec
[ 4986.42483452     0.             0.             0.             0.             0.
     0.        ]
Here's eyer
[  58.97334352  -17.37555188 -135.20751942 -259.90013857 -259.90013857
 -135.20751942  -17.37555188]
Here's idldt
[ 0.  0.  0.  0.  0.  0.  0.]
Here's vv
[ 4927.45149099    17.37555188   135.20751942   259.90013857   259.90013857
   135.20751942    17.37555188]
Here's minv
[[  9.86474830e+00  -4.97616806e-01  -5.05186032e+00  -9.86290404e+00
   -9.86290404e+00  -5.05186032e+00  -4.97616806e-01]
 [ -4.97616806e-01   1.78591381e+01  -1.30888235e+01   1.21048367e+00
   -1.63987444e-01   1.60717265e-02  -3.66475355e-01]
 [ -5.05186032e+00  -1.30888235e+01   3.03968192e+01  -8.83661356e+00
    6.25144554e+00   2.11795318e+00   1.60717265e-02]
 [ -9.86290404e+00   1.21048367e+00  -8.83661356e+00   3.77600522e+01
   -4.12579942e+00   6.25144554e+00  -1.63987444e-01]
 [ -9.86290404e+00  -1.63987444e-01   6.25144554e+00  -4.12579942e+00
    3.77600522e+01  -8.83661356e+00   1.21048367e+00]
 [ -5.05186032e+00   1.60717265e-02   2.11795318e+00   6.25144554e+00
   -8.83661356e+00   3.03968192e+01  -1.30888235e+01]
 [ -4.97616806e-01  -3.66475355e-01   1.60717265e-02  -1.63987444e-01
    1.21048367e+00  -1.30888235e+01   1.78591381e+01]]
array([ 8.21362667, -0.55654258, -4.1911408 , -8.03520525, -8.03520525,
       -4.1911408 , -0.55654258])
}}}

{{{id=27|
idot
///
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "_sage_input_97.py", line 10, in <module>
    exec compile(u'open("___code___.py","w").write("# -*- coding: utf-8 -*-\\n" + _support_.preparse_worksheet_cell(base64.b64decode("aWRvdA=="),globals())+"\\n"); execfile(os.path.abspath("___code___.py"))' + '\n', '', 'single')
  File "", line 1, in <module>
    
  File "/tmp/tmpyGVEA_/___code___.py", line 2, in <module>
    exec compile(u'idot' + '\n', '', 'single')
  File "", line 1, in <module>
    
NameError: name 'idot' is not defined
}}}

{{{id=22|
mm
///
array([[  1.07364394e+02,   2.78989073e+01,   2.36848067e+01,
          2.06557685e+01,   1.83188740e+01,   1.64388031e+01,
          1.48835427e+01],
       [  2.78989073e+01,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05,
          9.67626152e-06],
       [  2.36848067e+01,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05,
          1.09795070e-05],
       [  2.06557685e+01,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05],
       [  1.83188740e+01,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05],
       [  1.64388031e+01,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05],
       [  1.48835427e+01,   9.67626152e-06,   1.09795070e-05,
          1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05]])
}}}

{{{id=21|
mfull[2:9, 2:9]
///
array([[  3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05,   1.09795070e-05,   9.67626152e-06,
          8.64261251e-06],
       [  1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05,
          9.67626152e-06],
       [  1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05,
          1.09795070e-05],
       [  1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05],
       [  1.09795070e-05,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05],
       [  9.67626152e-06,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05],
       [  8.64261251e-06,   9.67626152e-06,   1.09795070e-05,
          1.27052998e-05,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05]])
}}}

{{{id=23|

///
}}}
