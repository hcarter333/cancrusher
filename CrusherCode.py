
{{{id=1|
import numpy as np
from numpy.linalg import inv
rstvty(temperature)=(4.5+1.133*10^-2*temperature)*10^-8
#setup the constants
perm0 = 1.25663706*1e-6
nmov = 6
nfix = 3

#need to change this for Pb
density = 2824.0

thickness = .00016

#Setup an array for coil dimensions
r = numpy.array([2, 2, 2, 2, 2, 2, 2, 2, 2], dtype=float)
#Here's the array for coil z coordinates
z = numpy.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], dtype=float)
#Here's the array for coil-can separtion
dr = numpy.array([.01, .01, .01, .01, .01, .01, .01, .01, .01], dtype=float)

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
dbcoilflux(rcoil, zc, rc, curren) = (2*perm0/(2*pi))*curren*((rcoil*rc)/(argm(rcoil, zc, rc)/t(rcoil, zc, rc)^0.5))*((2-argm(rcoil, zc, rc))*elliptic_kc(argm(rcoil, zc, rc))-(2*elliptic_ec(argm(rcoil, zc, rc))))
///
}}}

{{{id=11|

///
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

global mfull
mfull = find_mutual_inductance(mfull)
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
#now take mfull created above and use it to generate the reduced mutual induction matrix mm
#mm is dimensioned for only the moving coils plus one extra row and column
#ultimately we'll need to track the immediately previous version of mm
#in the original code this is called mmold, and here too!


#also create a function here that encapsulates the reduced array creation
#since this will need to be called multiple times
def make_reduced_matrix(mreduced, mfullwork):
    global mmold
    mmold = mm.copy()
    jmwork = numpy.array(range(nmov), dtype=float)
    #nt contains the fixed coil portion of the mutual inductance array
    #it's actually just a 
    ntwork = mfull[0:(nfix-1), 0:(nfix-1)].copy()
    nttwork = sum(ntwork)
    
    #jt is used to total the mutual inductance associated with each moveable coil
    #the result is placed in the jmwork array
    for i in range(0, 6):
        jtwork=mfull[i+nfix,0:nfix-1]
        jmwork[i]=sum(jtwork)
    
    #Here's why mm was dimensioned one bigger than it should have been
    #I don't have a good reason for this yet other than it came out of the old code
    #notice, however that we do overwrite the 0 row and column next
    mreduced = mfullwork[nfix-1:nmov+nfix, nfix-1:nmov+nfix].copy()
    #mm = mm*1e6

    for i in range(1, nmov+1):
        mreduced[i,0] = 1e6*jmwork[i-1]
        mreduced[0,i] = 1e6*jmwork[i-1]
    
    mreduced[0,0] = 1e6*nttwork

    return mreduced


mm = make_reduced_matrix(mm, mfull)
#mm is now the mutual inductance matrix in microHenrys
#The calculation portions of the initialize procedure are now completely implemented
#correctly???

delta = z[nmov+1] - z[nmov+2]


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

#setup the admittance matrix.  I've made this a single statement here, keep in mind that 
#there's a there is a loop treatment in the original code
adc = numpy.array([1/cap[0], 0, 0, 0, 0, 0, 0], dtype=float)
qq = numpy.array([V0[0]*cap[0], 0, 0, 0, 0, 0, 0], dtype=float)

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
bzeroline[0] = 0.0
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

{{{id=24|
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
    #first, we'll need some working arrays
    ccur = numpy.array(range(nmov+nfix), dtype=float)
    for i in range(0,nfix):
        ccur[i] = cur[0]
    for i in range(nfix,nfix+nmov):
        ccur[i] = cur[i-nfix+1]
    
    #variable names are copied from the original code
    nc1 = nmov+1
    nc11 = nmov
    idot = numpy.array(range(nc1*2), dtype=float)
    idot.shape = (nc1, 2)

    vv = numpy.array(range(nc1*2), dtype=float)
    vv.shape = (nc1, 2)
    
    #and here's the rather large array
    pcur = numpy.array(range(nc1*ntim), dtype=float)
    pcur.shape = (nc1, ntim)
    
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
    #c = q/v so quec looks like a voltage
    quec = qq*adc
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
        
        #sbz = sbz + bz
        
    #testing
    return cur
///
}}}

{{{id=14|
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

{{{id=19|
my_current
///
array([ -7.34696796e-09,   2.38495513e-03,   1.02871591e-03,
         7.62842744e-04,   7.62842744e-04,   1.02871591e-03,
         2.38495513e-03])
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
for kk in range(0,3):
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
    #
///
}}}

{{{id=22|
mm
///
array([[  1.07314394e+02,   2.78989073e+01,   2.36848067e+01,
          2.06557685e+01,   3.00000000e+06,   4.00000000e+06,
          5.00000000e+06],
       [  2.78989073e+01,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05,   1.27052998e-05,   1.09795070e-05,
          9.67626152e-06],
       [  2.36848067e+01,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05,   1.51936075e-05,   1.27052998e-05,
          1.09795070e-05],
       [  2.06557685e+01,   1.51936075e-05,   1.95128445e-05,
          3.41443527e-05,   1.95128445e-05,   1.51936075e-05,
          1.27052998e-05],
       [  3.00000000e+06,   1.27052998e-05,   1.51936075e-05,
          1.95128445e-05,   3.41443527e-05,   1.95128445e-05,
          1.51936075e-05],
       [  4.00000000e+06,   1.09795070e-05,   1.27052998e-05,
          1.51936075e-05,   1.95128445e-05,   3.41443527e-05,
          1.95128445e-05],
       [  5.00000000e+06,   9.67626152e-06,   1.09795070e-05,
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
