
{{{id=1|
import numpy
rstvty(temperature)=(4.5+1.133*10^-2*temperature)*10^-8
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
#In this cell, setup the flux calculation that mimics what goes in in dbcoil
#call it the same name as the original code for readability
#leave the unneccasary parts out though
#setup the constants
perm0 = 1.25663706*1e-6
t(rcoil, zc, rc) = (rc+rcoil)^2+zc^2
argm(rcoil, zc, rc)= 4*rcoil*rc/t(rcoil,zc,rc)
dbcoilflux(rcoil, zc, rc, curren) = (2*perm0/(2*pi))*curren*((rcoil*rc)/(argm(rcoil, zc, rc)/t(rcoil, zc, rc)^0.5))*((2-argm(rcoil, zc, rc))*elliptic_kc(argm(rcoil, zc, rc))-(2*elliptic_ec(argm(rcoil, zc, rc))))
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

#The mutual inductance array has to be calculated multiple times, so I'm putting 
#it into a Python funcion, (I hope)
def find_mutual_inductance(mfullarray):
    #first for the off-diagonal elements
    for i in range(0, 8):
       for j in range(i+1, 9):
            mfullarray[i,j] = dbcoilflux(r[i], z[j]-z[i], r[j]-dr[j], 1)
    
    #next for the diagonal elements
    for i in range(0, 9):
        mfullarray[i,i] = dbcoilflux(r[i], z[i]-z[i], r[i]-dr[i], 1)
   
    #finally copy the off diagonal elements to the other side of the array
    for i in range(1, 9):
        for j in range(0, i):
            mfullarray[i,j] = mfull[j,i]
            
    return mfullarray

#setting up the mutual inductance raw flux array
mfull = numpy.array(range(81), dtype=float)
mfull.shape=(9,9)

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
mfullcheck = numpy.array(range(81), dtype=float)
mfullcheck.shape=(9,9)
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

mm = numpy.array(range(49), dtype=float)
mm.shape=(7,7)
mmold = mm.copy()

#also create a function here that encapsulates the reduced array creation
#since this will need to be called multiple times
def make_reduced_matrix(mreduced, mfullwork):
    jmwork = numpy.array(range(6), dtype=float)
    #nt contains the fixed coil portion of the mutual inductance array
    #it's actually just a 
    ntwork = mfull[0:3-1, 0:3-1].copy()
    nttwork = sum(ntwork)
    
    #jt is used to total the mutual inductance associated with each moveable coil
    #the result is placed in the jmwork array
    for i in range(0, 6):
        jtwork=mfull[i+3,0:3-1]
        jmwork[i]=sum(jtwork)
    
    #Here's why mm was dimensioned one bigger than it should have been
    #I don't have a good reason for this yet other than it came out of the old code
    #notice, however that we do overwrite the 0 row and column next
    mreduced = mfullwork[2:9, 2:9].copy()
    #mm = mm*1e6

    for i in range(1, 7):
        mreduced[i,0] = 1e6*jmwork[i-1]
        mreduced[0,i] = 1e6*jmwork[i-1]
    
    mreduced[0,0] = 1e6*nttwork

    return mreduced


mm = make_reduced_matrix(mm, mfull)
#mm is now the mutual inductance matrix in microHenrys
#The calculation portions of the initialize procedure are now completely implemented
#correctly???

#in the following section, we'll need the resistances and capacitances that were pulled 
#from the input data file
#We'll also need a number of peripheral variables
#need to change this for Pb
density = 2824.0

delta = z[3+1] - z[3+2]

thickness = .00016

#Here's a clue to how the radius positions of coil and can should be worked out
#rcan = r[nfix+1]
#the indexing above seems odd
rcan = r[3+1]
mass = 2*pi*rcan*density*thickness*delta

#array stores the temerature of the can under each coil
temp = numpy.array([293.0, 293.0, 293.0, 293.0, 293.0, 293.0, 293.0], dtype=float)

#The following are the resistance, capacitance, inductance, admittance, 
#current, and charge arrays with one entry per can modelling coil
#They have been filled in here with the values from the model data file
#In the model data file, there are 7 values as opposed to the six I expected
#I neeed todo check on this
res = numpy.array([8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)

cap = numpy.array([0.00005, -1, -1, -1, -1, -1, -1, -1], dtype=float)

lckt = numpy.array([0.05, 0, 0, 0, 0, 0, 0], dtype=float)

cur = numpy.array([0.0, 0, 0, 0, 0, 0, 0], dtype=float)

V0 = numpy.array([5000.0, 0, 0, 0, 0, 0, 0], dtype=float)

#the next step is to fill in the res array with calls to rstvty
#This begins to make more sense after a fahion.  The one bigger is intentional, the 
#first value in the data file was real, all the others are filled in
for i in range(1, 7):
    res[i] = 1e3*2*pi*rcan*rstvty(temp[i])/(thickness*delta)

#now setup the variables for stepping through the simulation
ntim = 600
ddt = 0.00002
nchange = 20

#setup the admittance matrix.  I've made this a single statement here, keep in mind that 
#there's a there is a loop treatment in the original code
adc = numpy.array([1/cap[0], 0, 0, 0, 0, 0, 0], dtype=float)
qq = numpy.array([V0[0]*cap[0], 0, 0, 0, 0, 0, 0], dtype=float)

#Here are the arrays that will hold the data from the simulation
rstor = numpy.array(range(6*(ntim+1)), dtype=float)
rstor.shape = (6, ntim+1)

zstor = numpy.array(range(6*(ntim+1)), dtype=float)
zstor.shape = (6, ntim+1)

coilI = numpy.array(range(ntim+1), dtype=float)
bzero = numpy.array(range(ntim+1), dtype=float)
bzeroline = numpy.array(range(ntim+1), dtype=float)
heatenrg = numpy.array(range(ntim+1), dtype=float)
work = numpy.array(range(ntim+1), dtype=float)
enrgtot = numpy.array(range(ntim+1), dtype=float)
ptime = numpy.array(range(ntim+1), dtype=float)

dt = ddt
time = 0.0
aaa = 1.0
ntim1 = ntim-1
#This concludes all the setup of variables that is done by the initialize function
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
for jj in range(0, 3):
    rstor[jj,0] = r[jj+3]
    zstor[jj,0] = z[jj+3]
///
}}}

{{{id=14|

///
}}}

{{{id=19|
mm
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

{{{id=20|
#make sure that the reduced matrix function works properly
mmtest = numpy.array(range(49), dtype=float)
mmtest.shape=(7,7)

mmtest = make_reduced_matrix(mmtest, mfull)
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
    
    #now, find the mutual inductance
    mfull = find_mutual_inductance(mfull)
    #then, reduce the mutual inductance array again
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
