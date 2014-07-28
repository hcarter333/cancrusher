
{{{id=34|
import numpy
from numpy.linalg import inv
global mfull

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
#vr is the radial velocity of the can?  For some reason it included enough entries 
#for the can and the fixed coil
#vz is the z velocity, but it looks unused.  It looks like the axial force 
#implementation was left incomplete in the initial version
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

#Just changed the capacitance to 500 from 50
cap = numpy.array([0.0005, -1, -1, -1, -1, -1, -1], dtype=float)

lckt = numpy.array([0.05, 0, 0, 0, 0, 0, 0], dtype=float)

cur = numpy.array([0.0, 0, 0, 0, 0, 0, 0], dtype=float)
ccur = numpy.array(range(nmov+nfix), dtype=float)

#Changed to 3000 from 5000
V0 = numpy.array([3000.0, 0, 0, 0, 0, 0, 0], dtype=float)

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
coilOutTime = numpy.array(range((ntim+1)*2), dtype=float)
coilOutTime.shape = (ntim+1, 2)
bzero = numpy.array(range(ntim+1), dtype=float)
zeroline = numpy.array(range(ntim+1), dtype=float)
heatenrg = numpy.array(range(ntim+1), dtype=float)
work = numpy.array(range(ntim+1), dtype=float)
enrgtot = numpy.array(range(ntim+1), dtype=float)

#store the current time in microseconds
ptime = numpy.array(range(ntim+1), dtype=float)
denrg = 0.0
dheat = 0.0
dwork = 0.0

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
#this was missing the curren term next to perm0
#Also, the term containing perm0 was supposed to be multiplied, not added
dbcoilbz(rcoil, zc, rc, curren) = (bzfac(rcoil, zc, rc)*elliptic_ec(argm(rcoil, zc, rc))+elliptic_kc(argm(rcoil, zc, rc)))*((2*perm0*curren/(2*pi))/(t(rcoil, zc, rc)^0.5))

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

#The above was to provide the magnetic flux at one coil due to another for finding 
#the mutual inductance of a set of two coils.  The following models the coils and 
#calls the above functions to fill in the mutual induction matrix
#Utilize a python for loop here and I guess go ahead and use sage matrices, although 
#I'm not sure there's not a reason to use numpy matrices or plain old Python arrays
#The documentation for numpy arrays is the easiest for me to read, so I'm going that 
#direction.  I went ahead and put the import at the top of the worksheet

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
    
#Here are the arrays that will hold the data from the simulation
rstor = numpy.array(range(nmov*(ntim+1)), dtype=float)
rstor.shape = (nmov, ntim+1)

zstor = numpy.array(range(nmov*(ntim+1)), dtype=float)
zstor.shape = (nmov, ntim+1)

cntr = 0

ptime[0] = 0.0
coilI[0] = 0.0
bzero[0] = 0.0
zeroline[0] = 0.0
pctr = 0
heatenrg[0] = 0.0
enrgtot[0] = 0.0
work[0] = 0.0

    
    
#calls to initialize the code.  This should eventually be put into an initialize funciton
mfull = find_mutual_inductance(mfull)

mm = make_reduced_matrix(mm, mfull)

#load the rstor and zstor arrays with the initial positions of the moving coils
for jj in range(0, nfix):
    rstor[jj,0] = r[jj+nfix]
    zstor[jj,0] = z[jj+nfix]
///
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

{{{id=49|
#and inside the coil, being careful to avoid the point r = 0
plot(dbcoilflux(2, 0, rc, 1), 0.01, 2)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=39|
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
    global ptime
    global ccur
    #first, we'll need some working arrays
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
    #print "Here's cur"
    #print cur
    #print "Here's res"
    #print res
    #c = q/v so quec looks like a voltage
    quec = qq*adc
    #print "Here's qq"
    #print qq
    #print "Here's adc"
    #print adc
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
    #print "Here's quec"
    #print quec
    #print "Here's eyer"
    #print eyer
    #print "Here's idldt"
    #print idldt
    #print "Here's vv"
    #print vv
    #print "Here's minv"
    #print minv
    
    #what is this for?
    idot=np.dot(minv, vv)
    
    #Looks like we're making idot symmetric as well
    for i in range(1, (nmov/2)+1):
        iii=nmov+1-i
        idot[i]=(idot[i]+idot[iii])/2 
        idot[iii]=idot[i]
        
    cur=cur+idot*dt
    
    qq=qq-adc*cap*cur*dt
    coilI[cntr]=cur[0] 
    coilOutTime[cntr,0]= ptime[cntr]
    coilOutTime[cntr,1] = cur[0]
    
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
    #return cur
///
}}}

{{{id=26|
#This is the move_can function
#It returns the radial velocity of the can

def move_can():
    global z
    global r
    global ccur
    global dt
    global dwork
    
    brg = numpy.array(range(nfix + nmov), dtype=float)
    bzt = numpy.array(range(nfix + nmov), dtype=float)
    
    for i in range(1,nmov+1):
        sbz = 0.0
        imov = nfix + i - 1
        for j in range(0,nfix+nmov):
            if j==imov:
                 break
            zz = z[imov]-z[j]
            rr = r[imov]
            rso = r[j]
            curr = 1e3*ccur[j]
            #if rr is zero, the can has broken in half and all bets are off anyway           
            #if rr == 0:
            #    sbz = sbz + dbcoilbzrzero(rso, zz, rr, curr)
            #else:
            sbz = sbz + dbcoilbz(rso, zz, rr, curr)
                
        bzt[i+nfix-1] = sbz
    dwork = 0.0
    print "nmov"
    print nmov
    for i in range(1,(nmov/2)+1):
        #This looks like we're hitting opposite edges of the moving coils, the can, and 
        #working towards the center
        #print "i"
        #print i
        ii=nfix+i-1
        iii = (nfix + nmov) - i
        
        #Find force in kNewtons
        forcer = bzt[ii]*ccur[ii]*2*pi*r[ii]
        
        #Get the differential in velocity from the force and mass
        dvr = forcer*dt/mass
        vrnew = vr[ii]+dvr
        
        #get the new r position using the velocity
        rnew=r[ii]+vr[ii]*1e-3*dt
        print "rii"
        print r[ii]
        print "vrii"
        print vr[ii]
        print"dt"
        print dt
        
        #Find work in Joules
        dwork=dwork+2*forcer*vrnew*dt
        vr[ii]=vrnew
        r[ii] = rnew
        vr[iii]=vrnew
        r[iii]=rnew
    #print "completed loop"
///
}}}

{{{id=13|
#The simulation code lives in this cell
#currently, the full time range takes a while with ntim
#for kk in range(0,ntim):
for kk in range(0,599):
    #if the counter has advanced beyond nchange, then make the time step larger
    if cntr >= nchange:
        dt = ddt*10
    print cntr
    cntr = cntr + 1
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
    print "Time in microseconds"
    print ptime[cntr]
    compute_current()
    move_can()
    
    #track the heat and work for this step
    heatenrg[cntr]=heatenrg[cntr+1]+dheat
    work[cntr]=work[cntr-1]+dwork
    enrgtot[cntr]=enrgtot[cntr-1]+denrg
    for jj in range(0,nmov):
        jjmov = jj + nfix
        rstor[jj,kk+1] = r[jjmov]
        zstor[jj,kk+1] = z[jjmov]
///
WARNING: Output truncated!  
<html><a target='_new' href='/home/admin/4/cells/13/full_output.txt' class='file_link'>full_output.txt</a></html>



0
Time in microseconds
0.02
nmov
6
rii
0.0312
vrii
0.3
dt
0.0000200000000000000
rii
0.0312
vrii
0.4
dt
0.0000200000000000000
rii
0.0312
vrii
0.5
dt
0.0000200000000000000
1
Time in microseconds
0.04
nmov
6
rii
0.031200006
vrii
0.3
dt
0.0000200000000000000
rii
0.031200008
vrii
0.4
dt
0.0000200000000000000
rii
0.03120001
vrii
0.5
dt
0.0000200000000000000
2
Time in microseconds
0.06
nmov
6
rii
0.031200012
vrii
0.29996002604
dt
0.0000200000000000000
rii
0.031200016

...

-33675.5333951
dt
0.000200000000000000
rii
0.000740212325234
vrii
-3189717.40105
dt
0.000200000000000000
367
Time in microseconds
70.0
nmov
6
rii
0.0140248435912
vrii
-1.23570418117e+19
dt
0.000200000000000000
rii
-0.000221227583424
vrii
-3.10283778144e+20
dt
0.000200000000000000
rii
-0.637203267885
vrii
-2.97154711146e+22
dt
0.000200000000000000
368
Traceback (most recent call last):            dt = ddt*10
  File "", line 1, in <module>
    
  File "/tmp/tmpGBU5cS/___code___.py", line 5, in <module>
    exec compile(u'for kk in range(_sage_const_0 ,_sage_const_599 ):\n    #if the counter has advanced beyond nchange, then make the time step larger\n    if cntr >= nchange:\n        dt = ddt*_sage_const_10 \n    print cntr\n    cntr = cntr + _sage_const_1 \n    time = time + dt\n    #store the current time in microseconds\n    ptime[cntr] = time*_sage_const_1e3 \n    \n    #Even those these funcitons have been called in initialize, it\'s important to call them even on \n    #the first loop through here.  Otherwise, mmold winds up with junk in it.\n    #now, find the mutual inductance\n    global mfull\n    mfull = find_mutual_inductance(mfull)\n    #then, reduce the mutual inductance array again\n    global mm\n    mm = make_reduced_matrix(mm, mfull)\n    #now, finally, the first new simulation step, compute the currents\n    print "Time in microseconds"\n    print ptime[cntr]\n    compute_current()\n    move_can()\n    \n    #track the heat and work for this step\n    heatenrg[cntr]=heatenrg[cntr+_sage_const_1 ]+dheat\n    work[cntr]=work[cntr-_sage_const_1 ]+dwork\n    enrgtot[cntr]=enrgtot[cntr-_sage_const_1 ]+denrg\n    for jj in range(_sage_const_0 ,nmov):\n        jjmov = jj + nfix\n        rstor[jj,kk+_sage_const_1 ] = r[jjmov]\n        zstor[jj,kk+_sage_const_1 ] = z[jjmov]' + '\n', '', 'single')
  File "", line 15, in <module>
    
  File "/tmp/tmpltno8K/___code___.py", line 144, in find_mutual_inductance
    mfullarray[i,j] = dbcoilflux(r[i], z[j]-z[i], r[j]-dr[j], _sage_const_1 )
  File "expression.pyx", line 4361, in sage.symbolic.expression.Expression.__call__ (sage/symbolic/expression.cpp:21618)
  File "/home/sage/sage-6.2/local/lib/python2.7/site-packages/sage/symbolic/callable.py", line 477, in _call_element_
    return SR(_the_element.substitute(**d))
  File "expression.pyx", line 4212, in sage.symbolic.expression.Expression.substitute (sage/symbolic/expression.cpp:20868)
  File "/home/sage/sage-6.2/local/lib/python2.7/site-packages/sage/functions/special.py", line 505, in _eval_
    if self.name() in repr(s):
  File "/home/sage/sage-6.2/local/lib/python2.7/site-packages/sage/interfaces/maxima_abstract.py", line 1427, in __repr__
    r = P.get(self._name)
  File "/home/sage/sage-6.2/local/lib/python2.7/site-packages/sage/interfaces/maxima_lib.py", line 525, in get
    s = self.eval('%s;'%var)
  File "/home/sage/sage-6.2/local/lib/python2.7/site-packages/sage/interfaces/maxima_lib.py", line 423, in _eval_line
    if statement: result = ((result + '\n') if result else '') + max_to_string(maxima_eval("#$%s$"%statement))
  File "/home/sage/sage-6.2/local/lib/python2.7/site-packages/sage/interfaces/maxima_lib.py", line 260, in max_to_string
    return maxprint(s).python()[1:-1]
  File "ecl.pyx", line 784, in sage.libs.ecl.EclObject.__call__ (sage/libs/ecl.c:6640)
  File "ecl.pyx", line 355, in sage.libs.ecl.ecl_safe_apply (sage/libs/ecl.c:4528)
  File "c_lib.pyx", line 89, in sage.ext.c_lib.sig_raise_exception (sage/ext/c_lib.c:1094)
FloatingPointError: Floating point exception
}}}

{{{id=38|
r
///
array([  3.43000000e-02,   3.43000000e-02,   3.43000000e-02,
        -4.09642510e+19,   3.00055899e+19,  -1.65575866e+21,
        -1.65575866e+21,   3.00055899e+19,  -4.09642510e+19])
}}}

{{{id=35|
list_plot(coilOutTime[0:365, 0:2])
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=36|
r
///
array([  3.43000000e-02,   3.43000000e-02,   3.43000000e-02,
        -4.09221476e+19,   2.99745192e+19,  -1.65404165e+21,
        -1.65404165e+21,   2.99745192e+19,  -4.09221476e+19])
}}}

{{{id=40|
rstor[nfix, 1]
///
0.03120001
}}}

{{{id=41|
rstor[0:nmov, 1]
///
array([ 0.03120001,  0.03120001,  0.03120001,  0.03120001,  0.03120001,
        0.03120001])
}}}

{{{id=42|
show(list_plot(rstor[0:nmov, 1], plotjoined=True) + list_plot(rstor[0:nmov, 77], plotjoined=True) + list_plot(rstor[0:nmov, 138], plotjoined=True) + list_plot(rstor[0:nmov, 198], plotjoined=True) + list_plot(rstor[0:nmov, 258], plotjoined=True) + list_plot(rstor[0:nmov, 318], plotjoined=True))
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=43|
rstor[0:nmov, 3]
///
array([ 0.03120002,  0.03120002,  0.03120002,  0.03120002,  0.03120002,
        0.03120002])
}}}

{{{id=44|
rstor[0:nmov, 10]
///
array([ 0.03120003,  0.03119957,  0.0311987 ,  0.0311987 ,  0.03119957,
        0.03120003])
}}}

{{{id=46|
rstor[0:nmov, 100]
///
array([ 0.01660943,  0.00882953,  0.00515049,  0.00515049,  0.00882953,
        0.01660943])
}}}

{{{id=47|
rstor[0:nmov, 119]
///
array([  119.,   720.,  1321.,  1922.,  2523.,  3124.])
}}}

{{{id=48|
ptime
///
WARNING: Output truncated!  
<html><a target='_new' href='/home/admin/4/cells/48/full_output.txt' class='file_link'>full_output.txt</a></html>



array([  0.00000000e+00,   2.00000000e-02,   4.00000000e-02,
         6.00000000e-02,   8.00000000e-02,   1.00000000e-01,
         1.20000000e-01,   1.40000000e-01,   1.60000000e-01,
         1.80000000e-01,   2.00000000e-01,   2.20000000e-01,
         2.40000000e-01,   2.60000000e-01,   2.80000000e-01,
         3.00000000e-01,   3.20000000e-01,   3.40000000e-01,
         3.60000000e-01,   3.80000000e-01,   4.00000000e-01,
         6.00000000e-01,   8.00000000e-01,   1.00000000e+00,
         1.20000000e+00,   1.40000000e+00,   1.60000000e+00,
         1.80000000e+00,   2.00000000e+00,   2.20000000e+00,
         2.40000000e+00,   2.60000000e+00,   2.80000000e+00,
         3.00000000e+00,   3.20000000e+00,   3.40000000e+00,
         3.60000000e+00,   3.80000000e+00,   4.00000000e+00,
         4.20000000e+00,   4.40000000e+00,   4.60000000e+00,
         4.80000000e+00,   5.00000000e+00,   5.20000000e+00,
         5.40000000e+00,   5.60000000e+00,   5.80000000e+00,
         6.00000000e+00,   6.20000000e+00,   6.40000000e+00,
         6.60000000e+00,   6.80000000e+00,   7.00000000e+00,
         7.20000000e+00,   7.40000000e+00,   7.60000000e+00,
         7.80000000e+00,   8.00000000e+00,   8.20000000e+00,
         8.40000000e+00,   8.60000000e+00,   8.80000000e+00,
         9.00000000e+00,   9.20000000e+00,   9.40000000e+00,
         9.60000000e+00,   9.80000000e+00,   1.00000000e+01,
         1.02000000e+01,   1.04000000e+01,   1.06000000e+01,
         1.08000000e+01,   1.10000000e+01,   1.12000000e+01,
         1.14000000e+01,   1.16000000e+01,   1.18000000e+01,
         1.20000000e+01,   1.22000000e+01,   1.24000000e+01,
         1.26000000e+01,   1.28000000e+01,   1.30000000e+01,
         1.32000000e+01,   1.34000000e+01,   1.36000000e+01,
         1.38000000e+01,   1.40000000e+01,   1.42000000e+01,
         1.44000000e+01,   1.46000000e+01,   1.48000000e+01,
         1.50000000e+01,   1.52000000e+01,   1.54000000e+01,
         1.56000000e+01,   1.58000000e+01,   1.60000000e+01,
         1.62000000e+01,   1.64000000e+01,   1.66000000e+01,
         1.68000000e+01,   1.70000000e+01,   1.72000000e+01,
         1.74000000e+01,   1.76000000e+01,   1.78000000e+01,
         1.80000000e+01,   1.82000000e+01,   1.84000000e+01,
         1.86000000e+01,   1.88000000e+01,   1.90000000e+01,
         1.92000000e+01,   1.94000000e+01,   1.96000000e+01,
         1.98000000e+01,   2.00000000e+01,   2.02000000e+01,
         2.04000000e+01,   2.06000000e+01,   2.08000000e+01,
         2.10000000e+01,   2.12000000e+01,   2.14000000e+01,
         2.16000000e+01,   2.18000000e+01,   2.20000000e+01,
         2.22000000e+01,   2.24000000e+01,   2.26000000e+01,
         2.28000000e+01,   2.30000000e+01,   2.32000000e+01,
         2.34000000e+01,   2.36000000e+01,   2.38000000e+01,
         2.40000000e+01,   2.42000000e+01,   2.44000000e+01,
         2.46000000e+01,   2.48000000e+01,   2.50000000e+01,
         2.52000000e+01,   2.54000000e+01,   2.56000000e+01,
         2.58000000e+01,   2.60000000e+01,   2.62000000e+01,
         2.64000000e+01,   2.66000000e+01,   2.68000000e+01,
         2.70000000e+01,   2.72000000e+01,   2.74000000e+01,
         2.76000000e+01,   2.78000000e+01,   2.80000000e+01,
         2.82000000e+01,   2.84000000e+01,   2.86000000e+01,
         2.88000000e+01,   2.90000000e+01,   2.92000000e+01,
         2.94000000e+01,   2.96000000e+01,   2.98000000e+01,
         3.00000000e+01,   3.02000000e+01,   3.04000000e+01,
         3.06000000e+01,   3.08000000e+01,   3.10000000e+01,
         3.12000000e+01,   3.14000000e+01,   3.16000000e+01,

...

         4.23000000e+02,   4.24000000e+02,   4.25000000e+02,
         4.26000000e+02,   4.27000000e+02,   4.28000000e+02,
         4.29000000e+02,   4.30000000e+02,   4.31000000e+02,
         4.32000000e+02,   4.33000000e+02,   4.34000000e+02,
         4.35000000e+02,   4.36000000e+02,   4.37000000e+02,
         4.38000000e+02,   4.39000000e+02,   4.40000000e+02,
         4.41000000e+02,   4.42000000e+02,   4.43000000e+02,
         4.44000000e+02,   4.45000000e+02,   4.46000000e+02,
         4.47000000e+02,   4.48000000e+02,   4.49000000e+02,
         4.50000000e+02,   4.51000000e+02,   4.52000000e+02,
         4.53000000e+02,   4.54000000e+02,   4.55000000e+02,
         4.56000000e+02,   4.57000000e+02,   4.58000000e+02,
         4.59000000e+02,   4.60000000e+02,   4.61000000e+02,
         4.62000000e+02,   4.63000000e+02,   4.64000000e+02,
         4.65000000e+02,   4.66000000e+02,   4.67000000e+02,
         4.68000000e+02,   4.69000000e+02,   4.70000000e+02,
         4.71000000e+02,   4.72000000e+02,   4.73000000e+02,
         4.74000000e+02,   4.75000000e+02,   4.76000000e+02,
         4.77000000e+02,   4.78000000e+02,   4.79000000e+02,
         4.80000000e+02,   4.81000000e+02,   4.82000000e+02,
         4.83000000e+02,   4.84000000e+02,   4.85000000e+02,
         4.86000000e+02,   4.87000000e+02,   4.88000000e+02,
         4.89000000e+02,   4.90000000e+02,   4.91000000e+02,
         4.92000000e+02,   4.93000000e+02,   4.94000000e+02,
         4.95000000e+02,   4.96000000e+02,   4.97000000e+02,
         4.98000000e+02,   4.99000000e+02,   5.00000000e+02,
         5.01000000e+02,   5.02000000e+02,   5.03000000e+02,
         5.04000000e+02,   5.05000000e+02,   5.06000000e+02,
         5.07000000e+02,   5.08000000e+02,   5.09000000e+02,
         5.10000000e+02,   5.11000000e+02,   5.12000000e+02,
         5.13000000e+02,   5.14000000e+02,   5.15000000e+02,
         5.16000000e+02,   5.17000000e+02,   5.18000000e+02,
         5.19000000e+02,   5.20000000e+02,   5.21000000e+02,
         5.22000000e+02,   5.23000000e+02,   5.24000000e+02,
         5.25000000e+02,   5.26000000e+02,   5.27000000e+02,
         5.28000000e+02,   5.29000000e+02,   5.30000000e+02,
         5.31000000e+02,   5.32000000e+02,   5.33000000e+02,
         5.34000000e+02,   5.35000000e+02,   5.36000000e+02,
         5.37000000e+02,   5.38000000e+02,   5.39000000e+02,
         5.40000000e+02,   5.41000000e+02,   5.42000000e+02,
         5.43000000e+02,   5.44000000e+02,   5.45000000e+02,
         5.46000000e+02,   5.47000000e+02,   5.48000000e+02,
         5.49000000e+02,   5.50000000e+02,   5.51000000e+02,
         5.52000000e+02,   5.53000000e+02,   5.54000000e+02,
         5.55000000e+02,   5.56000000e+02,   5.57000000e+02,
         5.58000000e+02,   5.59000000e+02,   5.60000000e+02,
         5.61000000e+02,   5.62000000e+02,   5.63000000e+02,
         5.64000000e+02,   5.65000000e+02,   5.66000000e+02,
         5.67000000e+02,   5.68000000e+02,   5.69000000e+02,
         5.70000000e+02,   5.71000000e+02,   5.72000000e+02,
         5.73000000e+02,   5.74000000e+02,   5.75000000e+02,
         5.76000000e+02,   5.77000000e+02,   5.78000000e+02,
         5.79000000e+02,   5.80000000e+02,   5.81000000e+02,
         5.82000000e+02,   5.83000000e+02,   5.84000000e+02,
         5.85000000e+02,   5.86000000e+02,   5.87000000e+02,
         5.88000000e+02,   5.89000000e+02,   5.90000000e+02,
         5.91000000e+02,   5.92000000e+02,   5.93000000e+02,
         5.94000000e+02,   5.95000000e+02,   5.96000000e+02,
         5.97000000e+02,   5.98000000e+02,   5.99000000e+02,
         6.00000000e+02])
}}}

{{{id=50|
numpy.where(ptime > 6.00000000e+01)
///
(array([319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331,
       332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344,
       345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357,
       358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370,
       371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383,
       384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396,
       397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409,
       410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422,
       423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435,
       436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448,
       449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461,
       462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
       475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487,
       488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500,
       501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513,
       514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526,
       527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539,
       540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552,
       553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565,
       566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578,
       579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591,
       592, 593, 594, 595, 596, 597, 598, 599, 600]),)
}}}

{{{id=51|

///
}}}
