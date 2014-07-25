
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

{{{id=39|
#and inside the coil, being careful to avoid the point r = 0
plot(dbcoilflux(2, 0, rc, 1), 0.01, 2)
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
    for i in range(1,(nmov/2)+1):
        #This looks like we're hitting opposite edges of the moving coils, the can, and 
        #working towards the center
        ii=nfix+i-1
        iii = nfix + nmov - i
        #Find force in kNewtons
        forcer = bzt[ii]*ccur[ii]*2*pi*r[ii]
        #Get the differential in velocity from the force and mass
        dvr = forcer*dt/mass
        vrnew = vr[ii]+dvr
        #get the new r position using the velocity
        rnew=r[ii]+vr[ii]*1e-3*dt
        #Find work in Joules
        dwork=dwork+2*forcer*vrnew*dt
        vr[ii]=vrnew
        r[ii] = rnew
        vr[iii]=vrnew
        r[iii]=rnew
        return vr[5]
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
    compute_current()
    #move_can()
    
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
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58

...

548
549
550
551
552
553
554
555
556
557
558
559
560
561
562
563
564
565
566
567
568
569
570
571
572
573
574
575
576
577
578
579
580
581
582
583
584
585
586
587
588
589
590
591
592
593
594
595
596
597
598
599
Traceback (most recent call last):            dt = ddt*10
  File "", line 1, in <module>
    
  File "/tmp/tmpgu_TF5/___code___.py", line 5, in <module>
    exec compile(u"for kk in range(_sage_const_0 ,_sage_const_600 ):\n    #if the counter has advanced beyond nchange, then make the time step larger\n    if cntr >= nchange:\n        dt = ddt*_sage_const_10 \n    print cntr\n    cntr = cntr + _sage_const_1 \n    time = time + dt\n    #store the current time in microseconds\n    ptime[cntr] = time*_sage_const_1e3 \n    \n    #Even those these funcitons have been called in initialize, it's important to call them even on \n    #the first loop through here.  Otherwise, mmold winds up with junk in it.\n    #now, find the mutual inductance\n    global mfull\n    mfull = find_mutual_inductance(mfull)\n    #then, reduce the mutual inductance array again\n    global mm\n    mm = make_reduced_matrix(mm, mfull)\n    #now, finally, the first new simulation step, compute the currents\n    compute_current()\n    #move_can()\n    \n    #track the heat and work for this step\n    heatenrg[cntr]=heatenrg[cntr+_sage_const_1 ]+dheat\n    work[cntr]=work[cntr-_sage_const_1 ]+dwork\n    enrgtot[cntr]=enrgtot[cntr-_sage_const_1 ]+denrg\n    for jj in range(_sage_const_0 ,nmov):\n        jjmov = jj + nfix\n        rstor[jj,kk+_sage_const_1 ] = r[jjmov]\n        zstor[jj,kk+_sage_const_1 ] = z[jjmov]" + '\n', '', 'single')
  File "", line 24, in <module>
    
IndexError: index 601 is out of bounds for axis 0 with size 601
}}}

{{{id=38|
r
///
array([  3.43000000e-02,   3.43000000e-02,   3.43000000e-02,
        -2.75897764e+18,   3.12000000e-02,   3.12000000e-02,
         3.12000000e-02,   3.12000000e-02,  -2.75897764e+18])
}}}

{{{id=35|
list_plot(coilOutTime[0:599, 0:2])
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=36|
r
///
array([ 0.0343,  0.0343,  0.0343,  0.0312,  0.0312,  0.0312,  0.0312,
        0.0312,  0.0312])
}}}

{{{id=40|

///
}}}
