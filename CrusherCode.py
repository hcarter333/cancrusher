
{{{id=56|
import numpy as np
from numpy.linalg import inv

perm0 = 1.25663706*1e-6


#Simple Sage objects that live outside the class for now
rstvty(temperature)=(4.5+1.133*10^-2*temperature)*10^-8
#Returns the specific heat in J/Kg/K
specheat(temperature) = 4.18*(0.16+1.5e-4*temperature)*1000.0    

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


global mfull


class Crusher:

    def __init__(self):

        #Setup an array for coil dimensions
        self.r = numpy.array([3.43e-2, 3.43e-2, 3.43e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2], dtype=float)
    #Here's the array for coil z coordinates
        self.z = numpy.array([-2.67e-03, 0.0e-00, 2.67e-3,-6.67e-3, -4.0E-3, -1.33E-3, 1.33E-3, 4.0E-3, 6.67E-3], dtype=float)
        #Here's the array for coil-can separtion
        self.dr = numpy.array([1.33E-3, 1.33E-3, 1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3], dtype=float)
    





        #setup the constants
        self.nmov = 6
        self.nfix = 3

        #need to change this for Pb
        self.density = 2824.0

        self.thickness = .00016


        self.delta = abs(self.z[self.nmov+1] - self.z[self.nmov+2])

        #Here are two more arrays that are used, but I'm not sure what for yet
        #vr is the radial velocity of the can?  For some reason it included enough entries 
        #for the can and the fixed coil
        #vz is the z velocity, but it looks unused.  It looks like the axial force 
        #implementation was left incomplete in the initial version
        self.vr = numpy.array([0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)
        self.vz = numpy.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], dtype=float)

        #NOTE:  r, z, and vr will all be modified by move can

        #array stores the temerature of the can under each coil
        self.temp = numpy.array([293.0, 293.0, 293.0, 293.0, 293.0, 293.0, 293.0], dtype=float)

        #in the following, we'll need the resistances and capacitances that were pulled 
        #from the input data file
        #The following are the resistance, capacitance, inductance, admittance, 
        #current, and charge arrays with one entry per can modelling coil
        #They have been filled in here with the values from the model data file
        #In the model data file, there are 7 values as opposed to the six I expected
        #I need todo check on this
        self.res = numpy.array([8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)

        #Just changed the capacitance to 500 from 50
        self.cap = numpy.array([0.0004, -1, -1, -1, -1, -1, -1], dtype=float)

        self.lckt = numpy.array([0.05, 0, 0, 0, 0, 0, 0], dtype=float)

        self.cur = numpy.array([0.0, 0, 0, 0, 0, 0, 0], dtype=float)
        self.ccur = numpy.array(range(self.nmov+self.nfix), dtype=float)

        #Changed to 3000 from 5000
        self.V0 = numpy.array([3000.0, 0, 0, 0, 0, 0, 0], dtype=float)

        #setup the admittance matrix.  I've made this a single statement here, keep in mind that 
        #there's a there is a loop treatment in the original code
        self.adc = numpy.array([1/self.cap[0], 0, 0, 0, 0, 0, 0], dtype=float)
        self.qq = numpy.array([self.V0[0]*self.cap[0], 0, 0, 0, 0, 0, 0], dtype=float)

        #Here's a clue to how the radius positions of coil and can should be worked out
        #rcan = r[self.nfix+1]
        #the indexing above seems odd
        self.rcan = self.r[self.nfix+1]
        self.mass = 2*pi*self.rcan*self.density*self.thickness*self.delta

        #the next step is to fill in the res array with calls to rstvty
        #This begins to make more sense after a fahion.  The one bigger is intentional, the 
        #first value in the data file was real, all the others are filled in
        for i in range(1, self.nmov+1):
            self.res[i] = 1e3*2*pi*self.rcan*rstvty(self.temp[i])/(self.thickness*self.delta)


        #now setup the variables for stepping through the simulation
        self.ntim = 600
        self.ddt = 0.00002
        self.nchange = 20
        
        self.movecan = True

        self.coilI = numpy.array(range(self.ntim+1), dtype=float)
        self.coilOutTime = numpy.array(range((self.ntim+1)*2), dtype=float)
        self.coilOutTime.shape = (self.ntim+1, 2)
        self.bzero = numpy.array(range(self.ntim+1), dtype=float)
        self.zeroline = numpy.array(range(self.ntim+1), dtype=float)
        self.heatenrg = numpy.array(range(self.ntim+1), dtype=float)
        self.work = numpy.array(range(self.ntim+1), dtype=float)
        self.enrgtot = numpy.array(range(self.ntim+1), dtype=float)

        #store the current time in microseconds
        self.ptime = numpy.array(range(self.ntim+1), dtype=float)
        self.denrg = 0.0
        self.dheat = 0.0
        self.dwork = 0.0

        self.dt = self.ddt
        self.time = 0.0
        self.aaa = 1.0
        self.ntim1 = self.ntim-1

        #Here are the arrays that will hold the data from the simulation
        self.rstor = numpy.array(range(self.nmov*(self.ntim+1)), dtype=float)
        self.rstor.shape = (self.nmov, self.ntim+1)

        self.zstor = numpy.array(range(self.nmov*(self.ntim+1)), dtype=float)
        self.zstor.shape = (self.nmov, self.ntim+1)

        self.cntr = 0

        self.ptime[0] = 0.0
        self.coilI[0] = 0.0
        self.bzero[0] = 0.0
        self.zeroline[0] = 0.0
        self.pctr = 0
        self.heatenrg[0] = 0.0
        self.enrgtot[0] = 0.0
        self.work[0] = 0.0

    
    

        #load the rstor and zstor arrays with the initial positions of the moving coils
        for jj in range(0, self.nfix):
            self.rstor[jj,0] = self.r[jj+self.nfix]
            self.zstor[jj,0] = self.z[jj+self.nfix]

        #arrays to hold the reduced mutual inductance matrix
        self.mm = numpy.array(range((self.nmov+1)^2), dtype=float)
        self.mm.shape=(self.nmov+1,self.nmov+1)
        self.mmold = self.mm.copy()
        
        #setting up the mutual inductance raw flux array
        self.mfull = numpy.array(range((self.nmov+self.nfix)*(self.nmov+self.nfix)), dtype=float)
        self.mfull.shape=(self.nmov+self.nfix,self.nmov+self.nfix)
        
        #initialize the raw mutual inductance array
        self.mfull = self.find_mutual_inductance(self.mfull)

        self.mm = self.make_reduced_matrix(self.mm, self.mfull)


    #The mutual inductance array has to be calculated multiple times, so I'm putting 
    #it into a Python funcion, (I hope)
    def find_mutual_inductance(self, mfullarray):

        #first for the off-diagonal elements
        for i in range(0, self.nmov+self.nfix-1):
           for j in range(i+1, self.nmov+self.nfix):
                mfullarray[i,j] = dbcoilflux(self.r[i], self.z[j]-self.z[i], self.r[j]-self.dr[j], 1)
    
        #next for the diagonal elements
        for i in range(0, self.nmov+self.nfix):
            mfullarray[i,i] = dbcoilflux(self.r[i], self.z[i]-self.z[i], self.r[i]-self.dr[i], 1)
   
        #finally copy the off diagonal elements to the other side of the array
        #crap, there appears to be a bug below.  Moving mfull -> mfullarray
        for i in range(1, self.nmov+self.nfix):
            for j in range(0, i):
                mfullarray[i,j] = mfullarray[j,i]
            
        return mfullarray

    #also create a function here that encapsulates the reduced array creation
    #since this will need to be called multiple times
    def make_reduced_matrix(self, mreduced, mfullwork):
        self.mmold = self.mm.copy()
        jmwork = numpy.array(range(self.nmov), dtype=float)
        #nt contains the fixed coil portion of the mutual inductance array
        #it's actually just a 
        ntwork = self.mfull[0:(self.nfix), 0:(self.nfix)].copy()
        nttwork = sum(ntwork)
    
        #jt is used to total the mutual inductance associated with each moveable coil
        #the result is placed in the jmwork array
        for i in range(0, self.nmov):
            jtwork=self.mfull[i+self.nfix,0:self.nfix].copy()
            jmwork[i]=sum(jtwork)
    
        #Here's why mm was dimensioned one bigger than it should have been
        #I don't have a good reason for this yet other than it came out of the old code
        #notice, however that we do overwrite the 0 row and column next
        mreduced = mfullwork[self.nfix-1:self.nmov+self.nfix, self.nfix-1:self.nmov+self.nfix].copy()
        mreduced = 1e6*mreduced
        #print "mfullwork"
        #print mfullwork
        #print "mreduced"
        #print mreduced
        #mm = mm*1e6

        for i in range(1, self.nmov+1):
            mreduced[i,0] = 1e6*jmwork[i-1]
            mreduced[0,i] = 1e6*jmwork[i-1]
    
        mreduced[0,0] = 1e6*nttwork

        return mreduced


    #Calculate_current function moved up to the OO world
    def compute_current(self):

        #first, we'll need some working arrays
        for i in range(0,self.nfix):
            self.ccur[i] = self.cur[0]
    
        for i in range(self.nfix,self.nfix+self.nmov):
            self.ccur[i] = self.cur[i-self.nfix+1]
    
        #variable names are copied from the original code
        nc1 = self.nmov+1
        nc11 = self.nmov
        idot = numpy.array(range(nc1), dtype=float)
 
        vv = numpy.array(range(nc1), dtype=float)
    
        #and here's the rather large array
        #that in fact never gets used for anything in the code, so I'm going to 
        #comment it out.  It does look like it would make a good debug tool though
        #pcur = numpy.array(range(nc1*ntim), dtype=float)
        #pcur.shape = (nc1, ntim)
    
        yp = numpy.array(range(self.ntim), dtype=float)
        tim = numpy.array(range(self.ntim), dtype=float)
        tim[0] = 0
    
    
        #add the external inductance into the reduced mutual inductance array
        for i in range(0, self.nmov+1):
            self.mm[i,i] = self.mm[i,i]+self.lckt[i]
    
    
        #NOTE:  Here's why we left a copy of the original mm hanging around in mmold
        dmmdt = (self.mm-self.mmold)/self.dt
    
        #make a copy of mm
        mident = self.mm.copy
    
        #now, invert mm
        minv = inv(self.mm)
    
        #make the inverse symmetric
        for i in range(0, self.nmov):
            for j in range(i+1, self.nmov+1):
                minv[i,j] = (minv[i,j]+minv[j,i])/2.0
                minv[j,i] = minv[i,j]
            
        #now compute the currents in the coils
        #the following is an operation on matrices, but appears not to be a matrix multiply
        #which would be denoted in IDL by #  I'm guessing it just multiplies element by element
        #eyer looks like a voltage
        eyer = self.cur*self.res

        #c = q/v so quec looks like a voltage
        quec = self.qq*self.adc

        #The following actually is a matrix multiply
        #it's cute because it's not l*di/dt
        idldt = np.dot(dmmdt, self.cur)
    
        #Make the arrays symmetric as usual
        for i in range(1, (self.nmov/2)+1):
            iii=self.nmov+1-i
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
        for i in range(1, (self.nmov/2)+1):
            iii=self.nmov+1-i
            idot[i]=(idot[i]+idot[iii])/2 
            idot[iii]=idot[i]
        
        #This is an array operation
        self.cur=self.cur+idot*self.dt
    
        self.qq=self.qq-self.adc*self.cap*self.cur*self.dt
        self.coilI[self.cntr]=self.cur[0] 
        self.coilOutTime[self.cntr,0]= self.ptime[self.cntr]
        self.coilOutTime[self.cntr,1] = self.cur[0]
    
        denrg=(self.cur[0]*self.qq[0]*self.dt)/self.cap[0]
    
        self.zeroline[self.cntr]=0.0
    
        #this variable is actually local!
        sbz = 0.0
    
        for j in range(0, self.nfix+self.nmov):
            zz = self.z[j]
            rr = 0.0
            rso = self.r[j]
            curr = 1e3*self.ccur[j]
            #implement the portion of dbcoil that returns bz
            sbz = sbz + dbcoilbzrzero(rso, zz, rr, curr)
            #sbz = sbz + bz
        
        #This is the field in kGauss
        self.bzero[self.cntr] = 10*sbz
        tms = 1e3*self.time
    
        #Calculate the heat
        dheat = 0.0
        for i in range(1, self.nmov+1):
            tt=self.temp[i]
            sig = specheat(tt)
            enrg = self.cur[i]^2*self.res[i]*self.dt
            dheat = dheat + enrg
            dtemp = dheat/(sig*self.mass)
            self.temp[i] = self.temp[i] + dtemp
            tt = self.temp[i]
            rho = rstvty(tt)
            #Cool!  Updating the resistance based on the temperature change
            #in the can
            self.res[i] = 1e3*rho*2*pi*self.rcan/(self.thickness*self.delta)
            
    def move_can(self):

    
        brg = numpy.array(range(self.nfix + self.nmov), dtype=float)
        bzt = numpy.array(range(self.nfix + self.nmov), dtype=float)
    
        for i in range(1,self.nmov+1):
            sbz = 0.0
            imov = self.nfix + i - 1
            for j in range(0,self.nfix+self.nmov):
                if j==imov:
                     break
                zz = self.z[imov]-self.z[j]
                rr = self.r[imov]
                rso = self.r[j]
                curr = 1e3*self.ccur[j]
                #if rr is zero, the can has broken in half and all bets are off anyway           
                #if rr == 0:
                #    sbz = sbz + dbcoilbzrzero(rso, zz, rr, curr)
                #else:
                sbz = sbz + dbcoilbz(rso, zz, rr, curr)
                
            bzt[i+self.nfix-1] = sbz
        dwork = 0.0

        for i in range(1,(self.nmov/2)+1):
            #This looks like we're hitting opposite edges of the moving coils, the can, and 
            #working towards the center
            #print "i"
            #print i
            ii=self.nfix+i-1
            iii = (self.nfix + self.nmov) - i
        
            #Find force in kNewtons
            forcer = bzt[ii]*self.ccur[ii]*2*pi*self.r[ii]
            
            #Get the differential in velocity from the force and mass
            dvr = forcer*self.dt/self.mass
            vrnew = self.vr[ii]+dvr
        
            #get the new r position using the velocity
            rnew=self.r[ii]+self.vr[ii]*1e-3*self.dt

        
            #Find work in Joules
            dwork=dwork+2*forcer*vrnew*self.dt
            self.vr[ii]=vrnew
            self.r[ii] = rnew
            self.vr[iii]=vrnew
            self.r[iii]=rnew
        #print "completed loop"        
    
    def simulate(self, endtime):
        for kk in range(0,endtime):
            #if the counter has advanced beyond nchange, then make the time step larger
            if self.cntr >= self.nchange:
                self.dt = self.ddt*10
            print self.cntr
            self.cntr = self.cntr + 1
            self.time = self.time + self.dt
            #store the current time in microseconds
            self.ptime[self.cntr] = self.time*1e3
    
            #Even though these funcitons have been called in initialize, it's important to call them even on 
            #the first loop through here.  Otherwise, mmold winds up with junk in it.
            #now, find the mutual inductance

            self.mfull = self.find_mutual_inductance(self.mfull)
            #then, reduce the mutual inductance array again

            self.mm = self.make_reduced_matrix(self.mm, self.mfull)
            #now, finally, the first new simulation step, compute the currents
            print "Time in microseconds"
            print self.ptime[self.cntr]
            self.compute_current()
            if self.movecan:
                self.move_can()
    
            #track the heat and work for this step
            self.heatenrg[self.cntr]=self.heatenrg[self.cntr+1]+self.dheat
            self.work[self.cntr]=self.work[self.cntr-1]+self.dwork
            self.enrgtot[self.cntr]=self.enrgtot[self.cntr-1]+self.denrg
            for jj in range(0,self.nmov):
                jjmov = jj + self.nfix
                self.rstor[jj,kk+1] = self.r[jjmov]
                self.zstor[jj,kk+1] = self.z[jjmov]

    def set_movecan(self, move):
        self.movecan = move
#This concludes all the setup of variables that is done by the initialize function



#now take mfull created above and use it to generate the reduced mutual induction matrix mm
#mm is dimensioned for only the moving coils plus one extra row and column
#ultimately we'll need to track the immediately previous version of mm
#in the original code this is called mmold, and here too!




#The above was to provide the magnetic flux at one coil due to another for finding 
#the mutual inductance of a set of two coils.  The following models the coils and 
#calls the above functions to fill in the mutual induction matrix
#Utilize a python for loop here and I guess go ahead and use sage matrices, although 
#I'm not sure there's not a reason to use numpy matrices or plain old Python arrays
#The documentation for numpy arrays is the easiest for me to read, so I'm going that 
#direction.  I went ahead and put the import at the top of the worksheet

    
#calls to initialize the code.  This should eventually be put into an initialize function
#First, let's see if it reads in ok
#initialize()
///
}}}

{{{id=57|
#Cool, now, create an object
nm_crushtest = Crusher()
nm_crushtest.set_movecan(False)
nm_crushtest.simulate(372)

crushtest = Crusher()
crushtest.simulate(372)
///
WARNING: Output truncated!  
<html><a target='_new' href='/home/admin/4/cells/57/full_output.txt' class='file_link'>full_output.txt</a></html>



0
Time in microseconds
0.02
1
Time in microseconds
0.04
2
Time in microseconds
0.06
3
Time in microseconds
0.08
4
Time in microseconds
0.1
5
Time in microseconds
0.12
6
Time in microseconds
0.14
7
Time in microseconds
0.16
8
Time in microseconds
0.18
9
Time in microseconds
0.2
10
Time in microseconds
0.22
11
Time in microseconds
0.24
12
Time in microseconds
0.26
13
Time in microseconds
0.28
14
Time in microseconds
0.3
15
Time in microseconds
0.32
16
Time in microseconds
0.34
17
Time in microseconds
0.36
18
Time in microseconds
0.38
19
Time in microseconds

...

352
Time in microseconds
67.0
353
Time in microseconds
67.2
354
Time in microseconds
67.4
355
Time in microseconds
67.6
356
Time in microseconds
67.8
357
Time in microseconds
68.0
358
Time in microseconds
68.2
359
Time in microseconds
68.4
360
Time in microseconds
68.6
361
Time in microseconds
68.8
362
Time in microseconds
69.0
363
Time in microseconds
69.2
364
Time in microseconds
69.4
365
Time in microseconds
69.6
366
Time in microseconds
69.8
367
Time in microseconds
70.0
368
Time in microseconds
70.2
369
Time in microseconds
70.4
370
Time in microseconds
70.6
371
Time in microseconds
70.8
}}}

{{{id=34|
nomove = list_plot(nm_crushtest.coilOutTime[0:372, 0:2], color='red')
move = list_plot(crushtest.coilOutTime[0:372, 0:2], color='blue')
show(nomove + move)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=60|
show(list_plot(crushtest.rstor[0:crushtest.nmov, 1], plotjoined=True) + list_plot(crushtest.rstor[0:crushtest.nmov, 77], plotjoined=True) + list_plot(crushtest.rstor[0:crushtest.nmov, 138], plotjoined=True) + list_plot(crushtest.rstor[0:crushtest.nmov, 198], plotjoined=True) + list_plot(crushtest.rstor[0:crushtest.nmov, 258], plotjoined=True) + list_plot(crushtest.rstor[0:crushtest.nmov, 318], plotjoined=True) + list_plot(crushtest.rstor[0:crushtest.nmov, 372], plotjoined=True))
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=59|

///
}}}

{{{id=58|

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
        
    #This is an array operation
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
    #print "nmov"
    #print nmov
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
        #print "rii"
        #print r[ii]
        #print "vrii"
        #print vr[ii]
        #print"dt"
        #print dt
        
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
#Initialize the code
initialize()

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
Traceback (most recent call last):    
  File "", line 1, in <module>
    
  File "/tmp/tmpaWWE0L/___code___.py", line 4, in <module>
    initialize()
  File "/tmp/tmpeakUYR/___code___.py", line 138, in initialize
    mfull = find_mutual_inductance(mfull)
  File "/tmp/tmpeakUYR/___code___.py", line 181, in find_mutual_inductance
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
list_plot(coilOutTime[0:372, 0:2], color='red')
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=36|
save_no_move = coilOutTime[0:372, 0:2].copy()
///
}}}

{{{id=40|
show(list_plot(save_no_move[0:371, 0:2]) + list_plot(coilOutTime[0:371, 0:2], color = 'red'))
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=41|
rstor[0:nmov, 1]
///
array([ 0.03120001,  0.03120001,  0.03120001,  0.03120001,  0.03120001,
        0.03120001])
}}}

{{{id=42|
show(list_plot(rstor[0:nmov, 1], plotjoined=True) + list_plot(rstor[0:nmov, 77], plotjoined=True) + list_plot(rstor[0:nmov, 138], plotjoined=True) + list_plot(rstor[0:nmov, 198], plotjoined=True) + list_plot(rstor[0:nmov, 258], plotjoined=True) + list_plot(rstor[0:nmov, 318], plotjoined=True) + list_plot(rstor[0:nmov, 372], plotjoined=True))
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
S = sphere(size=.5, color='yellow')
///
}}}

{{{id=52|
show(S)
///
}}}

{{{id=53|
T = Cylindrical('height', ['radius', 'azimuth'])
r, theta, z = var('r theta z')
T.transform(radius=r, azimuth=theta, height=z)
(r*cos(theta), r*sin(theta), z)
plot3d(9-r^2, (r, 0, 3), (theta, 0, pi), transformation=T, cmap='hsv')
///
}}}

{{{id=54|
from sage.plot.colors import get_cmap
get_cmap('hsv')
///
<matplotlib.colors.LinearSegmentedColormap object at 0xc37eccc>
}}}

{{{id=55|

///
}}}
