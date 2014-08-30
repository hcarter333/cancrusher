import numpy as np
from numpy.linalg import inv

perm0 = 1.25663706*1e-6


#Simple Sage objects that live outside the class for now
#Moved to object
#rstvty(temperature)=(4.5+1.133*10^-2*temperature)*10^-8

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

#r direction magnetic field
#for the moment, never pass this method an r value of zero
brfac(rcoil, zc, rc) = ((rcoil^2+rc^2)+zc^2)/(((rcoil-rc)^2)+zc^2)
dbcoilbr(rcoil, zc, rc, curren) = (2*perm0*curren/(2*pi))*zc*(brfac(rcoil, zc, rc)*elliptic_ec(argm(rcoil, zc, rc)) - elliptic_kc(argm(rcoil, zc, rc)))/(t(rcoil, zc, rc)^0.5)/rc

#method for finding the magnetic field vs the z coponent
Bzradius(rcoil, z, curren) = dbcoilbz(rcoil, z, rcoil*(1-(z/rcoil)^2)^0.5, curren)

#Same thing for the radial component
Brradius(rcoil, z, curren) = dbcoilbr(rcoil, z, rcoil*(1-(z/rcoil)^2)^0.5, curren)

global mfull


class Crusher:

    def rstvty(self, Temperature):
        if self.superconduct == False:
            return (4.5+1.133*10^-2*Temperature)*10^-8
        
        return 0
        
    def __init__(self):

        #Make driven coil superconducting
        self.superconduct = False
        
        #Setup an array for coil dimensions
        self.r = np.array([3.43e-2, 3.43e-2, 3.43e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2], dtype=float)
        #Here's the array for coil z coordinates
        self.z = np.array([-2.67e-03, 0.0e-00, 2.67e-3,-6.67e-3, -4.0E-3, -1.33E-3, 1.33E-3, 4.0E-3, 6.67E-3], dtype=float)
        #Here's the array for coil-can separtion?
        self.dr = np.array([1.33E-3, 1.33E-3, 1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3, -1.33E-3], dtype=float)
    
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
        self.vr = np.array([0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)
        self.vz = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], dtype=float)

        #NOTE:  r, z, and vr will all be modified by move can

        #array stores the temerature of the can under each coil
        self.temp = np.array([293.0, 293.0, 293.0, 293.0, 293.0, 293.0, 293.0], dtype=float)

        #in the following, we'll need the resistances and capacitances that were pulled 
        #from the input data file
        #The following are the resistance, capacitance, inductance, admittance, 
        #current, and charge arrays with one entry per can modelling coil
        #They have been filled in here with the values from the model data file
        #In the model data file, there are 7 values as opposed to the six I expected
        #I need todo check on this
        self.res = np.array([8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)

        #Just changed the capacitance to 500 from 50
        self.cap = np.array([0.0004, -1, -1, -1, -1, -1, -1], dtype=float)

        self.lckt = np.array([0.05, 0, 0, 0, 0, 0, 0], dtype=float)

        self.cur = np.array([0.0, 0, 0, 0, 0, 0, 0], dtype=float)
        self.ccur = np.array(range(self.nmov+self.nfix), dtype=float)

        #Changed to 3000 from 5000
        self.V0 = np.array([3000.0, 0, 0, 0, 0, 0, 0], dtype=float)

        #setup the admittance matrix.  I've made this a single statement here, keep in mind that 
        #there's a there is a loop treatment in the original code
        self.adc = np.array([1/self.cap[0], 0, 0, 0, 0, 0, 0], dtype=float)
        self.qq = np.array([self.V0[0]*self.cap[0], 0, 0, 0, 0, 0, 0], dtype=float)

        #Here's a clue to how the radius positions of coil and can should be worked out
        #rcan = r[self.nfix+1]
        #the indexing above seems odd
        self.rcan = self.r[self.nfix+1]
        self.mass = 2*pi*self.rcan*self.density*self.thickness*self.delta

        #the next step is to fill in the res array with calls to rstvty
        #This begins to make more sense after a fahion.  The one bigger is intentional, the 
        #first value in the data file was real, all the others are filled in
        for i in range(1, self.nmov+1):
            self.res[i] = 1e3*2*pi*self.rcan*self.rstvty(self.temp[i])/(self.thickness*self.delta)


        #now setup the variables for stepping through the simulation
        self.ntim = 600
        self.ddt = 0.00002
        self.nchange = 20
        
        self.movecan = True

        self.coilI = np.array(range(self.ntim+1), dtype=float)
        #self.coilOutTime = np.array(range((self.ntim+1)*2), dtype=float)
        #self.coilOutTime.shape = (self.ntim+1, 2)
        self.coilOutTime = np.zeros((self.ntim+1, 2), dtype=float)
        
        #self.coilOutTemp= np.array(range((self.ntim+1)*2), dtype=float)        
        #self.coilOutTemp.shape = (self.ntim+1, 2)
        self.coilOutTemp = np.zeros((self.ntim+1, 2), dtype=float)
        
        self.bzero = np.array(range(self.ntim+1), dtype=float)
        self.zeroline = np.array(range(self.ntim+1), dtype=float)
        self.heatenrg = np.array(range(self.ntim+1), dtype=float)
        self.work = np.array(range(self.ntim+1), dtype=float)
        self.enrgtot = np.array(range(self.ntim+1), dtype=float)

        #store the current time in microseconds
        self.ptime = np.array(range(self.ntim+1), dtype=float)
        self.denrg = 0.0
        self.dheat = 0.0
        self.dwork = 0.0

        self.dt = self.ddt
        self.time = 0.0
        self.aaa = 1.0
        self.ntim1 = self.ntim-1

        #Here are the arrays that will hold the data from the simulation
        self.rstor = np.array(range(self.nmov*(self.ntim+1)), dtype=float)
        self.rstor.shape = (self.nmov, self.ntim+1)

        self.zstor = np.array(range(self.nmov*(self.ntim+1)), dtype=float)
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
        self.mm = np.array(range((self.nmov+1)^2), dtype=float)
        self.mm.shape=(self.nmov+1,self.nmov+1)
        self.mmold = self.mm.copy()
        
        #setting up the mutual inductance raw flux array
        self.mfull = np.array(range((self.nmov+self.nfix)*(self.nmov+self.nfix)), dtype=float)
        self.mfull.shape=(self.nmov+self.nfix,self.nmov+self.nfix)
        
        #initialize the raw mutual inductance array
        self.mfull = self.find_mutual_inductance(self.mfull)

        self.mm = self.make_reduced_matrix(self.mm, self.mfull)

    #Function to change the radius of the sample
    def setr(self, bigr):
        self.r[0] = bigr
        self.r[1] = bigr
        self.r[2] = bigr
        self.r[3] = bigr-.31e-2
        self.r[4] = bigr-.31e-2
        self.r[5] = bigr-.31e-2
        self.r[6] = bigr-.31e-2
        self.r[7] = bigr-.31e-2
        self.r[8] = bigr-.31e-2
        
    
    def setSuperconduct(self, Superconduct):
        self.superconduct = Superconduct
    
    #Function sets up the intitial temperature of all the movable coils
    #For now, it set them all to the same temperature
    def setTemp(self, Temperature):
        for i in range(0,7):
            self.temp[i] = Temperature


    #Set the initial voltage on the capacitor
    def setVnought(self, Voltage):
        self.V0[0] = Voltage
        #Now that the voltage is reset, make sure it effects the capacitor bank
        self.qq =[self.V0[0]*self.cap[0], 0, 0, 0, 0, 0, 0]
    
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

    #function to return the maximum current in the driving coil
    def get_max_current(self):
        return np.amax(self.coilOutTime[0:298, 1])
    
    #Return a plot of the magnetic field along the spherical sample as a function of z
    #This returns the z component only
    def get_sphere_field_z(self):
        return plot(Bzradius(self.r[0], z, self.get_max_current()*1000), 0.001, self.r[0] - (.01 * self.r[0]), axes_labels=['$meters$','$kGauss$'], legend_label = '$H radius = {0} cm$'.format(self.r[0]*100), title = '$Temperature=4.2K$', gridlines=True)

    #Return a plot of the magnetic field along the spherical sample as a function of z
    #This returns the z component only
    def get_sphere_field_r(self):
        return plot(Brradius(self.r[0], z, self.get_max_current()*1000), 0.001, self.r[0] - (.01 * self.r[0]), axes_labels=['$r meters$','$kGauss$'], legend_label = '$H radius = {0} cm$'.format(self.r[0]*100), title = '$Temperature=4.2K$', gridlines=True)
    
    #also create a function here that encapsulates the reduced array creation
    #since this will need to be called multiple times
    def make_reduced_matrix(self, mreduced, mfullwork):
        self.mmold = self.mm.copy()
        jmwork = np.array(range(self.nmov), dtype=float)
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
        idot = np.array(range(nc1), dtype=float)
 
        vv = np.array(range(nc1), dtype=float)
    
        #and here's the rather large array
        #that in fact never gets used for anything in the code, so I'm going to 
        #comment it out.  It does look like it would make a good debug tool though
        #pcur = np.array(range(nc1*ntim), dtype=float)
        #pcur.shape = (nc1, ntim)
    
        yp = np.array(range(self.ntim), dtype=float)
        tim = np.array(range(self.ntim), dtype=float)
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
        maxtemp = 0.0
        for i in range(1, self.nmov+1):
            tt=self.temp[i]
            sig = specheat(tt)
            enrg = self.cur[i]^2*self.res[i]*self.dt
            dheat = dheat + enrg
            dtemp = dheat/(sig*self.mass)
            self.temp[i] = self.temp[i] + dtemp
            
            #store the temp at each time
            #coilOutTemp
            #the implementatin is rather flawed, the ...
            #just fix it I suppose
            #Now, the maximum temperature at each time step is tracked
            if self.temp[i] > maxtemp:
                maxtemp = self.temp[i]
                self.coilOutTemp[self.cntr,0]= self.ptime[self.cntr]
                self.coilOutTemp[self.cntr,1] = maxtemp
            
            tt = self.temp[i]
            rho = self.rstvty(tt)
            #Cool!  Updating the resistance based on the temperature change
            #in the can
            self.res[i] = 1e3*rho*2*pi*self.rcan/(self.thickness*self.delta)
            
    def move_can(self):

    
        brg = np.array(range(self.nfix + self.nmov), dtype=float)
        bzt = np.array(range(self.nfix + self.nmov), dtype=float)
    
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
            #print self.cntr
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
            #print "Time in microseconds"
            #print self.ptime[self.cntr]
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
