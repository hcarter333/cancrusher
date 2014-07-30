
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
        self.r = np.array([3.43e-2, 3.43e-2, 3.43e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2, 3.12e-2], dtype=float)
    #Here's the array for coil z coordinates
        self.z = np.array([-2.67e-03, 0.0e-00, 2.67e-3,-6.67e-3, -4.0E-3, -1.33E-3, 1.33E-3, 4.0E-3, 6.67E-3], dtype=float)
        #Here's the array for coil-can separtion
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
            self.res[i] = 1e3*2*pi*self.rcan*rstvty(self.temp[i])/(self.thickness*self.delta)


        #now setup the variables for stepping through the simulation
        self.ntim = 600
        self.ddt = 0.00002
        self.nchange = 20
        
        self.movecan = True

        self.coilI = np.array(range(self.ntim+1), dtype=float)
        self.coilOutTime = np.array(range((self.ntim+1)*2), dtype=float)
        self.coilOutTime.shape = (self.ntim+1, 2)
        self.coilOutTemp= np.array(range((self.ntim+1)*2), dtype=float)
        self.coilOutTemp.shape = (self.ntim+1, 2)
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


    #Function sets up the intitial temperature of all the movable coils
    #For now, it set them all to the same temperature
    def setTemp(self, Temperature):
        for i in range(0,7):
            self.temp[i] = Temperature


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
            rho = rstvty(tt)
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
///
}}}

{{{id=66|
#This simulation compares the crusher currents and temperatures at room temp to 4.2 K
crushHe42 = Crusher()
crushHe42.setTemp(4.2)
crushHe42.simulate(298)
Sc = list_plot(crushHe42.coilOutTime[0:298, 0:2], axes_labels=['$\mu Seconds$','$kA$'], color='blue', legend_label = '$4.2 K$')

crushtest = Crusher()
crushtest.simulate(298)
Tc = list_plot(crushtest.coilOutTime[0:298, 0:2], color='red',  legend_label = '$293 K$')
show(Sc + Tc)

#Now plot the maximum temperatures of the can in each case
St = list_plot(crushHe42.coilOutTemp[0:298, 0:2], axes_labels=['$\mu Seconds$','$degrees K$'], color='blue', legend_label = '$4.2 K$')
Tt = list_plot(crushtest.coilOutTemp[0:298, 0:2], color='red',  legend_label = '$293 K$')
show(St + Tt)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
<html><font color='black'><img src='cell://sage1.png'></font></html>
}}}

{{{id=64|
#Now, let's look at the magnetic field along a constant spherical, not cylindrical 
#radius as the z goes from 0 to just less than the radius of the coil
Bzradius(rcoil, z, curren) = dbcoilbz(rcoil, z, rcoil*(1-(z/rcoil)^2)^0.5, curren)
Hc(z) = 0.760
Hcoil = plot(Bzradius(2.43e-2, z, 42000), 0.001, 2.42e-2, axes_labels=['$meters$','$kGauss$'], legend_label = '$H radius = 2.43 cm$')
Hcoil2 = plot(Bzradius(3.43e-2, z, 42000), 0.001, 3.42e-2, axes_labels=['$meters$','$kGauss$'], legend_label = '$H radius = 3.43 cm$')
Hcritical = plot(Hc(z), 0.001, 3.42e-2, color = 'red',  legend_label = '$Hc Pb$')
show(Hcoil + Hcoil2 + Hcritical)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=63|
#The following will hopefully keep me from having to retype the object name 
#every time it changes
gr = crushHe42

show(list_plot(gr.rstor[0:gr.nmov, 1], plotjoined=True) + list_plot(gr.rstor[0:gr.nmov, 77], plotjoined=True) + list_plot(gr.rstor[0:gr.nmov, 138], plotjoined=True) + list_plot(gr.rstor[0:gr.nmov, 198], plotjoined=True) + list_plot(gr.rstor[0:crushtest.nmov, 258], plotjoined=True) )
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=65|
#Cool, now, create an object
#nm_crushtest = Crusher()
#nm_crushtest.set_movecan(False)
#nm_crushtest.simulate(372)

crushtest = Crusher()
crushtest.simulate(298)
///
}}}

{{{id=57|
S = list_plot(crushHe42.coilOutTime[0:298, 0:2], axes_labels=['$\mu Seconds$','$kA$'], color='blue', legend_label = '$4.2 K$')
T = list_plot(crushtest.coilOutTime[0:298, 0:2], color='red',  legend_label = '$293 K$')
show(S + T)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=61|
nomove = list_plot(nm_crushtest.coilOutTime[0:372, 0:2], color='red')
move = list_plot(crushtest.coilOutTime[0:372, 0:2], color='blue')
show(nomove + move)
///
<html><font color='black'><img src='cell://sage0.png'></font></html>
}}}

{{{id=62|
#Run a simulation up to 193 steps,(about 3.5 uSeconds), and map the z = 0 magnetic field
magtest = Crusher()
magtest.simulate(372)
///
}}}

{{{id=34|
#magnet mapping code for peak current
#create a loop to walk around the circumference of the can and call dbcoilbz
#Load the results into an array with one result and two coordinates
///
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

///
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
