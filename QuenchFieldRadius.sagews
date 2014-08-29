load('./crushersim.sage')
#Set simulation time once
tickcount = 290
#This simulation calculates result for the coil size appropriate to the small 
#YBCO sample
crushR = Crusher()
crushR.setr(3.8e-2)
crushR.setVnought(5000)
crushR.set_movecan(False)
crushR.setTemp(4.2)
crushR.simulate(tickcount)
#Find the maximum current

#Now plot the results
Sc = list_plot(crushR.coilOutTime[0:tickcount, 0:2], axes_labels=['$\mu Seconds$','$kA$'], color='blue', legend_label = '$4.2 K$')
Tc = list_plot(crushR.coilOutTemp[0:tickcount - 3, 0:2], axes_labels=['$\mu Seconds$','$degrees K$'], color='blue', legend_label = '$4.2 K$')
Bc = crushR.get_sphere_field_z()
HcPb = plot(.670, 0, .038, color = 'red')
show(Sc)
show(Tc)
show(Bc + HcPb)
