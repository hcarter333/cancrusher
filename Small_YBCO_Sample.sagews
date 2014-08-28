load('./crushersim.sage')
#This simulation calculates result for the coil size appropriate to the small 
#YBCO sample
crushR = Crusher()
crushR.setr(0.8e-2)
crushR.set_movecan(False)
crushR.setTemp(4.2)
crushR.simulate(200)
#Find the maximum current

#Now plot the results
Sc = list_plot(crushR.coilOutTime[0:200, 0:2], axes_labels=['$\mu Seconds$','$kA$'], color='blue', legend_label = '$4.2 K$')
Tc = list_plot(crushR.coilOutTemp[0:197, 0:2], axes_labels=['$\mu Seconds$','$degrees K$'], color='blue', legend_label = '$4.2 K$')
Bc = crushR.get_sphere_field_z()
show(Sc)
show(Tc)
show(Bc)
