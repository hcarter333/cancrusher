load('./crushersim.sage')
#This simulation checks that setr doesn't break anything
#passes as of 2014/7/30 16:31 CST
crushBig = Crusher()
crushBig.setr(4.445e-2)
crushBig.setVnought(5000)
crushBig.setTemp(4.2)
crushBig.set_movecan(False)
crushBig.simulate(298)
Sc = list_plot(crushBig.coilOutTime[0:298, 0:2], axes_labels=['$\mu Seconds$','$kA$'], color='blue', legend_label = '$r=4.445cm$')

crushSmall = Crusher()
crushSmall.setr(1.9055e-2)
crushSmall.setVnought(5000)
crushSmall.setTemp(4.2)
crushSmall.set_movecan(False)
crushSmall.simulate(298)
#The currents should overlay each other
Tc = list_plot(crushSmall.coilOutTime[0:298, 0:2], color='red',  legend_label = '$r=1.9055cm$')
show(Sc + Tc)
np.where(crushSmall.coilOutTime[0:200, 1:2] == crushSmall.coilOutTime[0:200, 1:2].max())
