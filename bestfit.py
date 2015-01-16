
import numpy as np
import pylab as py
	

def meanstd(ell, dl):
	mean_l = sum(ell)/ell.size
	mean_dl= sum(dl)/dl.size
	dl2 = dl*dl
	var_dl = sum(dl2)/dl.size-mean_dl**2
	std_dl = np.sqrt(var_dl)
	stat   = np.array([mean_l, mean_dl, std_dl])
	return stat

def fitfunc(x, alpha, beta, gamma):
	trs = gamma * 145.                 #leq = 145
	if(x > trs):
		return (alpha+beta*x)
	else: 
		return (alpha+beta*trs)/trs*x
	

def chisq(xdata, ydata, sigma, alpha, beta, trs):
	ymodel = [fitfunc(x, alpha, beta, trs) for x in xdata]
	err    = [(ymodel[i]-ydata[i])**2/sigma[i]**2 for i in range(ydata.size)]
	return sum(err)

def bestfit(xdata, ydata, sigma):
	best_par = np.array([8., 1.e-3, 6.])
	best_chi = chisq(xdata, ydata, sigma, best_par[0], best_par[1], best_par[2])
	N_res = 20
	alpha = np.linspace(7., 9., N_res)
	beta  = np.linspace(0.5e-3, 2e-3, N_res)
	gamma = np.linspace(5.5,6.5 , N_res)
	for i in range(N_res):
		for j in range(N_res):
			for k in range(N_res):
				tmp = chisq(xdata, ydata, sigma, alpha[i], beta[j], gamma[k])
				if(tmp < best_chi):
					best_chi = tmp
					best_par = np.array([alpha[i], beta[j], gamma[k]])
	return best_par  
	

Nsamp = 100
Npeak = 18

tmp  = np.loadtxt('dell.dat')
dell = np.reshape(tmp, (Nsamp, Npeak, 6))
stat_tt = np.zeros((Npeak, 3))
stat_ee = np.zeros((Npeak, 3))
stat_te = np.zeros((Npeak, 3))
xdata = np.zeros((Npeak*3,3))
ydata = np.zeros((Npeak*3,3))
sigma = np.zeros((Npeak*3,3))
	
for i in range(Npeak):
	ell_tt = dell[:, i, 0]
	dl_tt  = dell[:, i, 1]
	ell_ee = dell[:, i, 2]
	dl_ee  = dell[:, i, 3]
	ell_te = dell[:, i, 4]
	dl_te  = dell[:, i, 5]

	stat_tt[i,:] = meanstd(ell_tt, dl_tt) 
	stat_ee[i,:] = meanstd(ell_ee, dl_ee) 
	stat_te[i,:] = meanstd(ell_te, dl_te) 

stat  = np.concatenate((stat_tt, stat_ee, stat_te), axis=0)
xdata = stat[:,0]
ydata = stat[:,1]
sigma = stat[:,2]

par = bestfit(xdata, ydata, sigma)
par_tt = bestfit(stat_tt[:,0], stat_tt[:,1], stat_tt[:,2])
par_ee = bestfit(stat_ee[:,0], stat_ee[:,1], stat_ee[:,2])
par_te = bestfit(stat_te[:,0], stat_te[:,1], stat_te[:,2])

#print par
#print par_tt
#print par_ee
#print par_te

ell  = np.linspace(10, 3000, 100)
dl   = [fitfunc(x, par[0],    par[1],    par[2]) for x in ell]
dltt = [fitfunc(x, par_tt[0], par_tt[1], par_tt[2]) for x in ell]
dlee = [fitfunc(x, par_ee[0], par_ee[1], par_ee[2]) for x in ell]
dlte = [fitfunc(x, par_te[0], par_te[1], par_te[2]) for x in ell]

py.plot(ell, dltt, 'r', label = 'TT')
#py.plot(ell, dlee,'b', label = 'EE')
#py.plot(ell, dlte, 'k', label = 'TE')
#py.plot(ell, dl, label = 'XX')
py.errorbar(stat_tt[:,0], stat_tt[:,1], yerr= stat_tt[:,2], fmt='r ')
#py.errorbar(stat_ee[:,0], stat_ee[:,1], yerr= stat_ee[:,2], fmt='b ')
#py.errorbar(stat_te[:,0], stat_te[:,1], yerr= stat_te[:,2], fmt='k ')
py.legend(loc='best')
py.xlabel(r'$\ell$', fontsize=20)
py.ylabel(r'$f(\ell)$', fontsize=20)
py.show()
