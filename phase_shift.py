
#Measure the phase shifts between two clxx


import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import camb4py

camb = camb4py.load('/home/zhenpan/camb/camb')

def run(As, ns, Tau, H0, omgb, omgc, nnu, Yp):
	result = camb(**{'get_scalar_cls':True, 'do_lensing':True, 'get_transfer':True, 'l_max_scalar': 6000,
		        'massless_neutrinos': nnu, 'massive_neutrinos':0, 'helium_fraction':Yp,
			'ombh2': omgb, 'omch2': omgc, 'hubble': H0, 're_optical_depth': Tau,
			'scalar_spectral_index(1)': ns,'scalar_amp(1)': As, 'pivot_scalar':0.05})
	Cls = result['scalar']
	Out = result['misc']
	return Cls, Out

def phase(phi, Delta):                                        #position of peaks and troughs
	N_inte = len(phi)
	xnew = phi
	ynew = Delta
	raw  = [i for i in range(1, N_inte-1) if (ynew[i]-ynew[i-1])*(ynew[i]-ynew[i+1]) > 0. ]
	posi = [raw[i] for i in range(1,len(raw)) if(raw[i]-raw[i-1] > 50)]
	phi_m   = refine(phi, Delta, posi)
	Delta_m = [ynew[i] for i in posi]
	return  phi_m, Delta_m

def refine(phi, Delta, posi):
	phi_m = np.zeros([len(posi)])
	for i in range(len(posi)):
		j = posi[i]
		xnew = [phi[j-1], phi[j], phi[j+1]]
		ynew = [Delta[j-1], Delta[j], Delta[j+1]]
		p = np.polyfit(xnew, ynew, 2)
		phi_m[i] = -1*p[1]/(2.*p[0])
	return phi_m
	

def analyse(As, ns, Tau, H0, omgb, omgc, nnu, Yp):           #return undamped spectrum and some backgroud viarables 
	Cls, Out = run(As, ns, Tau, H0, omgb, omgc, nnu, Yp)
	k_D = float(Out['k_D(zstar) Mpc'])
	r_s = float(Out['r_s(zstar)/Mpc'])
	theta_s = float(Out['100*theta'])
	theta_D = float(Out['100*theta_D'])
	zeq = float(Out[' '])
	keq = math.sqrt(2*(omgb+omgc)*(1+zeq))/3000.
	leq = keq*r_s/(theta_s/100)
	con = 0.65
	l_D = con*math.pi/(theta_D/100.) #k_D*r_s/(theta_s/100.)
	ell = [l for l in range(2,4000)]
	tt  = [Cls[l-2,1]*math.exp(2*(l/l_D)**1.18) for l in ell]
	ee  = [Cls[l-2,2]*math.exp(2*(l/l_D)**1.18) for l in ell]
	te  = [Cls[l-2,3]*math.exp(2*(l/l_D)**1.18) for l in ell]
	return ell, tt, ee, te, theta_s, theta_D, l_D, leq, zeq

def phase_shift(ell1, tt1, ell2, tt2):                     #phase shifts between two spectrums
	ell1_m, tt1_m = phase(ell1, tt1)
	ell2_m, tt2_m = phase(ell2, tt2)
	ell = [(ell1_m[i]+ell2_m[i])/2. for i in range(Npeak)]
	dell= [ell1_m[i]-ell2_m[i] for i in range(Npeak)]
	return ell, dell

Amp1, Ns1 , T_reio1, Hubble1, omgbh21, omgch21, Nnu1, Helium1 = np.loadtxt('params1a').T
Amp2, Ns2 , T_reio2, Hubble2, omgbh22, omgch22, Nnu2, Helium2 = np.loadtxt('params1b').T
Npeak = 18
Nsamp = 100
dell  = np.zeros((Nsamp, Npeak, 6)) 

for k in range(Nsamp):
	ell1, tt1, ee1, te1, ths1, thd1, ld1, leq1, zeq1 = analyse(Amp1[k], Ns1[k], T_reio1[k], Hubble1[k], omgbh21[k], omgch21[k], Nnu1[k], Helium1[k])
	ell2, tt2, ee2, te2, ths2, thd2, ld2, leq2, zeq2 = analyse(Amp2[k], Ns2[k], T_reio2[k], Hubble2[k], omgbh22[k], omgch22[k], Nnu2[k], Helium2[k])
	ell_tt, dell_tt = phase_shift(ell1, tt1, ell2, tt2)
	ell_ee, dell_ee = phase_shift(ell1, ee1, ell2, ee2)
	ell_te, dell_te = phase_shift(ell1, te1, ell2, te2)

	dell[k, :, 0] =  ell_tt
	dell[k, :, 1] = dell_tt
	dell[k, :, 2] =  ell_ee
	dell[k, :, 3] = dell_ee
	dell[k, :, 4] =  ell_te
	dell[k, :, 5] = dell_te
	
	tmp = np.reshape(dell,(Nsamp*Npeak, 6))
	np.savetxt('dell.dat', tmp)
 
	print k
	print  '100*theta_star:', ths1, ths2
	print  '100*theta_D:', thd1, thd2
	print  'z_eq:', zeq1, zeq2
	print  'l_D:', ld1,  ld2
	print 'l_eq:', leq1, leq2
	#plt.plot(ell_tt, dell_tt, 'o', label='TT')
	#plt.plot(ell_ee, dell_ee, 'o', label='EE')
	#plt.plot(ell_te, dell_te, 'o', label='TE')
	#plt.legend(loc = 'best')
	#plt.show()
