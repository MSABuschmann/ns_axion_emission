import sys
import numpy as np
import h5py
from scipy.interpolate import interp1d 

from iminuit import Minuit

EOS = np.array(['APR_EOS_Cat_APR_Cat_1.4', 'BSk22_test_BSk22_1.4', 'BSk24_test_BSk24_1.4', 'BSk25_test_BSk25_1.4', 'BSk26_test_BSk26_1.4'])
EOS_short = np.array(['APR', 'BSk22', 'BSk24', 'BSk25', 'BSk26']) 
SF = np.array(['0-0-0', 'SFB-0-0', 'SFB-0-T73', 'SFB-a-T73', 'SFB-b-T73', 'SFB-c-T73'])

def Load(eos, sf):
    APR = h5py.File('../output/' + eos +'_' +sf+'.hdf5','r')

    params = APR['params'][:]
    E_raw = APR['E'][:]
    alpha = APR['alpha'][:]

    pbf_1s0n = np.zeros((len(alpha), len(E_raw)))
    pbf_1s0p = np.zeros((len(alpha), len(E_raw)))
    pbf_3p2A = np.zeros((len(alpha), len(E_raw)))
    pbf_3p2B = np.zeros((len(alpha), len(E_raw)))
    bremsstrahlung_nn = np.zeros((len(alpha), len(E_raw)))
    bremsstrahlung_np = np.zeros((len(alpha), len(E_raw)))
    bremsstrahlung_pp = np.zeros((len(alpha), len(E_raw)))

    for i in range(len(alpha)):
        pbf_1s0n[i] = APR['1s0n_a'+str(i)][:]
        pbf_1s0p[i] = APR['1s0p_a'+str(i)][:]
        pbf_3p2A[i] = APR['3p2A_a'+str(i)][:]
        pbf_3p2B[i] = APR['3p2B_a'+str(i)][:]
        bremsstrahlung_nn[i] = APR['bremsstrahlung_nn_a'+str(i)][:]
        bremsstrahlung_np[i] = APR['bremsstrahlung_np_a'+str(i)][:]
        bremsstrahlung_pp[i] = APR['bremsstrahlung_pp_a'+str(i)][:]

    total = pbf_1s0n + pbf_1s0p + pbf_3p2A + pbf_3p2B + bremsstrahlung_nn + bremsstrahlung_np + bremsstrahlung_pp

    return E_raw, alpha, total, params[1]

all_E = []
all_alpha = []
all_total = []
all_rNS = []

for eos in EOS:
    tmp_E = []
    tmp_alpha = []
    tmp_total = []
    tmp_rNS = []
    for sf in SF:
        E, alpha, total, rNS  = Load(eos, sf)
        tmp_E.append(E)
        tmp_alpha.append(alpha)
        tmp_total.append(total)
        tmp_rNS.append(rNS)
        
    all_E.append(np.array(tmp_E))
    all_alpha.append(np.array(tmp_alpha))
    all_total.append(np.array(tmp_total))
    all_rNS.append(np.array(tmp_rNS))

all_E = np.array(all_E)
all_alpha = np.array(all_alpha)
all_total = np.array(all_total)
all_rNS = np.array(all_rNS)

file = '../notebooks/Analysis_Results.npz'
data = np.load(file)

E_NuStar = data['E']
Eval_Flux = data['Eval_Flux']
Eval_TS = data['Eval_TS']

f_ct = np.zeros(len(E_NuStar))
f_up = np.zeros(len(E_NuStar))
f_dn = np.zeros(len(E_NuStar))

for i in range(len(E_NuStar)):
    locs = np.where( Eval_TS[i,:] <= 1 )
    c = np.argmin(Eval_TS[i,:])
    f_ct[i] = Eval_Flux[i,c] * 1e-15 # erg / (cm^2 s keV)
    f_dn[i] = Eval_Flux[i,locs[0][0]] * 1e-15 # erg / (cm^2 s keV)
    f_up[i] = Eval_Flux[i,locs[0][-1]] * 1e-15 # erg / (cm^2 s keV)
    
bin_width_NuStar = np.ones(len(E_NuStar)) * 2.5
bin_width_NuStar[0] = 1

data1856 = np.load('../notebooks/1856.npz', allow_pickle=True )

f_XMM = data1856['lower']*1e-15 # erg / (cm^2 s keV)
error_low_XMM = data1856['central']*1e-15 # erg / (cm^2 s keV)
error_high_XMM = data1856['upper']*1e-15 # erg / (cm^2 s keV)

E_XMM = np.array([3,5,7])
bin_width_XMM = np.array([1,1,1])

# Skip the first two bins of NuStar (only >= 10 keV)
SkipNustar = 2
NXMM = len(E_XMM)

E = np.empty(len(E_XMM) + len(E_NuStar) - SkipNustar)
E[:NXMM] = E_XMM
E[NXMM:] = E_NuStar[SkipNustar:]

bins_min = np.empty(len(E))
bins_max = np.empty(len(E))
bins_min[:NXMM] = E_XMM - bin_width_XMM
bins_max[:NXMM] = E_XMM + bin_width_XMM
bins_min[NXMM:] = E_NuStar[SkipNustar:] - bin_width_NuStar[SkipNustar:]
bins_max[NXMM:] = E_NuStar[SkipNustar:] + bin_width_NuStar[SkipNustar:]
dE = bins_max - bins_min
bin_width = dE

# Data flux
sig = np.empty(len(E))
sig[:NXMM] = f_XMM
sig[NXMM:] = f_ct[SkipNustar:]

# Gaussian Error
error = np.empty(len(E))
error[:NXMM] = (error_high_XMM - error_low_XMM)/2.
error[NXMM:] = (f_up[SkipNustar:] - f_dn[SkipNustar:])/2.

# truncate data
sig = sig[NXMM:]
error = error[NXMM:]
E = E[NXMM:]
dE = dE[NXMM:]
bins_max = bins_max[NXMM:]
bins_min = bins_min[NXMM:]
bin_width = bin_width[NXMM:]

# Priors
order=np.array(['RX J1856.5-3754'])
NSd=np.array(  [123.])
NSdd=np.array( [13.])
NSB0=np.array( [1.5e13])*0.95
NSdB0=np.array([0.1e13])
NSdB0 = np.sqrt(NSdB0**2+(0.1*NSB0)**2)
NSTeff_infty= np.array([0.050])
NSdTeff_infty= np.array([0.014])
NSdalpha = np.sqrt((NSdTeff_infty/NSTeff_infty * 4.* 0.455)**2 + 0.3**2 + (0.25* 4.* 0.455)**2)
NSdlogB0 = ( np.log10(NSB0+NSdB0) -np.log10(NSB0-NSdB0) )/2.
NSdlogalpha = ( np.log10(1+NSdalpha) -np.log10(1-NSdalpha) )/2.
Bpole = np.array([2.9])*1e13

binned_flux = np.zeros((len(EOS), len(SF), len(all_alpha[0][0]), len(bins_min)))
alpha_arr = all_alpha[0][0]
for eos in range(len(EOS)):
    for sf in range(len(SF)):
        for a in range(len(alpha_arr)):
            for i in range(len(bins_min)):
                whs = np.where( ( (all_E[eos,sf]>bins_min[i]) & (all_E[eos,sf] <= bins_max[i])  ))[0]
                binned_flux[eos,sf,a,i] = np.sum(all_total[eos,sf,a][whs]*(all_E[eos,sf][1]-all_E[eos,sf][0]))/dE[i]

## for conversion prob
alpha_EM = 1/137.
Bc = 4.41e13 #G
G_to_GeV2 = 1.95e-20 #GeV^2
cm_to_IGeV = 1/(1.98e-14)

# copied from our old code
class conv_prob:
    def __init__(self,rNS,B0,ma,gagg):
        '''
        rNS in cm
        B0 in Gauss
        ma in eV
        gagg in GeV^{-1}
        '''
        self.rNS = rNS*cm_to_IGeV #GeV^{-1}
        self.B0 = B0*G_to_GeV2 #GeV^2
        self.Bc = Bc*G_to_GeV2
        self.ma = ma*1e-9 #GeV
        self.D_M =1/2. *gagg *self.B0 # GeV
        
        #self._calc()
        
    def _calc(self,Es):
        '''
        Es in GeV
        '''
        self.x = (7*alpha_EM/45./np.pi)**(1/6.)*(Es/self.ma*self.B0/self.Bc)**(1/3.)
        self.D_a = - self.ma**2/2./Es
        self.rag = self.x*self.rNS
        
    def return_prob(self,omega):
        '''
        omega in keV
        '''
        Es = omega*1e-6 #GeV
        self._calc(Es)
        whs_0 = np.where(np.abs(self.D_a*self.rag)>0.45)[0]
        whs_1 = np.where(np.abs(self.D_a*self.rag)<0.45)[0]
        prob = np.zeros(len(omega))
        
        common = (self.D_M*self.rNS**3/self.rag**2)**2
        if len(whs_0)>0:
            prob[whs_0]= (common*(np.pi/3./np.abs(self.D_a*self.rag)) * np.exp(6./5.*self.D_a*self.rag) )[whs_0]
        if len(whs_1)>0:
            prob[whs_1]= (common*(2.21816**2/5**(6/5.)/np.abs(self.D_a*self.rag)**(4/5.))  )[whs_1]
        return prob

parsec2cm = 3.086e+18

# in erg/(cm^2 s keV)
def Flux(eos, sf, ma, gagg, gaNN, d, B, alpha):
    N = len(E)
    rNS = all_rNS[eos, sf] * 10. # cm
    cp = conv_prob(rNS, B, ma, gagg)
    probs = cp.return_prob(E)

    flux = np.zeros(N)
    for i in range(N):
        flux[i] = np.interp(alpha, alpha_arr, binned_flux[eos, sf, :, i])
        
    flux /= (parsec2cm * d)**2 * 4. * np.pi
    flux *= (gaNN/1e-10)**2
    flux *= probs
    
    return flux

def LL_gaus(Flux):
    return -(Flux-sig)**2/error**2

def return_LLNS(eos, sf, ma, gagg, d, theta, alpha):
    eos = int(eos)
    sf = int(sf)
    axion_flux = np.sign(ma) * Flux(eos=eos, sf=sf, ma=np.abs(ma), gagg=gagg, gaNN=1e-10, d=d, B=Bpole[0]*np.sin(theta), alpha=alpha)
    LL = 0.5*np.sum(LL_gaus(axion_flux))
    LL += - (NSd[0] - d)**2/(2.*NSdd[0]**2) + np.log(np.sin(theta))/3. - (alpha-1.)**2/(2.*NSdalpha[0]**2)
    return -LL

Nma = 200
Ngagg = 200

ma_arr=np.logspace(-6, -1, Nma)
gagg_arr=np.logspace(-11, -5, Ngagg)

def minimize(eos, sf, ma_start, gagg_start, d_start, theta_start, alpha_start, gi, params):
    m = Minuit(return_LLNS, eos=eos, sf=sf, ma=ma_start, gagg=gagg_start, d=d_start, theta=theta_start, alpha=alpha_start)

    m.fixed["ma"] = True
    m.fixed["gagg"] = True
    m.fixed["eos"] = True
    m.fixed["sf"] = True

    m.errors["d"] = NSdd[0]
    m.errors["theta"] = 0.1
    m.errors["alpha"] = NSdalpha[0]

    m.limits["d"] = (60,2000)
    m.limits["theta"] = (1e-10,np.pi-1e-10)
    m.limits["alpha"] = (0.1, 5.0)

    m.errordef = 0.5

    m.migrad()
    BF_ma = m.values["ma"]
    BF_gagg = m.values["gagg"]
    BF_d = m.values["d"]
    BF_theta = m.values["theta"]
    BF_alpha = m.values["alpha"]
    LL_local = return_LLNS(eos, sf, BF_ma, BF_gagg, BF_d, BF_theta, BF_alpha)

    if LL_local < params[gi,-1] or params[gi, -1] == 0:
        params[gi] = np.array([BF_ma, BF_gagg ,BF_d, BF_theta ,BF_alpha, LL_local])

all_params = []
eos = int(sys.argv[1])
sf = int(sys.argv[2])

for mi in range(Nma):
    print(mi)
    params = np.zeros((Ngagg, 6))
    for gi in range(Ngagg):
        minimize(eos, sf, ma_arr[mi], gagg_arr[gi], NSd[0],        np.pi/2.,    1.0,            gi, params)
        minimize(eos, sf, ma_arr[mi], gagg_arr[gi], NSd[0],        np.pi/2.,    1.0+NSdalpha[0],gi, params)
        minimize(eos, sf, ma_arr[mi], gagg_arr[gi], NSd[0],        np.pi/2.,    1.0-NSdalpha[0],gi, params)
        minimize(eos, sf, ma_arr[mi], gagg_arr[gi], NSd[0],        np.pi/2.-0.2,1.0,            gi, params)
        minimize(eos, sf, ma_arr[mi], gagg_arr[gi], NSd[0],        np.pi/2.-0.5,1.0,            gi, params)
        minimize(eos, sf, ma_arr[mi], gagg_arr[gi], NSd[0]+NSdd[0],np.pi/2.,    1.0,            gi, params)
        minimize(eos, sf, ma_arr[mi], gagg_arr[gi], NSd[0]-NSdd[0],np.pi/2.,    1.0,            gi, params)
    all_params.append(params)

all_params = np.array([all_params])

archive = h5py.File('results/eom_'+str(eos)+'_'+str(sf)+'.h5', 'w')
archive.create_dataset('data', data=all_params)
archive.close()
