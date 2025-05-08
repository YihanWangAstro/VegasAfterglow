from VegasAfterglow import Emission, ModelParams, Setups
import matplotlib.pyplot as plt
import numpy as np

Jy2cgs = 1e-23

# Define some custom profiles f(theta) or f(theta, phi)
def eps_k_profile(theta):
    return 1e52 * np.exp(-theta**2 / (2 * 0.088**2))

def Gamma0_profile(theta):
    return 300 * np.exp(-theta**2 / (2 * 0.35**2))

def rho_profile(theta, phi, r):
    return 1e-20 / r / r / r

def m_profile(r):
    return 1e20 * r * r * r / 3

# Initialize configuration
config = Setups()

# Choose a default medium type (ism, wind)
# config.medium = "ism"

# or set a custom one:
config.medium = "custom"
config.rho_profile = rho_profile 
config.m_profile = m_profile

# Choose a default jet type (tophat, powerlaw, gaussian)
config.jet = "gaussian"

# or set a custom one:
# config.jet = "custom"
# config.eps_k_profile = eps_k_profile 
# config.Gamma0_profile = Gamma0_profile 

# Define arbitrary engine duration (if unset, while be calculated)
# config.T0 = 10

# Jet spreading
config.spreading = False

# Model parameters
param = ModelParams()
param.E_iso = 1e52
param.theta_c = 0.2
param.theta_v = 0.8
param.p = 2.2
param.eps_e = 1e-2
param.eps_B = 1e-4
param.n_ism = 1e-2
param.xi = 1

# Initialize the Emission object
emission = Emission(param, config)

# Generate the light curve (lc)
t_list = np.geomspace(1e2, 1e8, 100)
nu = 1e16
lc = emission.lc(nu=nu, t=t_list)
lc_r = emission.lc_r(nu=nu, t=t_list)

# Generate the spectrum (spec)
t = 100
nu_list = np.geomspace(1e8, 1e18, 100)
spec = emission.spec(nu=nu_list, t=t)

fig, ax = plt.subplots(1, 2, figsize=(8, 4))
ax[0].loglog(t_list, np.array(lc) * nu * Jy2cgs, label='Forward Shock')
ax[0].loglog(t_list, np.array(lc_r) * nu * Jy2cgs, label='Reverse Shock')
ax[0].set_xlabel(r'$t$ [s]')
ax[0].set_ylabel(r'$F$ [cgs]')

ax[1].loglog(nu_list, np.array(spec) * nu * Jy2cgs)
ax[1].set_xlabel(r'$\nu$ [Hz]')
ax[1].set_ylabel(r'$F$ [cgs]')

ax[0].grid(True)
ax[1].grid(True)
ax[0].legend()

plt.tight_layout()
plt.savefig('emission.pdf')
plt.show()