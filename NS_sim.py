
def ns_sim(p, alpha, rs, b_0, e_in_en, g_in, r_in, theta_in, en_loss, xi_obs=None, att='False'):
    '''
    :param p: Period in s of the pulsar
    :param alpha: Angle between magnetic and rotation axis
    :param rs: Stellar radio
    :param b_0: Magnetic field in star surface
    :param e_in_en: Initial electron energy
    :param g_in: Initial lorentz factor for electrons
    :param r_in: Initial radial position of the electron
    :param theta_in: Initial colatitude for the electron
    :param en_loss: Energy that electrons lose each step
    :param xi_obs: Observation angle in galactic coordinates
    :param att: Boolean value which allows to visualize attenuation effects in each step
    :return: sn: summed spectrum for each step, f_en: photon energies (may need to be adjusted manually depending on
    initial conditions), sc: sum of binned spectrum for each step, binned_en: binned photon energies, photon_number: sum
    of the emitted photons in each step, phi: phase in which photons are emitted in each step, xi: galactical polar
    coordinate of the emitted photons in each step.
    '''
    # We import our libraries
    import numpy as np
    import scipy.special as sp
    from scipy.integrate import quad
    # We define our constants
    e = 1.6e-19  # C
    me = 9.1e-31  # kg
    mu = 4*np.pi*1e-7  # H/m
    c = 299792458  # m/s
    h = 4.135667696*1e-15  # eV/s
    hbar = h/(2*np.pi)

    if en_loss >= 0:
        return print('Energy loss rate must be negative')
    # We calculate geometry parameters
    rlc = c*p/(2*np.pi)  # Radio of the light cylinder
    theta_pc = np.arcsin(np.sqrt(rs/rlc))  # Polar cap angle
    rpc = theta_pc*rs  # Polar cap radius
    
    def curv_r(r, theta):  # Curvature radius function
        k = r/(np.sin(theta))**2
        return k*((np.sin(theta))**4+(np.sin(2*theta))**2)/((np.sin(theta))**4+2*(np.sin(2*theta))**2-2*(np.sin(theta))**2*np.cos(2*theta))

    def curv_tan(theta):  # Local tangent function
        return (3*(np.cos(theta))**2-1)/(3*np.sin(theta)*np.cos(theta))

    # Electron acceleration
    def E_par(height):
        return 1e9*np.exp(-height/rpc)

    # Curvature radiation maximum
    def Ecr(gamma, r, theta):
        return 3/2*c*hbar*gamma**3/curv_r(r, theta)

    def Fsyn(x):
        return x * quad(lambda y: sp.kv(5 / 3, y), x, np.inf)[0]

    # Curvature radiation
    def dNdE(E, gamma, r, theta):
        return (np.sqrt(3)*e**2*gamma)/(h*E*curv_r(r, theta))*Fsyn(E/Ecr(gamma, r, theta))

    # Energy loss rate
    def gdot(gamma, r, theta):
        return -2/3*e**2/(me*c)*(gamma**4)/((curv_r(r, theta))**2)

    # Absorption via magnetic field
    def B(m, r, theta):  # theta must be in rad
        return mu*m/(4*np.pi*r**3)*np.sqrt(1+3*(np.sin(np.pi-theta))**2)
    mag = b_0 * rs**3 * 4*np.pi/(mu * np.sqrt(1+3*(np.sin(np.pi-theta_in))**2))  # we calculate the intensity of
    # magnetic dipole given B on the NS surface
    bcrit = 4.4 * 10**13  # in gauss

    def e_esc(r, theta):
        return (0.1 * bcrit * 0.01) / (theta * b_0) * (r / rs)**(5/2)
    # The algorithm: given initial E and position (r, theta)
    en = np.copy(e_in_en)
    e_r = np.copy(r_in)
    e_theta = np.copy(theta_pc)
    k = e_r/(np.sin(e_theta))**2

    if en == None:
        g = g_in
        if g == None:
            return print('No initial electron energy or Lorentz factor given')
    else:
        g = en/0.510998950  # electron energy must be in MeV
        g_in = np.copy(g)
    f_en = np.logspace(-3, 4, 1000)  # we initialize photon energies
    sn = np.zeros(len(f_en))  # 0's list in which we store the summed photon spectrum
    hist = np.histogram(f_en, bins=f_en)  # we bin the data
    bins = hist[1]
    binned_en = []
    photon_number = []
    phi=[]
    xi=[]
    beta = np.linspace(0, 2*np.pi, 10000)
    for n in range(len(bins) - 1):  # this can be done just once since f_en is the same in each step
        binned_en.append((bins[n] + bins[n + 1]) / 2)
    sc = np.zeros(len(binned_en))  # 0's list in which store counts for binned spectrum
    while e_r < rlc:  # we track the electron until it reaches the light cylinder
        vn = []  # empty list in which we store spectrum in each step
        for n in f_en:  # we calculate photon spectrum
            vn.append(dNdE(n*10**9, g, e_r, e_theta)*10**9)  # This gives dNdE in eV^-1, so we multiply by 10^9 to
            # get GeV^-1)
        # vn = vn / vn[-1]  # we normalize the counts
        hist = np.histogram(f_en, bins=f_en, weights=vn)  # we bin the data
        counts = hist[0]
        v = np.copy(counts)
        esc = e_esc(e_r, e_theta)
        print('Escape energy in this step:', '{:.2f}'.format(esc))
        # we calculate attenuation for each spectrum
        for n in range(len(binned_en)-1):
            long_bin = (binned_en[n+1] - binned_en[n]) / 2 #semi-size of each bin
            if binned_en[n]-long_bin < esc < binned_en[n]+long_bin:
                p_abs = esc/(binned_en[n]+long_bin)  # probability to be absorbed
                counts[n] = counts[n]*(1-p_abs)  # number of photons after the absorption
                for d in range(n+1, len(binned_en)):
                    counts[d] = 0
        sc += counts
        sn += vn
        z = np.copy(counts)
        if att == 'True':
            for n in range(len(binned_en)):
                v[n] = v[n] * binned_en[n] ** 2
                z[n] = z[n]*binned_en[n]**2
            plt.plot(binned_en, v, label='No absorption')
            plt.plot(binned_en, z, '-.', label='With absorption')
            plt.grid(True)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Photon energy E [GeV]')
            plt.ylabel(r'$E^2 dN/dE [MeV s^{-1}]$')
            plt.legend()
            plt.show()
        if xi_obs == None:
            for n in range(len(beta)):
                photon_number.append(np.sum(counts))  # this number of photons is emitted in the following direction
                xi.append(np.arccos((np.cos(alpha)*np.cos(e_theta)+np.sin(alpha)*np.cos(e_theta)*np.cos(np.pi-beta[n]))))
                phi.append(np.arccos((np.cos(e_theta)-np.cos(alpha)*np.cos(xi[n]))/(np.sin(alpha)*np.sin(xi[n]))))
        if xi_obs != None:
            phi.append(np.arccos((np.cos(e_theta)-np.cos(alpha)*np.cos(xi_obs))/(np.sin(alpha)*np.sin(xi_obs))))
            photon_number.append(sum(counts))

        ds = en_loss * c * g / gdot(g, e_r, e_theta)  # step for fixed energy loos rate
        if ds > 0.01*rs:
            ds = 0.01*rs
        dr = ds * (1+(np.tan(e_theta))**2/4)**(-1/2)  # step for radius
        e_r += dr  # we update the parameters to the next step
        e_theta = np.arcsin(np.sqrt(e_r/k))
        g += ds*gdot(g, e_r, e_theta)/c
        print('Progress:', '{:.2f}'.format(e_r/rlc*100), '%')
    return sn, f_en, sc, binned_en, photon_number, phi, xi

import numpy as np
import matplotlib.pyplot as plt


summed, photon_en, sc, binned_en, photon_number, phi, xi = ns_sim(p=2.3*10**(-3), xi_obs=65*np.pi/180, alpha=60*np.pi/180, rs=10e3, b_0=10**9, e_in_en=None, g_in=2e7, r_in=10e3, theta_in=15*np.pi/180, en_loss=-0.01)
hist2 = np.histogram(photon_en, 100, weights=summed)
counts = hist2[0]
b_photon_en = hist2[1]

scaled_sc = []
scaled_counts = []

for n in range(len(binned_en)):
    scaled_sc.append(binned_en[n]**2 * 10**3 * sc[n])

plt.figure(1)
plt.grid(True)
plt.plot(binned_en, scaled_sc, '-k')
plt.title('Summed emitted photon spectrum')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Photon energy E [GeV]')
plt.ylabel(r'$E^2 dN/dE [MeV s^{-1}]$')
