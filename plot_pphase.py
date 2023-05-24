def plot_pphase(file, phase=None, energy=None, min_energy=None, max_energy=None, bin_size=0.001953125):
    '''
    :param file:
    :param phase:
    :param energy:
    :param min_energy:
    :param max_energy:
    :param bin_size:
    :return:
    '''
    # Importamos las librerias necesarias
    from astropy.io import fits
    from glob import glob
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    if phase == None or energy == None: # Si no se aportan las listas de energía y fase se obtienen del archivo
        archivo_fits = fits.open(file)  # open file
        data_list = archivo_fits[1].data
        energy = []
        phase = []
        pbar1 = tqdm(range(len(data_list)))
        for m in pbar1:
            energy.append((data_list[m])[0])
            phase.append((data_list[m])[-1])  # la pulse_phase debe estar al final de la lista
            pbar1.set_description(f'Extracting data')

    if min_energy == None:
        min_energy = min(energy)
        print('Min energy is', min_energy)
    if max_energy == None:
        max_energy = max(energy)
        print('Max energy is', max_energy)

    ## Realizamos los cortes en la fase según la energía indicada
    cut_phase = []
    pbar2 = tqdm(range(len(energy)))
    for m in pbar2:
        if float(min_energy) <= energy[m] <= float(max_energy):
            cut_phase.append(phase[m])
            pbar2.set_description(f'Applying energy cuts')

    # Pintamos el resultado
    n_bin = int(np.floor(1 / bin_size))
    plt.figure(1)
    (counts, bins, patches) = plt.hist(cut_phase, bins=n_bin, histtype='step', color='k')
    print('Number of photons in this energy range is', sum(counts))
    plt.grid()
    plt.xlabel("Pulse Phase")
    plt.ylabel("Counts")
    plt.ylim([min(counts)*0.9, max(counts)*1.1])
    plt.title('Pulsar pulse phase for' + ' ' + str(np.floor(min_energy)) + ' MeV < E <' + ' ' + str(np.floor(max_energy)) + ' ' +'MeV')

    phase2 = []

    for m in range(len(bins) - 1):
        phase2.append((bins[m] + bins[m + 1]) / 2)

    idx = np.argmax(counts)
    maxphase1 = phase2[idx]

    print(maxphase1)

    plt.figure(2)
    plt.plot(phase2, counts, 'k')
    plt.plot(phase2+np.ones(len(phase2)), counts, 'k')
    plt.grid()
    plt.xlabel("Pulse Phase", fontsize=12)
    plt.ylabel("Counts",fontsize=12)
    plt.xlim([0, 2])
    plt.ylim([min(counts) * 0.9, max(counts) * 1.1])
    plt.title('Vela pulsar light curve for E > 0.1 GeV', fontsize=18)

    plt.show()


plot_pphase('vela_merged2.fits', min_energy=100)
