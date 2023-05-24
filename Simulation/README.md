# Simulation code
This code simulates pulsar emission using polar cap model. 
Absorption via magnetic field is calculated using an aproximated probability method based on escape energy of outgoing photons.
Simulation step is calculated from the energy that electrons lose, but sometimes this happens to be very large. Using an if condition one can modify this step in stellar radius units (we recommend 0.001 stellar radius for good results in light curves, but this may lead to large computational costs)
