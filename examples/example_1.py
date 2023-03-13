from mcmphoton.classes import Simulation, Layer, TissueModel, Source, Point3D, Plotter

if __name__ == "__main__":
    # Case a:
    sim_size = 0.5
    sim = Simulation(sim_size=sim_size, multithreaded=True)

    # Simulation init
    # refractive_index, scattering_coefficient, absorption_coefficient, anisotropy_factor, z_range
    # If the absorption and scattering coefficient are very small the plots look odd.
    # Probably because the MFP becomes very long and photons just travel in a straight line without any interaction
    n = 1.1
    mu_s = 10
    mu_a = 1
    g = 0.7

    layer1 = Layer(n, mu_s, mu_a, g, (0, 0.5))

    tissue_model = TissueModel()
    tissue_model.add_layer(layer1)
    sim.add_tissue_model(tissue_model)

    source_position = Point3D(0, 0, -0.5)
    direction = Point3D(0, 0, 1)
    source = Source(sim, source_position, direction, number_photons=1000)
    sim.add_source(source)

    # Simulation run
    sim.start()
    sim.simulate(timing=True)

    # Plotting
    plotter = Plotter(sim)
    plotter.plot_photon_paths()
    plotter.plot_photon_paths(heatmap=True)