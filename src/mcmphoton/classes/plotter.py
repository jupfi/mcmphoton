import matplotlib.pyplot as plt
import numpy as np


class Plotter():
    """
    A class to plot the results of a simulation.

    Attributes:
    simulation (Simulation): An instance of the simulation class.
    """

    def __init__(self, simulation):
        """
        The constructor for the Plotter class.

        Args:
            simulation (Simulation): An instance of the simulation class.
        """
        self.simulation = simulation

    def plot_photon_paths(self, heatmap=False, debug=False):
        """
        Plot the photon paths for a simulation.

        Args:
            heatmap (bool): If True, plots a heatmap of the photon absorption points.
            debug (bool): If True, marks the points where scattering or boundary crossing occurs for each photon.
                      Use this only for one photon since the plot gets very loaded otherwise.
        """
        fig, ax = plt.subplots(dpi=256)
        colors = ["green", "red", "blue", "yellow"]

        # Formatting
        ax.invert_yaxis()
        ax.set_xlim(-self.simulation.sim_size, self.simulation.sim_size)
        ax.set_xlabel("Radius in cm")
        ax.set_ylabel("z in cm")

        # Visualize Layers
        layers = self.simulation.tissue_model.layers
        for i, layer in enumerate(layers):
            # Skip visualization of first and last layers
            if i == (len(layers) - 1) or i == 0:
                continue
            ax.axhspan(layer.z_range[0], layer.z_range[1],
                       alpha=0.15, color=colors[i - 1 % len(colors)])

        sources = self.simulation.sources
        rs = []
        zs = []
        for i, source in enumerate(sources):
            # Visualize Sources
            r, z = source.get2D_position()
            ax.plot(r, z, "o", label="Source", color=colors[-i % len(colors)])

            # Plot photon paths
            photons = source.photons

            for photon in photons:
                points = photon.get2D_path_history()
                r = points[:, 0]
                z = points[:, 1]

                # Save all absorption points for every photon
                if heatmap:
                    rs.append(r[-1])
                    zs.append(z[-1])

                # Mark points when scattering or boundary crossing occurs for debugging
                elif debug:
                    ax.plot(r, z, color=colors[-i %
                            len(colors)], linewidth=0.1)
                    ax.plot(r, z, 'o', markersize=0.75,
                            color=colors[-i % len(colors)])

                # Plot photon paths and absorption point
                else:
                    ax.plot(r, z, color=colors[-i %
                            len(colors)], linewidth=0.1)
                    ax.plot(r[-1], z[-1], 'x', color=colors[-i % len(colors)])

            # Plot heatmap
            if heatmap:
                heatmap, xedges, yedges = np.histogram2d(rs, zs, bins=40)
                extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

                ax.imshow(heatmap.T, extent=extent, origin='lower')
                ax.invert_yaxis()
                ax.set_xlim(-self.simulation.sim_size,
                            self.simulation.sim_size)

        plt.show()
