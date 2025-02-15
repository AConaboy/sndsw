import numpy as np
import matplotlib.pyplot as plt

# Constants
fontsize=16
num_pixels = 14331  # Total number of pixels
pde = 0.5  # Photon Detection Efficiency (50%)
num_photons = np.linspace(0, 200000, 200, dtype=int)  # Photon counts from 0 to 1,000,000

# Function to simulate photon hits and saturation
def simulate_mppc_saturation(num_photons, num_pixels, pde):
    fired_pixels = []  # List to store the fired pixels for each nphoton burst
    
    for photons in num_photons:
        # Adjust for photon detection efficiency
        effective_photons = int(photons * pde)
        
        # Generate random hits for the effective photons across the pixels
        hits = np.random.randint(0, num_pixels, effective_photons)
        
        # Count unique pixels hit (saturation occurs when pixels are hit more than once)
        unique_pixels = np.unique(hits)
        fired_pixels.append(len(unique_pixels))
    
    return np.array(fired_pixels)

# Run the simulation
fired_pixels = simulate_mppc_saturation(num_photons, num_pixels, pde)

# Create the plot
fig, ax1 = plt.subplots(figsize=(8,6))

# Plot detected photons vs. total photons
ax1.tick_params(axis='both', which='major', labelsize=fontsize)
ax1.plot(num_photons, fired_pixels/num_pixels, marker='o', label="Fraction of pixels hit")
ax1.plot(num_photons, np.minimum(num_photons * pde / num_pixels, 1), 'r--', label=r"Ideal ($n_{\gamma} * PDE / n_{pixels}$)")
ax1.axhline(0.5, label='50% pixels fired', color='black', linestyle='--')

# Linear scale for the primary x-axis
# ax1.set_xlim(0, 1_000_000)
# ax1.set_ylim(0, num_pixels)
ax1.set_xlabel("Total photons fired", fontsize=fontsize)
ax1.set_ylabel("Fraction of pixels hit", fontsize=fontsize)

# Add secondary x-axis: Fraction of pixels hit
def photon_to_fraction(x):
    return (x / num_pixels)

def fraction_to_photon(x):
    return (x * num_pixels)

secay = ax1.secondary_yaxis('right', functions=(fraction_to_photon, photon_to_fraction))
secay.set_ylabel("Total fired pixels", fontsize=fontsize)
secay.tick_params(axis='y', which='major', labelsize=fontsize)
secay.grid(axis='y', linestyle="--", linewidth=0.5)

# Add title, legend, and grid
plt.title("Hamamatsu S14160 6050HS saturation curve"+"\n PDE("+r"$\lambda = 450$ nm) = 50%", fontsize=fontsize)
ax1.legend(fontsize=fontsize)
ax1.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.tight_layout()
plt.savefig('./saturationsimulation.png')
