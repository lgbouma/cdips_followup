from cdips_followup.spectools import remove_cosmic_rays
import numpy as np
import matplotlib.pyplot as plt

def main() -> None:
    """Demonstrate the cosmic ray removal routine using synthetic data."""
    np.random.seed(42)
    wavelength = np.linspace(4000, 5000, 1000)
    flux = (np.exp(-0.5 * ((wavelength - 4500) / 50)**2) +
            0.1 * np.random.normal(size=wavelength.size))
    num_cr = 10
    indices = np.random.choice(wavelength.size, num_cr, replace=False)
    flux[indices] += 5

    corrected_flux = remove_cosmic_rays(wavelength, flux)

    fig, ax = plt.subplots()
    ax.plot(wavelength, flux, label='Original Spectrum')
    ax.plot(wavelength, corrected_flux, label='Corrected Spectrum')
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Flux')
    ax.legend()
    fig.tight_layout()
    fig.savefig('spectrum_cosmic_ray_correction.png')


if __name__ == '__main__':
    main()

