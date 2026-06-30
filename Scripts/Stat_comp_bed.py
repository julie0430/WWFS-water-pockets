import numpy as np
import rasterio
import matplotlib.pyplot as plt
from scipy.stats import linregress, norm


def compare_bedrock_tifs(measured_tif, modelled_tif):
  
    # Read rasters
    with rasterio.open(measured_tif) as src:
        measured = src.read(1).astype(float)
        measured[measured == src.nodata] = np.nan

    with rasterio.open(modelled_tif) as src:
        modelled = src.read(1).astype(float)
        modelled[modelled == src.nodata] = np.nan

    # Keep valid pixels only
    mask = np.isfinite(measured) & np.isfinite(modelled)

    measured = measured[mask]
    modelled = modelled[mask]


    # Errors
    errors = modelled - measured

    # Statistics
    rmse = np.sqrt(np.mean(errors**2))
    mad = np.mean(np.abs(errors))
    bias = np.mean(errors)

    slope, intercept, r_value, _, _ = linregress(
        modelled,
        measured
    )

    r_squared = r_value**2

    print(f"RMSE : {rmse:.2f} m")
    print(f"MAD  : {mad:.2f} m")
    print(f"Bias : {bias:.2f} m")
    print(f"R²   : {r_squared:.3f}")

    # ==================================================
    # Scatter plot
    # ==================================================

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(modelled,measured,s=2,alpha=0.3)

    min_val = min(modelled.min(),measured.min())

    max_val = max(modelled.max(), measured.max())

    # 1:1 line
    ax.plot([min_val, max_val],[min_val, max_val],'k--',label='1:1')

    # Regression line
    x_line = np.linspace(min_val, max_val, 100)

    ax.plot(x_line,slope*x_line + intercept,'r',label=f'R² = {r_squared:.3f}')

    ax.text( 0.02, 0.98,
        f'RMSE = {rmse:.2f} m\n'
        f'MAD = {mad:.2f} m\n'
        f'Bias = {bias:.2f} m\n'
        f'R² = {r_squared:.3f}',
        transform=ax.transAxes,va='top')

    ax.set_xlabel("Bedrock modélisé [m]")
    ax.set_ylabel("Bedrock mesuré [m]")
    ax.legend()
    ax.grid(True)
    ax.set_aspect('equal')

    plt.title("Comparaison des bedrocks")

    plt.tight_layout()
    plt.show()

    # ==================================================
    # Histograms
    # ==================================================

    bins = 50
    # Absolute errors
    absolute_errors = np.abs(errors)
    absolute_errors = absolute_errors[np.isfinite(absolute_errors)]

    mu_abs = np.mean(absolute_errors)

    plt.figure(figsize=(8, 5))

    plt.hist(absolute_errors, bins=bins, alpha=0.6, density=True, edgecolor='black', label='Erreurs absolues')

    plt.axvline(mu_abs,color='red', linestyle='--',label=f'Moyenne = {mu_abs:.2f} m')

    plt.xlabel("Erreur absolue |Modélisé - Mesuré| [m]")
    plt.ylabel("Densité")
    plt.title("Distribution des erreurs absolues du bedrock")
    plt.legend()
    plt.xlim(0,300)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


    # Relative errors
    relative_errors = errors / measured * 100
    relative_errors = relative_errors[np.isfinite(relative_errors)]

    # Optional filter to remove extreme values
    relative_errors = relative_errors[
        (relative_errors > -90) & (relative_errors < 500)
    ]

    mu_rel, sigma_rel = norm.fit(relative_errors)

    plt.figure(figsize=(8, 5))

    print("Nombre de valeurs :", len(relative_errors))
    print("Nombre de valeurs uniques :", len(np.unique(relative_errors)))
    print("Min :", relative_errors.min())
    print("Max :", relative_errors.max())

    plt.hist(relative_errors,bins=bins, alpha=0.6,density=True,edgecolor='black',label='Erreurs relatives')

    x_rel = np.linspace(relative_errors.min(),relative_errors.max(), 100)

    plt.plot( x_rel,norm.pdf(x_rel, mu_rel, sigma_rel),linestyle='--',label=f'Loi normale, σ = {sigma_rel:.1f} %')

    plt.axvline(mu_rel,color='red',linestyle='--',label=f'Moyenne = {mu_rel:.1f} %')

    plt.xlabel("Erreur relative (Modélisé - Mesuré) / Mesuré [%]")
    plt.ylabel("Densité")
    plt.title("Distribution des erreurs relatives du bedrock")
    plt.legend()
    plt.xlim(-10,10)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    



compare_bedrock_tifs(
    "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Surf_BonnePierre_2024_11m.tif",
    "T:/RTM/06_ROGP/01_SUIVI_PAPROG/02-ACTIONS_EN_COURS_ONF/2026/4_stage-poches-eau/06_notes/WWFS/Input/Bed_BonnePierre_2025_11m.tif"
                    )     






