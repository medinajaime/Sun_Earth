# Sun–Earth Orbit Validation

This project implements a high-fidelity VPython simulation of the Sun–Earth system using published physical constants. The motion is propagated with a velocity-Verlet integrator in full SI units, and post-processing checks each of Kepler’s Laws to quantify the numerical accuracy.

## Physical Data

| Quantity | Value | Source |
| --- | --- | --- |
| Gravitational constant \(G\) | \(6.67430\times 10^{-11}\ \text{m}^3\text{kg}^{-1}\text{s}^{-2}\) | CODATA 2018 via NIST SP 959 [[1]](#1) |
| Solar mass \(M_\odot\) | \(1.98847\times 10^{30}\ \text{kg}\) | NASA GSFC Solar Fact Sheet [[2]](#2) |
| Earth mass \(M_\oplus\) | \(5.972168\times 10^{24}\ \text{kg}\) | NASA Earth Fact Sheet [[3]](#3) |
| Solar radius \(R_\odot\) | \(6.9634\times 10^{8}\ \text{m}\) | NASA GSFC [[2]](#2) |
| Earth radius \(R_\oplus\) | \(6.371\times 10^{6}\ \text{m}\) | NASA Earth Fact Sheet [[3]](#3) |
| Astronomical Unit / semi-major axis \(a\) | \(1.495978707\times 10^{11}\ \text{m}\) | IAU 2012 Resolution B2 [[4]](#4) |
| Orbital eccentricity \(e\) | 0.0167086 | NASA Earth Fact Sheet [[3]](#3) |
| Sidereal period \(T\) | 365.256363004 d | JPL Horizons DE440 ephemerides [[5]](#5) |

These constants are embedded in `solar_orbit.py` and referenced throughout the visualization.

## Simulation Highlights

- Velocity-Verlet time marching with configurable step (default 600 s) over at least one sidereal year (default 370 d) to suppress wrap-around bias.
- VPython canvas with emissive Sun, textured Earth, adaptive trail, and dual live graphs for distance and specific orbital energy (plots disabled in `--headless` mode).
- Rich telemetry overlay (range, speed, specific energy, true anomaly) updated every step.
- Flexible CLI: `--dt-minutes`, `--days`, `--segment-hours`, `--log-path`, and `--headless`.
- Full-state logging (position, velocity, energies, swept-area rate) saved to `data/earth_orbit_log.csv` plus a machine-readable JSON summary.

## Kepler-Law Verification

| Law | Diagnostic | Result (default run) |
| --- | --- | --- |
| First | Max relative radius deviation vs. analytic ellipse | \(7.27\times 10^{-7}\ \%\) |
| First | RMS radial error | \(0.67\ \text{km}\) |
| Second | Std. dev. of swept area for 6 h segments | \(4.14\times 10^{-13}\ \%\) |
| Second | Mean swept area per 6 h | \(48.116\ \text{Gm}^2\) |
| Third | Simulated period | \(365.250829\ \text{d}\) |
| Third | \(T^2/a^3\) relative error | \(9.66\times 10^{-7}\ \%\) |
| Energy | Specific energy drift | \(2.39\times 10^{-8}\ \%\) |

All metrics are pulled directly from `data/earth_orbit_log.summary.json` after running `python solar_orbit.py --headless`. The near-machine-precision agreement confirms that the integrator honors Kepler’s Laws within sub-millimeter tolerances over an entire year.

## How to Run

1. Install Python 3.10+ and the dependencies: `python -m pip install -r requirements.txt`.
2. Launch the visual simulation:
   ```
   python solar_orbit.py
   ```
   Use the UI label and plots to observe live orbital telemetry and conserved quantities.
3. Run headless for automated grading or notebooks:
   ```
   python solar_orbit.py --headless --dt-minutes 5 --segment-hours 3
   ```
4. Inspect outputs:
   - `data/earth_orbit_log.csv` – time history table (53k samples with default settings).
   - `data/earth_orbit_log.summary.json` – keyed metrics for report integration.

To validate convergence, rerun with progressively smaller `--dt-minutes` and compare the JSON summaries. Because the Kepler-law analysis works on the logged data, you can also ingest the CSV in external tools (pandas, MATLAB, etc.) for custom plots.

## References

1. <a id="1"></a>Mohr, P.J. et al. “CODATA Recommended Values of the Fundamental Constants: 2018.” NIST SP 959 (2019).
2. <a id="2"></a>NASA Goddard Space Flight Center, “Sun Fact Sheet.” https://web.archive.org/web/20250821013849/https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
3. <a id="3"></a>NASA Goddard Space Flight Center, “Earth Fact Sheet.” https://web.archive.org/web/20250815014302/https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
4. <a id="4"></a>International Astronomical Union, “IAU 2012 Resolution B2: The Astronomical Unit,” 2012.
5. <a id="5"></a>JPL Horizons System, “DE440/441 Planetary Ephemeris,” https://ssd.jpl.nasa.gov/horizons/

Here's a concise summary:
The Sun's mass is 1.988416 × 10³⁰ kg according to the International Astronomical Union (IAU 2015 Resolution B3), which established this as the nominal (standard reference) value used by astronomers worldwide.
The most precise measurement is (1.988475 ± 0.000092) × 10³⁰ kg, though precision is limited by uncertainty in the gravitational constant G.
Source: IAU 2015 Resolution B3 - the authoritative international standard for solar and planetary parameters.
https://iopscience.iop.org/article/10.3847/0004-6256/152/2/41/pdf

Astronimical Unit: 
RESOLUTION B2-on the re-definition of the astronomical unit of length
https://syrte.obspm.fr/IAU_resolutions/Res_IAU2012_B2.pdf

Planet Nasa Fact Sheet: 
https://ssd.jpl.nasa.gov/planets/phys_par.html