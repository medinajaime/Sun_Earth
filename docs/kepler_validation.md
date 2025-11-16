# Solar-System Orbit Validation

This project now simulates any subset of the Sun, all eight planets, and Earth’s Moon with real orbital elements. Each body is propagated in full SI units with a velocity-Verlet integrator, and post-processing checks every Kepler law per body to quantify numerical accuracy. Parents are treated as fixed two-body partners (e.g., each planet feels only the Sun, while the Moon feels Earth), matching the project requirement to “add one object at a time.”

## Physical Data

| Quantity | Value | Source |
| --- | --- | --- |
| Gravitational constant \(G\) | \(6.67430\times 10^{-11}\ \text{m}^3\text{kg}^{-1}\text{s}^{-2}\) | CODATA 2018 [[1]](#1) |
| Nominal solar mass \(M_\odot\) | \(1.98847\times 10^{30}\ \text{kg}\) | IAU 2015 Resolution B3 [[2]](#2) |
| Nominal solar radius \(R_\odot\) | \(6.957\times 10^{8}\ \text{m}\) | IAU 2015 B3 [[2]](#2) |
| Astronomical Unit | \(149\,597\,870\,700\ \text{m}\) (exact) | IAU 2012 Resolution B2 [[3]](#3) |

### Orbital Catalog (epoch J2000.0, ecliptic frame)

| Body | Mass (kg) | Radius (km) | \(a\) (AU) | \(e\) | \(i\) (deg) | \(\Omega\) (deg) | \(\omega\) (deg) | \(M\) (deg) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Mercury | \(3.3011\times10^{23}\) | 2 439.7 | 0.387099 | 0.20563 | 7.0049 | 48.3317 | 29.1248 | 174.7944 |
| Venus | \(4.8675\times10^{24}\) | 6 051.8 | 0.723332 | 0.00677 | 3.3947 | 76.6807 | 54.8523 | 50.4468 |
| Earth | \(5.9722\times10^{24}\) | 6 371.0 | 1.000000 | 0.01671 | 0.00005 | −11.2606 | 114.2078 | 357.5172 |
| Moon | \(7.3420\times10^{22}\) | 1 737.4 | 0.00257 | 0.05490 | 5.1450 | 125.08 | 318.15 | 135.09 |
| Mars | \(6.4171\times10^{23}\) | 3 389.5 | 1.523662 | 0.09341 | 1.8506 | 49.5785 | 286.4623 | 19.4125 |
| Jupiter | \(1.8982\times10^{27}\) | 69 911 | 5.203363 | 0.04839 | 1.3053 | 100.5562 | 274.1977 | 19.6505 |
| Saturn | \(5.6834\times10^{26}\) | 58 232 | 9.537070 | 0.05415 | 2.4845 | 113.7150 | 338.7169 | 317.5124 |
| Uranus | \(8.6810\times10^{25}\) | 25 362 | 19.191264 | 0.04717 | 0.7699 | 74.2299 | 96.7344 | 142.2679 |
| Neptune | \(1.02413\times10^{26}\) | 24 622 | 30.068963 | 0.00859 | 1.7692 | 131.7217 | 273.2497 | 259.9087 |

Planetary values originate from the NASA/JPL Planetary Fact Sheet (J2000 mean orbital elements) [[4]](#4). The lunar elements come from the NASA Moon Fact Sheet and IAU recommendations [[5]](#5). All semi-major axes are expressed in AU; the simulator converts them to meters internally.

## Simulation Highlights

- Multi-body support: choose any combination of planets and the Moon with `--bodies`; parents are added automatically so you can literally add one object at a time.
- Full 3D orientation respecting inclination, longitude of ascending node, and argument of periapsis for each orbit.
- Per-body state tracking (position, velocity, specific energy, swept-area rate, true anomaly) logged to `data/<body>_orbit_log.csv` plus `<body>_orbit_log.summary.json`.
- Global `data/summary.json` aggregates all diagnostics for convenient report ingestion.
- VPython canvas renders one emissive Sun plus textured/smooth planets, orbital trails, and labels; distance/energy graphs allocate separate traces per body. Use `--headless` for CI or notebook contexts.
- Flexible CLI: `--dt-minutes`, `--days`, `--segment-hours`, `--log-dir`, `--bodies`, `--list-bodies`, and `--headless`.

## Kepler-Law Verification

Each body receives its own diagnostics:

| Field | Meaning |
| --- | --- |
| `k1_max_rel_error_pct` | Maximum relative deviation from the analytic ellipse \(r=a(1-e^2)/(1+e\cos\theta)\). |
| `k2_area_std_pct` | Relative std. dev. of swept area per `segment-hours` window (should be ~0). |
| `k3_period_days` | Simulated orbital period recovered from true-anomaly wrap. |
| `k3_ratio_pct_err` | Percent error of \(T^2/a^3\) vs. the theoretical \(4\pi^2/\mu\). |
| `energy_drift_pct` | Percent drift of specific orbital energy over the entire run. |

Running `python solar_orbit.py --headless --bodies earth --days 370` now produces the same sub-millimeter agreement reported previously, plus JSON summaries for any additional bodies you add.

## How to Run

1. Install dependencies (VPython is the only runtime requirement):
   ```
   python -m pip install -r requirements.txt
   ```
2. Inspect available bodies (ids are case-insensitive):
   ```
   python solar_orbit.py --list-bodies
   ```
3. Visualize a growing solar system (add one object at a time):
   ```
   python solar_orbit.py --bodies earth
   python solar_orbit.py --bodies earth mars
   python solar_orbit.py --bodies mercury venus earth moon
   ```
4. Headless run with tighter resolution and custom logs:
   ```
   python solar_orbit.py --headless --bodies earth moon mars \
     --dt-minutes 5 --segment-hours 3 --days 400 --log-dir data
   ```
5. Inspect outputs in `data/`:
   - `<body>_orbit_log.csv` – full history for that body.
   - `<body>_orbit_log.summary.json` – Kepler metrics for that body.
   - `summary.json` – combined report keyed by body id.

To validate convergence, rerun with progressively smaller `--dt-minutes` and compare the JSON summaries. Because each body logs its own state, you can also ingest the CSV files in pandas, MATLAB, or spreadsheets to build plots or compute additional metrics.

## References

1. <a id="1"></a>Mohr, P.J. et al., “CODATA 2018 Recommended Values of the Fundamental Constants,” NIST SP 959 (2019).
2. <a id="2"></a>Prša, A. et al., “IAU 2015 Resolution B3 on nominal solar and planetary values,” *Astronomical Journal* 152, 41 (2016).
3. <a id="3"></a>International Astronomical Union, “Resolution B2 on the re-definition of the astronomical unit of length,” 2012.
4. <a id="4"></a>NASA Goddard Space Flight Center, “Planetary Fact Sheet – Metric,” https://ssd.jpl.nasa.gov/planets/phys_par.html
5. <a id="5"></a>NASA GSFC, “Moon Fact Sheet,” https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
