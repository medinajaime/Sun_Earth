"""
Sun–Earth orbital simulation with VPython visualization and Kepler-law diagnostics.

Physical constants sourced from:
- CODATA 2018 gravitational constant via NIST SP 959 (G = 6.67430e-11 m^3 kg^-1 s^-2).
- NASA GSFC Solar System Exploration: solar mass 1.98847e30 kg, solar radius 6.9634e8 m.
- NASA Earth Fact Sheet: Earth mass 5.972168e24 kg, Earth radius 6.371e6 m, orbital eccentricity 0.0167086.
- IAU 2012 resolution B2: Astronomical Unit (AU) = 149,597,870,700 m and Earth's semi-major axis 1.495978707e11 m.
- JPL Horizons (DE440): Sidereal orbital period 365.256363004 days.

Run `python solar_orbit.py` for full VPython experience or pass `--headless` to skip
graphics while still generating logs and metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from vpython import arrow, canvas, color, curve, graph, gcurve, label, rate, sphere, vector, cross, mag, norm


# --- Fundamental constants (SI units) ---
G = 6.67430e-11  # CODATA 2018
M_SUN = 1.988416e30  # kg, NASA GSFC
M_EARTH = 5.972168e24  # kg, NASA Earth Fact Sheet
SUN_RADIUS = 6.957e8  # m, NASA
EARTH_RADIUS = 6.371e6  # m, NASA
AU = 1.495978707e11  # m, IAU 2012
EARTH_ECCENTRICITY = 0.01671  # NASA / JPL
EARTH_SEMI_MAJOR_AXIS = AU  # m
SIDEREAL_PERIOD_DAYS = 365.256363004  # JPL Horizons

MU = G * (M_SUN + M_EARTH)  # Standard gravitational parameter (m^3 / s^2)


@dataclass
class SimulationConfig:
    dt: float = 600.0  # seconds
    total_days: float = 370.0
    log_path: Path = Path("data/earth_orbit_log.csv")
    segment_hours: float = 6.0

    @property
    def total_time(self) -> float:
        return self.total_days * 86400.0

    @property
    def segment_duration(self) -> float:
        return self.segment_hours * 3600.0


class SunEarthSimulation:
    """VPython-based velocity-Verlet integrator for the Sun–Earth two-body problem."""

    def __init__(self, config: SimulationConfig, visual: bool = True) -> None:
        self.config = config
        self.visual = visual
        self.time_s = 0.0
        self.display_scale = AU  # meters per VPython unit
        self.records: List[Dict[str, float]] = []

        self.position = vector(EARTH_SEMI_MAJOR_AXIS * (1 - EARTH_ECCENTRICITY), 0, 0)
        v_peri = math.sqrt(
            MU * (1 + EARTH_ECCENTRICITY) / (EARTH_SEMI_MAJOR_AXIS * (1 - EARTH_ECCENTRICITY))
        )
        self.velocity = vector(0, v_peri, 0)

        self.scene = None
        #self.star_outer = None
        #self.star_inner = None
        #self.starfield = None
        self.sun = None
        self.earth = None
        self.orbit_trail = None
        self.distance_curve = None
        self.energy_curve = None
        self.stats_label = None
        self.radial_line = None
        self.velocity_arrow = None
        self.peri_marker = None
        self.aphelion_marker = None
        self.peri_line = None
        self.aphelion_line = None

        if self.visual:
            self._setup_scene()

    # --- Scene helpers -----------------------------------------------------
    def _setup_scene(self) -> None:
        self.scene = canvas(
            title="Sun–Earth Orbit (Real Units, Velocity-Verlet)",
            width=1000,
            height=700,
            background=color.black,
        )

        sun_radius_vis = max(0.2, SUN_RADIUS / self.display_scale * 30)
        earth_radius_vis = max(0.04, EARTH_RADIUS / self.display_scale * 500)
        star_texture = "https://i.postimg.cc/VYw7JrfL/2k-stars-milky-way.jpg"
        self.scene.forward
        background_dir = -norm(self.scene.forward)
        BACKGROUND_DISTANCE = 200 * AU    # far behind everything
        BACKGROUND_RADIUS   = 180 * AU    # large enough to cover FOV


        
        """self.star_outer = sphere(
            pos=vector(0,0,0),
            radius=1.25,
            texture=star_texture,
            emissive=True,
            shininess=0,
        )

        self.star_inner = sphere(
            pos=vector(0,0,0),
            radius=1.249,
            texture=star_texture,
            emissive=True,
            shininess=0,
            opacity=0.9999
        )

        self.starfield = sphere(
            pos = background_dir * BACKGROUND_DISTANCE,
            radius = BACKGROUND_RADIUS,
            texture = "https://i.postimg.cc/VYw7JrfL/2k-stars-milky-way.jpg",
            shininess = 0,
            emissive = True,
            opacity = 1
        )
        """



        self.sun = sphere(
            pos=vector(0, 0, 0), 
            radius=sun_radius_vis, 
            color=color.orange, 
            emissive=True,
            texture="https://i.postimg.cc/WjBQw4SW-/2k-sun.jpg")
            
        self.earth = sphere(
            pos=self._to_visual(self.position),
            radius=earth_radius_vis,
            make_trail=False,
            texture="https://i.postimg.cc/ncRX2TcM/2k-earth-daymap.jpg",
        )


        self.orbit_trail = curve(color=color.cyan, radius=0.01)
        self.radial_line = curve(color=color.white, radius=0.004)
        self.radial_line.append(pos=vector(0, 0, 0))
        self.radial_line.append(pos=self._to_visual(self.position))
        self.velocity_arrow = arrow(
            pos=self._to_visual(self.position),
            axis=self._velocity_axis(self.velocity),
            color=color.green,
            shaftwidth=0.02,
            opacity=0.9,
        )
        self._setup_orbital_markers(earth_radius_vis)

        scene_range = 1.3  # in AU units
        self.scene.range = scene_range
        self.scene.forward = vector(-1, -0.4, -1.3)

        dist_graph = graph(
            title="Earth–Sun Distance",
            xtitle="Days",
            ytitle="Distance (AU)",
            fast=False,
            align="left",
            width=600,
            height=250,
        )
        energy_graph = graph(
            title="Specific Orbital Energy",
            xtitle="Days",
            ytitle="Energy (MJ/kg)",
            fast=False,
            align="right",
            width=600,
            height=250,
        )
        self.distance_curve = gcurve(color=color.cyan, graph=dist_graph)
        self.energy_curve = gcurve(color=color.yellow, graph=energy_graph)
        self.stats_label = label(
            pos=vector(-scene_range * 0.8, scene_range * 1.1, 0),
            text="Initializing...",
            color=color.white,
            box=False,
            height=16,
        )

    def _to_visual(self, vec: vector) -> vector:
        return vec / self.display_scale

    def _velocity_axis(self, velocity: vector) -> vector:
        speed = mag(velocity)
        if speed == 0:
            return vector(0, 0, 0)
        length_m = 0.12 * AU
        direction = norm(velocity)
        return self._to_visual(direction * length_m)

    def _setup_orbital_markers(self, earth_radius_vis: float) -> None:
        marker_radius = max(earth_radius_vis * 0.5, 0.02)
        self.perihelion_position = vector(EARTH_SEMI_MAJOR_AXIS * (1 - EARTH_ECCENTRICITY), 0, 0)
        self.aphelion_position = vector(-EARTH_SEMI_MAJOR_AXIS * (1 + EARTH_ECCENTRICITY), 0, 0)
        self.peri_marker = sphere(
            pos=self._to_visual(self.perihelion_position),
            radius=marker_radius,
            color=color.red,
            emissive=True,
        )
        self.aphelion_marker = sphere(
            pos=self._to_visual(self.aphelion_position),
            radius=marker_radius,
            color=color.blue,
            emissive=True,
        )
        self.peri_line = self._create_dashed_line(vector(0, 0, 0), self.perihelion_position, color.red)
        self.aphelion_line = self._create_dashed_line(vector(0, 0, 0), self.aphelion_position, color.blue)

    def _create_dashed_line(self, start: vector, end: vector, line_color: vector) -> List[curve]:
        segments = []
        total_segments = 18
        span = end - start
        for i in range(total_segments):
            frac_start = i / total_segments
            frac_end = min((i + 0.5) / total_segments, 1.0)
            p1 = start + span * frac_start
            p2 = start + span * frac_end
            segment = curve(color=line_color, radius=0.002)
            segment.append(pos=self._to_visual(p1))
            segment.append(pos=self._to_visual(p2))
            segments.append(segment)
        return segments

    # --- Physics core ------------------------------------------------------
    def _gravity(self, position: vector) -> vector:
        r = mag(position)
        return -MU * position / r**3

    def _log_state(self, pos: vector, vel: vector, accel: vector) -> None:
        r = mag(pos)
        speed = mag(vel)
        true_anomaly = math.atan2(pos.y, pos.x)
        specific_energy = 0.5 * speed**2 - MU / r
        area_rate = 0.5 * mag(cross(pos, vel))
        self.records.append(
            {
                "time_s": self.time_s,
                "x_m": pos.x,
                "y_m": pos.y,
                "z_m": pos.z,
                "vx_m_s": vel.x,
                "vy_m_s": vel.y,
                "vz_m_s": vel.z,
                "distance_m": r,
                "speed_m_s": speed,
                "true_anomaly_rad": true_anomaly,
                "specific_energy_J_kg": specific_energy,
                "area_rate_m2_s": area_rate,
                "acc_mag_m_s2": mag(accel),
            }
        )

        if self.visual:
            day = self.time_s / 86400.0
            self.distance_curve.plot(day, r / AU)
            self.energy_curve.plot(day, specific_energy / 1e6)
            self.orbit_trail.append(pos=self._to_visual(pos))
            self.earth.pos = self._to_visual(pos)
            if self.radial_line:
                self.radial_line.modify(0, vector(0, 0, 0))
                self.radial_line.modify(1, self._to_visual(pos))
            if self.velocity_arrow:
                self.velocity_arrow.pos = self._to_visual(pos)
                self.velocity_arrow.axis = self._velocity_axis(vel)
            self.stats_label.text = (
                f"Day {day:7.2f}\n"
                f"r = {r / AU:0.6f} AU\n"
                f"|v| = {speed / 1000:0.3f} km/s\n"
                f"E/m = {specific_energy / 1e6:0.3f} MJ/kg\n"
                f"θ = {math.degrees(true_anomaly):0.2f}°"
            )

    def step(self) -> None:
        dt = self.config.dt
        position = self.position
        velocity = self.velocity

        acceleration = self._gravity(position)
        position_next = position + velocity * dt + 0.5 * acceleration * dt**2
        next_acceleration = self._gravity(position_next)
        velocity_next = velocity + 0.5 * (acceleration + next_acceleration) * dt

        self.position = position_next
        self.velocity = velocity_next
        self.time_s += dt

        self._log_state(position_next, velocity_next, next_acceleration)

    def run(self) -> List[Dict[str, float]]:
        total_steps = int(self.config.total_time // self.config.dt)
        for _ in range(total_steps):
            if self.visual:
                rate(240)
            self.step()
        return self.records


# --- Analysis utilities -----------------------------------------------------
def _extract(records: List[Dict[str, float]], key: str) -> List[float]:
    return [row[key] for row in records]


def _unwrap_angles(values: List[float]) -> List[float]:
    if not values:
        return []
    unwrapped = [values[0]]
    for angle in values[1:]:
        value = angle
        prev = unwrapped[-1]
        delta = value - prev
        while delta <= -math.pi:
            value += 2 * math.pi
            delta = value - prev
        while delta > math.pi:
            value -= 2 * math.pi
            delta = value - prev
        unwrapped.append(value)
    return unwrapped


def compute_kepler_first(records: List[Dict[str, float]]) -> Dict[str, float]:
    theta = _unwrap_angles(_extract(records, "true_anomaly_rad"))
    theta_mod = [(angle % (2 * math.pi)) for angle in theta]
    radius_sim = _extract(records, "distance_m")
    radius_analytic = [
        EARTH_SEMI_MAJOR_AXIS * (1 - EARTH_ECCENTRICITY**2) / (1 + EARTH_ECCENTRICITY * math.cos(angle))
        for angle in theta_mod
    ]
    errors = [abs(rs - ra) / ra for rs, ra in zip(radius_sim, radius_analytic) if ra != 0]
    if radius_sim:
        rms = math.sqrt(sum((rs - ra) ** 2 for rs, ra in zip(radius_sim, radius_analytic)) / len(radius_sim))
    else:
        rms = float("nan")
    return {
        "max_rel_error_pct": max(errors) * 100 if errors else float("nan"),
        "r_rms_error_km": rms / 1000.0,
    }


def compute_kepler_second(records: List[Dict[str, float]], config: SimulationConfig) -> Dict[str, float]:
    dt = config.dt
    steps_per_segment = max(1, int(round(config.segment_duration / dt)))
    total_segments = len(records) // steps_per_segment
    if total_segments == 0:
        return {"segment_hours": config.segment_hours, "std_pct": float("nan"), "mean_area_Gm2": float("nan")}
    area_rate = _extract(records, "area_rate_m2_s")
    areas = []
    for seg_idx in range(total_segments):
        start = seg_idx * steps_per_segment
        end = start + steps_per_segment
        segment_area = sum(area_rate[start:end]) * dt
        areas.append(segment_area)
    mean_area = statistics.fmean(areas)
    std_area = statistics.pstdev(areas) if len(areas) > 1 else 0.0
    rel_std = (std_area / mean_area) * 100 if mean_area != 0 else float("nan")
    return {
        "segment_hours": config.segment_hours,
        "mean_area_Gm2": mean_area / 1e18,
        "std_pct": rel_std,
    }


def compute_kepler_third(records: List[Dict[str, float]]) -> Dict[str, float]:
    theta = _unwrap_angles(_extract(records, "true_anomaly_rad"))
    times = _extract(records, "time_s")
    if not theta:
        return {"period_days": float("nan"), "ratio_pct_err": float("nan")}
    target = 2 * math.pi
    for idx, angle in enumerate(theta):
        if angle >= target:
            if idx == 0:
                break
            theta_before, theta_after = theta[idx - 1], angle
            time_before, time_after = times[idx - 1], times[idx]
            fraction = (target - theta_before) / (theta_after - theta_before)
            period = time_before + fraction * (time_after - time_before)
            ratio_sim = (period**2) / (EARTH_SEMI_MAJOR_AXIS**3)
            ratio_theory = 4 * math.pi**2 / MU
            pct_err = (ratio_sim - ratio_theory) / ratio_theory * 100
            return {"period_days": period / 86400.0, "ratio_pct_err": pct_err}
    return {"period_days": float("nan"), "ratio_pct_err": float("nan")}


def analyze_records(records: List[Dict[str, float]], config: SimulationConfig) -> Dict[str, float]:
    k1 = compute_kepler_first(records)
    k2 = compute_kepler_second(records, config)
    k3 = compute_kepler_third(records)
    energy = _extract(records, "specific_energy_J_kg")
    if energy:
        energy_drift = (max(energy) - min(energy)) / abs(statistics.fmean(energy))
    else:
        energy_drift = float("nan")
    summary = {
        "samples": len(records),
        "dt_s": config.dt,
        "total_days": config.total_days,
        "k1_max_rel_error_pct": k1["max_rel_error_pct"],
        "k1_r_rms_error_km": k1["r_rms_error_km"],
        "k2_segment_hours": k2["segment_hours"],
        "k2_area_std_pct": k2["std_pct"],
        "k2_mean_area_Gm2": k2["mean_area_Gm2"],
        "k3_period_days": k3["period_days"],
        "k3_ratio_pct_err": k3["ratio_pct_err"],
        "energy_drift_pct": energy_drift * 100,
    }
    return summary


def save_outputs(records: List[Dict[str, float]], summary: Dict[str, float], config: SimulationConfig) -> None:
    config.log_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(records[0].keys()) if records else []
    with config.log_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    summary_path = config.log_path.with_suffix(".summary.json")
    with summary_path.open("w", encoding="utf-8") as fp:
        json.dump(summary, fp, indent=2)

    print(f"\nSaved {len(records)} samples to {config.log_path}")
    print(f"Wrote summary metrics to {summary_path}")


def print_summary(summary: Dict[str, float]) -> None:
    print("\n=== Kepler Verification Summary ===")
    print(f"Samples: {summary['samples']:,}")
    print(f"dt: {summary['dt_s']:.0f} s over {summary['total_days']:.1f} days")
    print(
        f"Kepler I  : max radius deviation = {summary['k1_max_rel_error_pct']:.4f}% "
        f"(RMS {summary['k1_r_rms_error_km']:.2f} km)"
    )
    print(
        f"Kepler II : area std = {summary['k2_area_std_pct']:.4f}% "
        f"over {summary['k2_segment_hours']:.1f} h segments (mean area {summary['k2_mean_area_Gm2']:.3f} Gm^2)"
    )
    print(
        f"Kepler III: simulated period = {summary['k3_period_days']:.6f} days, "
        f"T^2/a^3 error = {summary['k3_ratio_pct_err']:.6f}%"
    )
    print(f"Energy drift across run: {summary['energy_drift_pct']:.4f}%")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="High-fidelity VPython Sun–Earth orbit simulation.")
    parser.add_argument("--dt-minutes", type=float, default=10.0, help="Integrator step in minutes (default: 10).")
    parser.add_argument("--days", type=float, default=370.0, help="Total simulated days (default: 370).")
    parser.add_argument(
        "--segment-hours",
        type=float,
        default=6.0,
        help="Segment duration for Kepler II area comparison (default: 6 h).",
    )
    parser.add_argument(
        "--log-path",
        type=Path,
        default=Path("data/earth_orbit_log.csv"),
        help="CSV path for logged samples (default: data/earth_orbit_log.csv).",
    )
    parser.add_argument(
        "--headless",
        action="store_true",
        help="Skip VPython graphics (useful for automated testing or headless environments).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = SimulationConfig(
        dt=args.dt_minutes * 60.0,
        total_days=args.days,
        log_path=args.log_path,
        segment_hours=args.segment_hours,
    )
    sim = SunEarthSimulation(config=config, visual=not args.headless)
    df = sim.run()
    summary = analyze_records(df, config)
    save_outputs(df, summary, config)
    print_summary(summary)


if __name__ == "__main__":
    main()

