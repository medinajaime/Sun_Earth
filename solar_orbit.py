"""
Multi-body Solar System orbital simulation with VPython visualization and Kepler-law
diagnostics.

Physical constants sourced from:
- CODATA 2018 gravitational constant via NIST SP 959 (G = 6.67430e-11 m^3 kg^-1 s^-2).
- IAU 2015 Resolution B3 nominal solar/planetary parameters (masses, radii).
- NASA / JPL Planetary Fact Sheets (orbital elements referenced to the J2000 ecliptic).
- IAU 2012 Resolution B2: Astronomical Unit = 149,597,870,700 m.

Run `python solar_orbit.py --bodies earth` for the original Sun–Earth experience or pass
additional body identifiers (e.g., `--bodies mercury earth mars`) to add more planets.
Use `--list-bodies` to see the catalog. Include `--headless` to skip graphics while
still generating per-body logs and metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from vpython import (
    arrow,
    button,
    canvas,
    checkbox,
    color,
    curve,
    graph,
    gcurve,
    label,
    rate,
    slider,
    sphere,
    vector,
    wtext,
    cross,
    dot,
    mag,
)


# --- Fundamental constants (SI units) ---
G = 6.67430e-11  # CODATA 2018, m^3 kg^-1 s^-2
M_SUN = 1.98847e30  # kg, IAU 2015 nominal solar mass
M_EARTH = 5.972168e24  # kg, NASA Earth Fact Sheet
SUN_RADIUS = 6.957e8  # m, IAU 2015 nominal solar radius
EARTH_RADIUS = 6.371e6  # m, NASA
AU = 1.495978707e11  # m, IAU 2012 Resolution B2


@dataclass
class OrbitElements:
    """Keplerian elements referenced to the J2000 ecliptic."""

    semi_major_axis_m: float
    eccentricity: float
    inclination_deg: float
    ascending_node_deg: float
    arg_periapsis_deg: float
    mean_longitude_deg: Optional[float] = None
    mean_anomaly_deg: Optional[float] = None

    def __post_init__(self) -> None:
        if self.mean_anomaly_deg is None:
            if self.mean_longitude_deg is None:
                raise ValueError(
                    "Either mean_anomaly_deg or mean_longitude_deg must be provided."
                )
            peri_long = (self.ascending_node_deg + self.arg_periapsis_deg) % 360.0
            self.mean_anomaly_deg = (self.mean_longitude_deg - peri_long) % 360.0


@dataclass
class BodyDefinition:
    body_id: str
    name: str
    parent_id: Optional[str]
    mass_kg: float
    radius_m: float
    color_rgb: Tuple[float, float, float]
    texture_url: Optional[str]
    orbit: Optional[OrbitElements] = None
    trail_color_rgb: Optional[Tuple[float, float, float]] = None
    axial_tilt_deg: float = 0.0
    rotation_period_hours: Optional[float] = None


@dataclass
class BodyState:
    definition: BodyDefinition
    parent: Optional["BodyState"]
    mu: Optional[float]
    position: vector
    velocity: vector
    records: List[Dict[str, float]] = field(default_factory=list)
    sphere: Optional[sphere] = None
    trail: Optional[curve] = None
    label: Optional[label] = None
    spin_axis: vector = field(default_factory=lambda: vector(0, 0, 1))
    spin_rate_rad_s: float = 0.0
    axis_arrow: Optional[arrow] = None
    energy_reference: Optional[float] = None
    angular_reference: Optional[float] = None

    @property
    def absolute_position(self) -> vector:
        return self.position

    @property
    def absolute_velocity(self) -> vector:
        return self.velocity

    def relative_state(self) -> Tuple[vector, vector]:
        if self.parent:
            return self.position - self.parent.position, self.velocity - self.parent.velocity
        return self.position, self.velocity


BODY_DEFINITIONS: Dict[str, BodyDefinition] = {
    "sun": BodyDefinition(
        body_id="sun",
        name="Sun",
        parent_id=None,
        mass_kg=M_SUN,
        radius_m=SUN_RADIUS,
        color_rgb=(1.0, 0.65, 0.0),
        texture_url="https://i.postimg.cc/WjBQw4SW-/2k-sun.jpg",
        trail_color_rgb=(1.0, 0.5, 0.0),
        axial_tilt_deg=7.25,
        rotation_period_hours=609.12,
    ),
    "mercury": BodyDefinition(
        body_id="mercury",
        name="Mercury",
        parent_id="sun",
        mass_kg=3.3011e23,
        radius_m=2.4397e6,
        color_rgb=(0.7, 0.7, 0.7),
        texture_url="https://i.postimg.cc/x1trzQhw/2k-mercury.jpg",
        trail_color_rgb=(0.8, 0.8, 0.8),
        axial_tilt_deg=0.03,
        rotation_period_hours=1407.5,
        orbit=OrbitElements(
            semi_major_axis_m=0.38709893 * AU,
            eccentricity=0.20563069,
            inclination_deg=7.00487,
            ascending_node_deg=48.33167,
            arg_periapsis_deg=29.12478,
            mean_longitude_deg=252.25084,
        ),
    ),
    "venus": BodyDefinition(
        body_id="venus",
        name="Venus",
        parent_id="sun",
        mass_kg=4.8675e24,
        radius_m=6.0518e6,
        color_rgb=(1.0, 0.8, 0.6),
        texture_url="https://i.postimg.cc/3xnMvhzh/2k-venus-surface.jpg",
        trail_color_rgb=(1.0, 0.85, 0.65),
        axial_tilt_deg=177.36,
        rotation_period_hours=-5832.5,
        orbit=OrbitElements(
            semi_major_axis_m=0.72333199 * AU,
            eccentricity=0.00677323,
            inclination_deg=3.39471,
            ascending_node_deg=76.68069,
            arg_periapsis_deg=54.85229,
            mean_longitude_deg=181.97973,
        ),
    ),
    "earth": BodyDefinition(
        body_id="earth",
        name="Earth",
        parent_id="sun",
        mass_kg=M_EARTH,
        radius_m=EARTH_RADIUS,
        color_rgb=(0.3, 0.6, 1.0),
        texture_url="https://i.postimg.cc/ncRX2TcM/2k-earth-daymap.jpg",
        trail_color_rgb=(0.2, 0.8, 1.0),
        axial_tilt_deg=23.44,
        rotation_period_hours=23.934,
        orbit=OrbitElements(
            semi_major_axis_m=1.00000011 * AU,
            eccentricity=0.01671022,
            inclination_deg=0.00005,
            ascending_node_deg=-11.26064,
            arg_periapsis_deg=114.20783,
            mean_longitude_deg=100.46435,
        ),
    ),
    "moon": BodyDefinition(
        body_id="moon",
        name="Moon",
        parent_id="earth",
        mass_kg=7.342e22,
        radius_m=1.7374e6,
        color_rgb=(0.8, 0.8, 0.75),
        texture_url=None,
        trail_color_rgb=(0.9, 0.9, 0.8),
        axial_tilt_deg=6.68,
        rotation_period_hours=655.728,
        orbit=OrbitElements(
            semi_major_axis_m=3.844e8,
            eccentricity=0.0549,
            inclination_deg=5.145,
            ascending_node_deg=125.08,
            arg_periapsis_deg=318.15,
            mean_longitude_deg=218.3166,
        ),
    ),
    "mars": BodyDefinition(
        body_id="mars",
        name="Mars",
        parent_id="sun",
        mass_kg=6.4171e23,
        radius_m=3.3895e6,
        color_rgb=(1.0, 0.4, 0.2),
        texture_url="https://i.postimg.cc/Cxmyq0XX/2k-mars.jpg",
        trail_color_rgb=(1.0, 0.5, 0.3),
        axial_tilt_deg=25.19,
        rotation_period_hours=24.623,
        orbit=OrbitElements(
            semi_major_axis_m=1.52366231 * AU,
            eccentricity=0.09341233,
            inclination_deg=1.85061,
            ascending_node_deg=49.57854,
            arg_periapsis_deg=286.4623,
            mean_longitude_deg=355.45332,
        ),
    ),
    "jupiter": BodyDefinition(
        body_id="jupiter",
        name="Jupiter",
        parent_id="sun",
        mass_kg=1.8982e27,
        radius_m=6.9911e7,
        color_rgb=(1.0, 0.8, 0.5),
        texture_url="https://i.postimg.cc/zG0mgN6c/2k-jupiter.jpg",
        trail_color_rgb=(1.0, 0.85, 0.55),
        axial_tilt_deg=3.13,
        rotation_period_hours=9.925,
        orbit=OrbitElements(
            semi_major_axis_m=5.20336301 * AU,
            eccentricity=0.04839266,
            inclination_deg=1.3053,
            ascending_node_deg=100.55615,
            arg_periapsis_deg=274.1977,
            mean_longitude_deg=34.40438,
        ),
    ),
    "saturn": BodyDefinition(
        body_id="saturn",
        name="Saturn",
        parent_id="sun",
        mass_kg=5.6834e26,
        radius_m=5.8232e7,
        color_rgb=(0.9, 0.8, 0.6),
        texture_url="https://i.postimg.cc/T3C8DxB8/2k-saturn.jpg",
        trail_color_rgb=(0.95, 0.85, 0.65),
        axial_tilt_deg=26.73,
        rotation_period_hours=10.7,
        orbit=OrbitElements(
            semi_major_axis_m=9.53707032 * AU,
            eccentricity=0.0541506,
            inclination_deg=2.48446,
            ascending_node_deg=113.71504,
            arg_periapsis_deg=338.7169,
            mean_longitude_deg=49.94432,
        ),
    ),
    "uranus": BodyDefinition(
        body_id="uranus",
        name="Uranus",
        parent_id="sun",
        mass_kg=8.6810e25,
        radius_m=2.5362e7,
        color_rgb=(0.6, 0.9, 0.9),
        texture_url="https://i.postimg.cc/vZ3RVyCF/2k-uranus.jpg",
        trail_color_rgb=(0.65, 0.95, 0.95),
        axial_tilt_deg=97.77,
        rotation_period_hours=-17.24,
        orbit=OrbitElements(
            semi_major_axis_m=19.19126393 * AU,
            eccentricity=0.04716771,
            inclination_deg=0.76986,
            ascending_node_deg=74.22988,
            arg_periapsis_deg=96.73436,
            mean_longitude_deg=313.23218,
        ),
    ),
    "neptune": BodyDefinition(
        body_id="neptune",
        name="Neptune",
        parent_id="sun",
        mass_kg=1.02413e26,
        radius_m=2.4622e7,
        color_rgb=(0.4, 0.6, 1.0),
        texture_url="https://i.postimg.cc/52ncCJTr/2k-neptune.jpg",
        trail_color_rgb=(0.5, 0.7, 1.0),
        axial_tilt_deg=28.32,
        rotation_period_hours=16.11,
        orbit=OrbitElements(
            semi_major_axis_m=30.06896348 * AU,
            eccentricity=0.00858587,
            inclination_deg=1.76917,
            ascending_node_deg=131.72169,
            arg_periapsis_deg=273.24966,
            mean_longitude_deg=304.88003,
        ),
    ),
}


@dataclass
class SimulationConfig:
    dt: float = 600.0  # seconds
    total_days: float = 370.0
    log_dir: Path = Path("data")
    segment_hours: float = 6.0
    body_ids: Sequence[str] = field(default_factory=lambda: ["earth"])
    dt_min: float = 60.0
    dt_max: float = 3600.0
    adaptive_dt: bool = False
    adaptive_distance_fraction: float = 0.02
    adaptive_step_scale: float = 1.4

    @property
    def total_time(self) -> float:
        return self.total_days * 86400.0

    @property
    def segment_duration(self) -> float:
        return self.segment_hours * 3600.0

    def __post_init__(self) -> None:
        self.dt_min = max(1e-6, self.dt_min)
        self.dt_max = max(self.dt_min, self.dt_max)
        self.adaptive_step_scale = max(1.0, self.adaptive_step_scale)


def solve_kepler(mean_anomaly_rad: float, eccentricity: float, tol: float = 1e-10, max_iter: int = 50) -> float:
    """Solve Kepler's equation for the eccentric anomaly."""

    if eccentricity < 0.8:
        E = mean_anomaly_rad
    else:
        E = math.pi
    for _ in range(max_iter):
        f = E - eccentricity * math.sin(E) - mean_anomaly_rad
        fp = 1 - eccentricity * math.cos(E)
        delta = f / fp
        E -= delta
        if abs(delta) < tol:
            break
    return E


def perifocal_to_inertial(vec: vector, elements: OrbitElements) -> vector:
    """Rotate a vector from the perifocal frame to the inertial J2000 ecliptic frame."""

    Omega = math.radians(elements.ascending_node_deg)
    inc = math.radians(elements.inclination_deg)
    argp = math.radians(elements.arg_periapsis_deg)

    cO = math.cos(Omega)
    sO = math.sin(Omega)
    ci = math.cos(inc)
    si = math.sin(inc)
    cw = math.cos(argp)
    sw = math.sin(argp)

    r11 = cO * cw - sO * sw * ci
    r12 = -cO * sw - sO * cw * ci
    r13 = sO * si
    r21 = sO * cw + cO * sw * ci
    r22 = -sO * sw + cO * cw * ci
    r23 = -cO * si
    r31 = sw * si
    r32 = cw * si
    r33 = ci

    x = r11 * vec.x + r12 * vec.y + r13 * vec.z
    y = r21 * vec.x + r22 * vec.y + r23 * vec.z
    z = r31 * vec.x + r32 * vec.y + r33 * vec.z
    return vector(x, y, z)


def elements_to_state(elements: OrbitElements, mu: float) -> Tuple[vector, vector]:
    """Return Cartesian state vectors (position, velocity) from orbital elements."""

    e = elements.eccentricity
    a = elements.semi_major_axis_m
    mean_anomaly_rad = math.radians(elements.mean_anomaly_deg or 0.0)

    eccentric_anomaly = solve_kepler(mean_anomaly_rad, e)
    cos_E = math.cos(eccentric_anomaly)
    sin_E = math.sin(eccentric_anomaly)
    r = a * (1 - e * cos_E)

    true_anomaly = 2 * math.atan2(math.sqrt(1 + e) * math.sin(eccentric_anomaly / 2), math.sqrt(1 - e) * math.cos(eccentric_anomaly / 2))
    position_pf = vector(r * math.cos(true_anomaly), r * math.sin(true_anomaly), 0)

    velocity_factor = math.sqrt(mu * a) / r
    velocity_pf = vector(-sin_E, math.sqrt(1 - e**2) * cos_E, 0) * velocity_factor

    return perifocal_to_inertial(position_pf, elements), perifocal_to_inertial(velocity_pf, elements)


def spin_axis_from_tilt(tilt_deg: float) -> vector:
    tilt_rad = math.radians(tilt_deg)
    axis = vector(math.sin(tilt_rad), 0.0, math.cos(tilt_rad))
    if mag(axis) == 0:
        return vector(0, 0, 1)
    return axis / mag(axis)


def compute_orbital_elements(rel_pos: vector, rel_vel: vector, mu: float) -> Dict[str, float]:
    r = mag(rel_pos)
    v = mag(rel_vel)
    if r == 0 or mu == 0:
        return {}
    h_vec = cross(rel_pos, rel_vel)
    h = mag(h_vec)
    k_vec = vector(0, 0, 1)
    n_vec = cross(k_vec, h_vec)
    n = mag(n_vec)
    e_vec = ((rel_pos * (v**2 - mu / r)) - (rel_vel * dot(rel_pos, rel_vel))) / mu
    e = mag(e_vec)
    energy = 0.5 * v**2 - mu / r
    semi_major_axis = float("inf") if energy == 0 else -mu / (2 * energy)

    def clamp(value: float) -> float:
        return max(-1.0, min(1.0, value))

    inclination = math.degrees(math.acos(clamp(h_vec.z / h))) if h != 0 else 0.0
    raan = 0.0
    if n != 0:
        raan = math.degrees(math.acos(clamp(n_vec.x / n)))
        if n_vec.y < 0:
            raan = (360.0 - raan) % 360.0
    arg_periapsis = 0.0
    if n != 0 and e > 1e-8:
        arg_periapsis = math.degrees(math.acos(clamp(dot(n_vec, e_vec) / (n * e))))
        if e_vec.z < 0:
            arg_periapsis = (360.0 - arg_periapsis) % 360.0
    period = float("nan")
    if semi_major_axis > 0 and math.isfinite(semi_major_axis):
        period = 2 * math.pi * math.sqrt(semi_major_axis**3 / mu)
    return {
        "semi_major_axis_m": semi_major_axis,
        "eccentricity": e,
        "inclination_deg": inclination,
        "raan_deg": raan,
        "arg_periapsis_deg": arg_periapsis,
        "period_s": period,
    }


class SolarSystemSimulation:
    """VPython-based velocity-Verlet integrator for configurable Solar-System subsets."""

    def __init__(self, config: SimulationConfig, visual: bool = True) -> None:
        self.config = config
        self.visual = visual
        self.time_s = 0.0
        self.current_dt = self.config.dt
        self.display_scale = AU
        self.scene = None
        self.info_label = None
        self.distance_graph = None
        self.energy_graph = None
        self.angular_graph = None
        self.distance_curves: Dict[str, gcurve] = {}
        self.energy_curves: Dict[str, gcurve] = {}
        self.angular_curves: Dict[str, gcurve] = {}
        self.paused = False
        self.time_scale_factor = 1.0
        self.trails_enabled = True
        self.play_button = None
        self.speed_slider = None
        self.speed_label = None
        self.trail_checkbox = None
        self.auto_checkbox = None
        self.camera_targets: List[BodyState] = []
        self.camera_target_index = 0
        self.camera_mode = "free"
        self.auto_rescale = True
        self.zoom_factor = 1.0
        self.manual_range = 1.0

        self.bodies: List[BodyState] = []
        self.body_map: Dict[str, BodyState] = {}
        self._build_body_states()
        self.camera_targets = self._camera_target_candidates()
        self._log_all_states()  # capture epoch data

        if self.visual:
            self._setup_scene()
            self._setup_controls()

    # --- Initialization helpers -------------------------------------------
    def _build_body_states(self) -> None:
        requested = set()
        for body_id in self.config.body_ids:
            if body_id not in BODY_DEFINITIONS:
                raise ValueError(f"Unknown body '{body_id}'. Use --list-bodies to inspect catalog.")
            cursor = body_id
            while cursor:
                if cursor in requested:
                    break
                requested.add(cursor)
                cursor = BODY_DEFINITIONS[cursor].parent_id

        ordered: List[str] = []
        visited: set[str] = set()

        def visit(body_id: str) -> None:
            if body_id in visited:
                return
            definition = BODY_DEFINITIONS[body_id]
            if definition.parent_id:
                visit(definition.parent_id)
            ordered.append(body_id)
            visited.add(body_id)

        for body_id in requested:
            visit(body_id)

        for body_id in ordered:
            definition = BODY_DEFINITIONS[body_id]
            parent_state = self.body_map.get(definition.parent_id)
            mu = None
            position_rel = vector(0, 0, 0)
            velocity_rel = vector(0, 0, 0)
            if definition.parent_id and definition.orbit:
                parent_mass = parent_state.definition.mass_kg if parent_state else M_SUN
                mu = G * (definition.mass_kg + parent_mass)
                position_rel, velocity_rel = elements_to_state(definition.orbit, mu)
            abs_position = position_rel
            abs_velocity = velocity_rel
            if parent_state:
                abs_position = parent_state.position + position_rel
                abs_velocity = parent_state.velocity + velocity_rel
            spin_axis = spin_axis_from_tilt(definition.axial_tilt_deg)
            spin_rate = 0.0
            if definition.rotation_period_hours:
                period_s = definition.rotation_period_hours * 3600.0
                if period_s != 0:
                    spin_rate = (2 * math.pi) / period_s
            body_state = BodyState(
                definition=definition,
                parent=parent_state,
                mu=mu,
                position=abs_position,
                velocity=abs_velocity,
                spin_axis=spin_axis,
                spin_rate_rad_s=spin_rate,
            )
            self.body_map[body_id] = body_state
            self.bodies.append(body_state)

    def _camera_target_candidates(self) -> List[BodyState]:
        targets = [self.body_map[body_id] for body_id in self.config.body_ids if body_id in self.body_map]
        return targets or list(self.bodies)

    # --- Scene helpers -----------------------------------------------------
    def _setup_scene(self) -> None:
        self.scene = canvas(
            title="Solar-System Orbits (Velocity-Verlet, real units)",
            width=1200,
            height=750,
            background=color.black,
        )

        max_a = max(
            (state.definition.orbit.semi_major_axis_m for state in self.bodies if state.definition.orbit),
            default=AU,
        )
        self.scene.range = max(1.5, 1.2 * (max_a / self.display_scale))
        self.scene.forward = vector(-1.2, -0.6, -1.0)
        self.scene.up = vector(0, 0, 1)

        self.distance_graph = graph(
            title="Distance vs. Time",
            xtitle="Days",
            ytitle="Distance (AU)",
            fast=False,
            align="left",
            width=600,
            height=280,
        )
        self.energy_graph = graph(
            title="Energy Drift",
            xtitle="Days",
            ytitle="ΔEnergy (%)",
            fast=False,
            align="right",
            width=600,
            height=280,
        )
        self.angular_graph = graph(
            title="Angular Momentum Drift",
            xtitle="Days",
            ytitle="Δ|h| (%)",
            fast=False,
            align="left",
            width=600,
            height=280,
        )

        for state in self.bodies:
            color_vec = vector(*state.definition.color_rgb)
            trail_color_vec = (
                vector(*state.definition.trail_color_rgb)
                if state.definition.trail_color_rgb
                else color_vec
            )
            scale_factor = 40 if state.definition.body_id == "sun" else 400
            min_radius = 0.15 if state.definition.body_id == "sun" else 0.04
            radius_vis = max(
                min_radius,
                state.definition.radius_m / self.display_scale * scale_factor,
            )
            state.sphere = sphere(
                pos=self._to_visual(state.absolute_position),
                radius=radius_vis,
                color=color_vec,
                texture=state.definition.texture_url,
                emissive=True if state.definition.body_id == "sun" else False,
                shininess=0.2,
                make_trail=True,
            )
            axis_dir = state.spin_axis if mag(state.spin_axis) else vector(0, 0, 1)
            axis_dir = axis_dir / mag(axis_dir)
            axis_vec = axis_dir * (state.sphere.radius * 2.2)
            state.sphere.axis = axis_vec
            state.trail = curve(color=trail_color_vec, radius=0.01, retain=2000)
            state.label = label(
                pos=self._to_visual(state.absolute_position),
                text=state.definition.name,
                color=color_vec,
                box=False,
                height=12,
            )
            state.axis_arrow = arrow(
                pos=self._to_visual(state.absolute_position),
                axis=axis_vec,
                color=color.white,
                shaftwidth=state.sphere.radius * 0.25,
                opacity=0.4,
            )
            if self.distance_graph and state.definition.orbit:
                self.distance_curves[state.definition.body_id] = gcurve(
                    graph=self.distance_graph,
                    color=trail_color_vec,
                    label=state.definition.name,
                )
            if self.energy_graph and state.definition.orbit:
                self.energy_curves[state.definition.body_id] = gcurve(
                    graph=self.energy_graph,
                    color=color_vec,
                )
            if self.angular_graph and state.definition.orbit:
                self.angular_curves[state.definition.body_id] = gcurve(
                    graph=self.angular_graph,
                    color=trail_color_vec,
                )

        self.info_label = label(
            pos=vector(-self.scene.range * 0.9, self.scene.range * 1.05, 0),
            text="Initializing...",
            color=color.white,
            box=False,
            height=16,
        )
        self.manual_range = self.scene.range
        self.scene.bind("keydown", self._handle_keypress)

    def _setup_controls(self) -> None:
        if not self.visual or not self.scene:
            return
        self.scene.append_to_caption("\n\nControls:\n")
        self.play_button = button(text="Pause", bind=self._toggle_pause_button)
        self.scene.append_to_caption("  Speed: ")
        self.speed_slider = slider(
            min=0.1,
            max=5.0,
            value=self.time_scale_factor,
            step=0.1,
            length=200,
            bind=self._set_speed_from_slider,
        )
        self.speed_label = wtext(text=f"{self.time_scale_factor:0.1f}x")
        self.scene.append_to_caption("\n")
        self.trail_checkbox = checkbox(text="Show trails", checked=self.trails_enabled, bind=self._set_trails_checked)
        self.scene.append_to_caption("   ")
        self.auto_checkbox = checkbox(text="Auto-rescale", checked=self.auto_rescale, bind=self._set_auto_rescale_checked)
        self.scene.append_to_caption("\n(Keyboard: space=pause, f=follow, c/n=cycle, +/-=zoom, a=auto, t=trails)\n")

    def _handle_keypress(self, evt) -> None:
        if not self.visual:
            return
        raw_key = getattr(evt, "key", "")
        key = raw_key.lower()
        if "+" in key:
            key = key.split("+")[-1]
        if key == "f":
            self._toggle_camera_mode()
        elif key in {"c", "n"}:
            self._cycle_camera_target(1)
        elif key in {"b", "p"}:
            self._cycle_camera_target(-1)
        elif key in {"=", "+"}:
            self._adjust_zoom(0.9)
        elif key in {"-", "_"}:
            self._adjust_zoom(1.1)
        elif key == "a":
            self._toggle_auto_rescale()
        elif key in {" ", "space"}:
            self._set_paused(not self.paused)
        elif key == "t":
            self._set_trails_enabled(not self.trails_enabled)

    def _set_paused(self, value: bool) -> None:
        self.paused = value
        if self.play_button:
            self.play_button.text = "Resume" if self.paused else "Pause"

    def _toggle_pause_button(self, _ctrl) -> None:
        self._set_paused(not self.paused)

    def _set_speed_from_slider(self, slider_ctrl) -> None:
        self.time_scale_factor = slider_ctrl.value
        if self.speed_label:
            self.speed_label.text = f"{self.time_scale_factor:0.1f}x"

    def _set_trails_enabled(self, enabled: bool) -> None:
        self.trails_enabled = enabled
        for state in self.bodies:
            if state.trail:
                state.trail.visible = self.trails_enabled
                if not self.trails_enabled:
                    state.trail.clear()
        if self.trail_checkbox:
            self.trail_checkbox.checked = self.trails_enabled

    def _set_trails_checked(self, checkbox_ctrl) -> None:
        self._set_trails_enabled(bool(checkbox_ctrl.checked))

    def _set_auto_rescale(self, enabled: bool) -> None:
        self.auto_rescale = enabled
        if self.auto_rescale:
            self.zoom_factor = 1.0
        else:
            if self.scene:
                self.manual_range = self.scene.range
        if self.auto_checkbox:
            self.auto_checkbox.checked = self.auto_rescale

    def _set_auto_rescale_checked(self, checkbox_ctrl) -> None:
        self._set_auto_rescale(bool(checkbox_ctrl.checked))

    def _toggle_camera_mode(self) -> None:
        self.camera_mode = "follow" if self.camera_mode == "free" else "free"
        if self.camera_mode == "follow":
            self._update_camera(force=True)

    def _cycle_camera_target(self, direction: int) -> None:
        if not self.camera_targets:
            return
        self.camera_target_index = (self.camera_target_index + direction) % len(self.camera_targets)
        if self.camera_mode == "follow":
            self._update_camera(force=True)

    def _toggle_auto_rescale(self) -> None:
        self._set_auto_rescale(not self.auto_rescale)

    def _adjust_zoom(self, factor: float) -> None:
        if self.auto_rescale:
            self.zoom_factor = max(0.2, min(5.0, self.zoom_factor * factor))
        else:
            self.manual_range = max(0.1, min(100.0, self.manual_range * factor))
            if self.scene:
                self.scene.range = self.manual_range
        if self.camera_mode == "follow":
            self._update_camera(force=self.auto_rescale is False)

    def _current_camera_target(self) -> Optional[BodyState]:
        if not self.camera_targets:
            return None
        self.camera_target_index %= len(self.camera_targets)
        return self.camera_targets[self.camera_target_index]

    def _to_visual(self, vec_world: vector) -> vector:
        return vec_world / self.display_scale

    # --- Simulation core ---------------------------------------------------
    def _log_all_states(self) -> None:
        for state in self.bodies:
            if state.definition.orbit and state.mu:
                rel_pos, rel_vel = state.relative_state()
                self._log_state(state, rel_pos, rel_vel, 0.0)

    def _next_dt(self, remaining: float = float("inf")) -> float:
        if not self.config.adaptive_dt:
            dt = min(self.config.dt, remaining)
            self.current_dt = dt
            return dt
        base_dt = self.current_dt or self.config.dt
        proposals: List[float] = []
        for state in self.bodies:
            if not state.definition.orbit or not state.mu:
                continue
            rel_pos, rel_vel = state.relative_state()
            distance = mag(rel_pos)
            speed = mag(rel_vel)
            if distance <= 0 or speed <= 0:
                continue
            proposals.append(self.config.adaptive_distance_fraction * distance / speed)
        dt_candidate = min(proposals) if proposals else base_dt
        dt_candidate = max(self.config.dt_min, min(self.config.dt_max, dt_candidate))
        max_scale = self.config.adaptive_step_scale
        dt_candidate = min(dt_candidate, base_dt * max_scale)
        dt_candidate = max(dt_candidate, base_dt / max_scale)
        dt = min(dt_candidate, remaining)
        self.current_dt = dt
        return dt

    def _update_camera(self, force: bool = False) -> None:
        if not self.visual or not self.scene:
            return
        target = self._current_camera_target()
        if self.camera_mode != "follow" or target is None:
            return
        center = self._to_visual(target.absolute_position)
        self.scene.center = center
        if self.auto_rescale:
            rel_pos, _ = target.relative_state()
            distance = mag(rel_pos)
            if distance <= 0:
                distance = target.definition.radius_m * 5
            base_range = max(0.15, 1.5 * (distance / self.display_scale))
            desired_range = base_range * self.zoom_factor
        else:
            desired_range = self.manual_range
        if force:
            self.scene.range = desired_range
        else:
            smoothing = 0.2
            self.scene.range = (1 - smoothing) * self.scene.range + smoothing * desired_range

    def _compute_accelerations(self, positions: List[vector]) -> List[vector]:
        accelerations: List[vector] = []
        for idx, pos in enumerate(positions):
            accel = vector(0, 0, 0)
            for other_idx, other_pos in enumerate(positions):
                if other_idx == idx:
                    continue
                rel = other_pos - pos
                dist_sq = rel.x**2 + rel.y**2 + rel.z**2
                if dist_sq == 0:
                    continue
                inv_r3 = 1.0 / (dist_sq * math.sqrt(dist_sq))
                accel += rel * (G * self.bodies[other_idx].definition.mass_kg * inv_r3)
            accelerations.append(accel)
        return accelerations

    def _log_state(self, state: BodyState, rel_pos: vector, rel_vel: vector, step_dt: float) -> None:
        abs_pos = state.absolute_position
        abs_vel = state.absolute_velocity
        r = mag(rel_pos)
        speed = mag(rel_vel)
        elements: Dict[str, float] = {}
        energy_drift_pct = 0.0
        angular_drift_pct = 0.0
        if state.mu:
            h_vec = cross(rel_pos, rel_vel)
            h_mag = mag(h_vec)
            area_rate = 0.5 * h_mag
            specific_energy = 0.5 * speed**2 - state.mu / r
            true_anomaly = self._true_anomaly(rel_pos, rel_vel, state.mu, state.definition.orbit.eccentricity, state.definition.orbit.semi_major_axis_m)
            elements = compute_orbital_elements(rel_pos, rel_vel, state.mu)
            if state.energy_reference is None:
                state.energy_reference = specific_energy
            if state.angular_reference is None:
                state.angular_reference = h_mag
            if state.energy_reference not in (None, 0):
                energy_drift_pct = ((specific_energy - state.energy_reference) / abs(state.energy_reference)) * 100.0
            if state.angular_reference not in (None, 0):
                angular_drift_pct = ((h_mag - state.angular_reference) / state.angular_reference) * 100.0
        else:
            specific_energy = 0.0
            area_rate = 0.0
            true_anomaly = 0.0

        state.records.append(
            {
                "time_s": self.time_s,
                "dt_s": step_dt,
                "x_m": abs_pos.x,
                "y_m": abs_pos.y,
                "z_m": abs_pos.z,
                "vx_m_s": abs_vel.x,
                "vy_m_s": abs_vel.y,
                "vz_m_s": abs_vel.z,
                "distance_m": r,
                "speed_m_s": speed,
                "true_anomaly_rad": true_anomaly,
                "specific_energy_J_kg": specific_energy,
                "area_rate_m2_s": area_rate,
                "energy_drift_pct": energy_drift_pct,
                "angular_momentum_drift_pct": angular_drift_pct,
                "semi_major_axis_m_live": elements.get("semi_major_axis_m", float("nan")),
                "eccentricity_live": elements.get("eccentricity", float("nan")),
                "inclination_deg_live": elements.get("inclination_deg", float("nan")),
                "raan_deg_live": elements.get("raan_deg", float("nan")),
                "arg_periapsis_deg_live": elements.get("arg_periapsis_deg", float("nan")),
                "period_days_live": (elements.get("period_s", float("nan")) / 86400.0)
                if math.isfinite(elements.get("period_s", float("nan")))
                else float("nan"),
            }
        )

        if self.visual:
            day = self.time_s / 86400.0
            if state.definition.body_id in self.distance_curves:
                self.distance_curves[state.definition.body_id].plot(day, r / AU)
            if state.definition.body_id in self.energy_curves:
                self.energy_curves[state.definition.body_id].plot(day, energy_drift_pct)
            if state.definition.body_id in self.angular_curves:
                self.angular_curves[state.definition.body_id].plot(day, angular_drift_pct)
            if self.trails_enabled and state.trail:
                state.trail.append(pos=self._to_visual(abs_pos))
            if state.sphere:
                state.sphere.pos = self._to_visual(abs_pos)
                if state.spin_rate_rad_s and step_dt > 0:
                    rotation_angle = state.spin_rate_rad_s * step_dt
                    if rotation_angle:
                        state.sphere.rotate(
                            angle=rotation_angle,
                            axis=state.spin_axis,
                            origin=state.sphere.pos,
                        )
            if state.label:
                state.label.pos = self._to_visual(abs_pos)
            if state.axis_arrow and state.sphere:
                axis_dir = state.spin_axis if mag(state.spin_axis) else vector(0, 0, 1)
                axis_dir = axis_dir / mag(axis_dir)
                axis_vec = axis_dir * (state.sphere.radius * 2.2)
                state.axis_arrow.pos = self._to_visual(abs_pos)
                state.axis_arrow.axis = axis_vec

        if self.visual and self.info_label:
            primary = self.config.body_ids[0] if self.config.body_ids else "earth"
            if primary in self.body_map and self.body_map[primary].records:
                latest = self.body_map[primary].records[-1]
                target = self._current_camera_target()
                camera_line = f"Cam: {self.camera_mode}"
                if target:
                    camera_line += f" → {target.definition.name}"
                camera_line += " (space:pause f:follow c/n:cycle +/-:zoom a:auto t:trails)"
                e_live = latest.get("eccentricity_live", float("nan"))
                a_live = latest.get("semi_major_axis_m_live", float("nan"))
                incl_live = latest.get("inclination_deg_live", float("nan"))
                raan_live = latest.get("raan_deg_live", float("nan"))
                argp_live = latest.get("arg_periapsis_deg_live", float("nan"))
                period_live = latest.get("period_days_live", float("nan"))
                self.info_label.text = (
                    f"t = {self.time_s/86400.0:8.2f} days\n"
                    f"{self.body_map[primary].definition.name}: {latest['distance_m']/AU:0.3f} AU\n"
                    f"|v| = {latest['speed_m_s']/1000:0.3f} km/s\n"
                    f"e={e_live:.4f}  a={a_live/AU if math.isfinite(a_live) else float('nan'):.3f} AU\n"
                    f"i={incl_live:.2f}°  Ω={raan_live:.2f}°  ω={argp_live:.2f}°\n"
                    f"P={period_live:.3f} d\n"
                    f"{camera_line}"
                )

    def _true_anomaly(
        self,
        rel_pos: vector,
        rel_vel: vector,
        mu: float,
        eccentricity: float,
        semi_major_axis: float,
    ) -> float:
        r = mag(rel_pos)
        if eccentricity == 0:
            return math.atan2(rel_pos.y, rel_pos.x)
        cos_nu = ((semi_major_axis * (1 - eccentricity**2) / r) - 1) / eccentricity
        cos_nu = max(-1.0, min(1.0, cos_nu))
        vr = dot(rel_pos, rel_vel) / r
        h_mag = mag(cross(rel_pos, rel_vel))
        if mu == 0 or eccentricity == 0 or h_mag == 0:
            sin_nu = 0.0
        else:
            sin_nu = (vr * h_mag) / (mu * eccentricity)
        sin_nu = max(-1.0, min(1.0, sin_nu))
        return math.atan2(sin_nu, cos_nu)
    def step(self, dt_override: Optional[float] = None) -> None:
        if dt_override is not None:
            dt = dt_override
            self.current_dt = dt
        else:
            dt = self._next_dt()
        if dt <= 0.0:
            return
        current_positions = [state.position for state in self.bodies]
        accelerations = self._compute_accelerations(current_positions)
        new_positions: List[vector] = []
        for state, accel in zip(self.bodies, accelerations):
            new_positions.append(
                    state.position
                    + state.velocity * dt
                    + accel * (0.5 * dt**2)
)

        next_accelerations = self._compute_accelerations(new_positions)
        for idx, state in enumerate(self.bodies):
            accel = accelerations[idx]
            next_accel = next_accelerations[idx]
            state.position = new_positions[idx]
            state.velocity = state.velocity + (accel + next_accel) * (0.5 * dt)

        self.time_s += dt
        for state in self.bodies:
            if state.definition.orbit and state.mu:
                rel_pos, rel_vel = state.relative_state()
                self._log_state(state, rel_pos, rel_vel, dt)
        if self.visual:
            self._update_camera()

    def run(self) -> Dict[str, List[Dict[str, float]]]:
        while self.time_s < self.config.total_time:
            if self.visual:
                rate(240)
                if self.paused:
                    continue
            remaining = self.config.total_time - self.time_s
            base_dt = self._next_dt(remaining)
            dt = max(self.config.dt_min, min(self.config.dt_max, base_dt * self.time_scale_factor))
            dt = min(dt, remaining)
            if dt <= 0:
                break
            self.step(dt_override=dt)
        return {state.definition.body_id: state.records for state in self.bodies if state.definition.orbit}


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


def compute_kepler_first(
    records: List[Dict[str, float]], semi_major_axis: float, eccentricity: float
) -> Dict[str, float]:
    theta = _unwrap_angles(_extract(records, "true_anomaly_rad"))
    theta_mod = [(angle % (2 * math.pi)) for angle in theta]
    radius_sim = _extract(records, "distance_m")
    radius_analytic = [
        semi_major_axis * (1 - eccentricity**2) / (1 + eccentricity * math.cos(angle))
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


def compute_kepler_second(
    records: List[Dict[str, float]], config: SimulationConfig
) -> Dict[str, float]:
    if not records:
        return {"segment_hours": config.segment_hours, "std_pct": float("nan"), "mean_area_Gm2": float("nan")}
    segment_seconds = config.segment_duration
    if segment_seconds <= 0:
        return {"segment_hours": config.segment_hours, "std_pct": float("nan"), "mean_area_Gm2": float("nan")}
    areas: List[float] = []
    current_time = records[0]["time_s"]
    next_boundary = current_time + segment_seconds
    accum_area = 0.0
    for idx in range(1, len(records)):
        rate = records[idx]["area_rate_m2_s"]
        prev_time = records[idx - 1]["time_s"]
        end_time = records[idx]["time_s"]
        cursor = max(current_time, prev_time)
        step_duration = max(0.0, end_time - prev_time)
        remaining = step_duration
        while remaining > 0:
            time_to_boundary = next_boundary - cursor
            if time_to_boundary <= 0:
                next_boundary += segment_seconds
                continue
            slice_dt = min(remaining, time_to_boundary)
            accum_area += rate * slice_dt
            cursor += slice_dt
            remaining -= slice_dt
            if cursor >= next_boundary - 1e-9:
                areas.append(accum_area)
                accum_area = 0.0
                next_boundary += segment_seconds
        current_time = cursor
    if accum_area > 0:
        areas.append(accum_area)
    if not areas:
        return {"segment_hours": config.segment_hours, "std_pct": float("nan"), "mean_area_Gm2": float("nan")}
    mean_area = statistics.fmean(areas)
    std_area = statistics.pstdev(areas) if len(areas) > 1 else 0.0
    rel_std = (std_area / mean_area) * 100 if mean_area != 0 else float("nan")
    return {
        "segment_hours": config.segment_hours,
        "mean_area_Gm2": mean_area / 1e18,
        "std_pct": rel_std,
    }


def compute_kepler_third(
    records: List[Dict[str, float]], semi_major_axis: float, mu: float
) -> Dict[str, float]:
    theta = _unwrap_angles(_extract(records, "true_anomaly_rad"))
    times = _extract(records, "time_s")
    if not theta:
        return {"period_days": float("nan"), "ratio_pct_err": float("nan")}
    start_angle = theta[0]
    target = start_angle + 2 * math.pi
    for idx, angle in enumerate(theta):
        if angle >= target:
            if idx == 0:
                break
            theta_before, theta_after = theta[idx - 1], angle
            time_before, time_after = times[idx - 1], times[idx]
            fraction = (target - theta_before) / (theta_after - theta_before)
            period = time_before + fraction * (time_after - time_before)
            ratio_sim = (period**2) / (semi_major_axis**3)
            ratio_theory = 4 * math.pi**2 / mu
            pct_err = (ratio_sim - ratio_theory) / ratio_theory * 100
            return {"period_days": period / 86400.0, "ratio_pct_err": pct_err}
    return {"period_days": float("nan"), "ratio_pct_err": float("nan")}


def analyze_records(
    records_by_body: Dict[str, List[Dict[str, float]]],
    config: SimulationConfig,
    bodies: Sequence[BodyState],
) -> Dict[str, Dict[str, float]]:
    summaries: Dict[str, Dict[str, float]] = {}
    body_lookup = {state.definition.body_id: state for state in bodies}
    for body_id, records in records_by_body.items():
        state = body_lookup[body_id]
        orbit = state.definition.orbit
        if not orbit or state.mu is None:
            continue
        k1 = compute_kepler_first(records, orbit.semi_major_axis_m, orbit.eccentricity)
        k2 = compute_kepler_second(records, config)
        k3 = compute_kepler_third(records, orbit.semi_major_axis_m, state.mu)
        energy = _extract(records, "specific_energy_J_kg")
        dt_samples = [row.get("dt_s", config.dt) for row in records if row.get("dt_s", 0.0) > 0]
        mean_dt = statistics.fmean(dt_samples) if dt_samples else config.dt
        if energy:
            energy_drift = (max(energy) - min(energy)) / abs(statistics.fmean(energy))
        else:
            energy_drift = float("nan")
        summaries[body_id] = {
            "name": state.definition.name,
            "samples": len(records),
            "dt_s": mean_dt,
            "dt_mode": "adaptive" if config.adaptive_dt else "fixed",
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
    return summaries


def save_outputs(
    records_by_body: Dict[str, List[Dict[str, float]]],
    summary_by_body: Dict[str, Dict[str, float]],
    config: SimulationConfig,
) -> None:
    config.log_dir.mkdir(parents=True, exist_ok=True)
    combined_summary: Dict[str, Dict[str, float]] = {}
    for body_id, records in records_by_body.items():
        log_path = config.log_dir / f"{body_id}_orbit_log.csv"
        fieldnames = list(records[0].keys()) if records else []
        with log_path.open("w", newline="", encoding="utf-8") as fp:
            writer = csv.DictWriter(fp, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        summary = summary_by_body.get(body_id, {})
        summary_path = log_path.with_suffix(".summary.json")
        with summary_path.open("w", encoding="utf-8") as fp:
            json.dump(summary, fp, indent=2)
        combined_summary[body_id] = summary
        print(f"Saved {len(records)} samples to {log_path}")
        print(f"Wrote summary metrics to {summary_path}")

    combined_path = config.log_dir / "summary.json"
    with combined_path.open("w", encoding="utf-8") as fp:
        json.dump(combined_summary, fp, indent=2)
    print(f"Aggregated summary written to {combined_path}")


def print_summary(summary_by_body: Dict[str, Dict[str, float]]) -> None:
    print("\n=== Kepler Verification Summary ===")
    if not summary_by_body:
        print("No orbiting bodies were analyzed.")
        return
    for body_id, summary in summary_by_body.items():
        print(f"\n-- {summary.get('name', body_id.title())} --")
        print(f"Samples: {summary['samples']:,}")
        mode = summary.get("dt_mode", "fixed")
        print(f"dt ({mode}): {summary['dt_s']:.0f} s over {summary['total_days']:.1f} days")
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
    parser = argparse.ArgumentParser(description="High-fidelity VPython Solar-System orbit simulation.")
    parser.add_argument("--dt-minutes", type=float, default=10.0, help="Integrator step in minutes (default: 10).")
    parser.add_argument(
        "--adaptive-dt",
        action="store_true",
        help="Enable adaptive timestep controller to refine periapsis resolution.",
    )
    parser.add_argument(
        "--dt-min-minutes",
        type=float,
        default=1.0,
        help="Minimum timestep (minutes) when adaptive dt is enabled (default: 1).",
    )
    parser.add_argument(
        "--dt-max-minutes",
        type=float,
        default=60.0,
        help="Maximum timestep (minutes) when adaptive dt is enabled (default: 60).",
    )
    parser.add_argument(
        "--distance-fraction",
        type=float,
        default=0.02,
        help="Max fraction of primary distance a body may travel per step (default: 0.02).",
    )
    parser.add_argument(
        "--step-scale",
        type=float,
        default=1.4,
        help="Limit multiplier applied to timestep growth/shrink per step (default: 1.4).",
    )
    parser.add_argument("--days", type=float, default=370.0, help="Total simulated days (default: 370).")
    parser.add_argument(
        "--segment-hours",
        type=float,
        default=6.0,
        help="Segment duration for Kepler II area comparison (default: 6 h).",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path("data"),
        help="Directory for per-body CSV/JSON logs (default: data/).",
    )
    parser.add_argument(
        "--bodies",
        nargs="+",
        default=["mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune"],
        help="Body identifiers to simulate (default: earth). Parents are added automatically.",
    )
    parser.add_argument(
        "--list-bodies",
        action="store_true",
        help="List the available body identifiers and exit.",
    )
    parser.add_argument(
        "--headless",
        action="store_true",
        help="Skip VPython graphics (useful for automated testing or headless environments).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.list_bodies:
        print("Available bodies (parent in parentheses):")
        for body_id, definition in BODY_DEFINITIONS.items():
            parent = definition.parent_id or "-"
            orbit_note = "static" if definition.orbit is None else f"orbits {parent}"
            print(f"- {body_id:<9} : {definition.name} ({orbit_note})")
        return

    requested_bodies = [body.lower() for body in args.bodies]
    config = SimulationConfig(
        dt=args.dt_minutes * 60.0,
        total_days=args.days,
        log_dir=args.log_dir,
        segment_hours=args.segment_hours,
        body_ids=requested_bodies,
        dt_min=args.dt_min_minutes * 60.0,
        dt_max=args.dt_max_minutes * 60.0,
        adaptive_dt=args.adaptive_dt,
        adaptive_distance_fraction=args.distance_fraction,
        adaptive_step_scale=args.step_scale,
    )
    sim = SolarSystemSimulation(config=config, visual=not args.headless)
    records = sim.run()
    summary = analyze_records(records, config, sim.bodies)
    save_outputs(records, summary, config)
    print_summary(summary)


if __name__ == "__main__":
    main()

