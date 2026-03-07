from __future__ import annotations

from pydantic import BaseModel, Field


class ObservationGroup(BaseModel):
    legend: str = "data"
    x_unit: str | None = None
    y_unit: str = "mJy"
    text: str = ""
    visible: bool = True


class SharedParams(BaseModel):
    d_L_mpc: float = 100.0
    theta_obs: float = 0.0
    flux_unit: str = "mJy"
    time_unit: str = "s"
    jet_type: str = "Top-hat"
    theta_c: float = 0.1
    E_iso: float = 1e52
    Gamma0: float = 300.0
    spreading: bool = False
    duration: float = 1.0
    k_e: float = 2.0
    k_g: float = 2.0
    theta_w: float = 0.3
    E_iso_w: float = 1e51
    Gamma0_w: float = 100.0
    medium_type: str = "ISM"
    n_ism: float = 1.0
    A_star: float = 0.1
    k_m: float = 2.0
    eps_e: float = 0.1
    eps_B: float = 1e-3
    p: float = 2.3
    xi_e: float = 1.0
    ssc: bool = False
    kn: bool = False
    enable_rvs: bool = False
    eps_e_r: float = 0.1
    eps_B_r: float = 1e-3
    p_r: float = 2.3
    xi_e_r: float = 1.0
    rvs_ssc: bool = False
    rvs_kn: bool = False
    num_t: int = 100
    res_phi: float = 0.1
    res_theta: float = 0.25
    res_t: float = 10.0


class LightCurveRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    frequencies_input: str = "1e9, R, 1keV"
    t_min: float = 1.0
    t_max: float = 1e8
    selected_instruments: list[str] = Field(default_factory=list)
    observation_groups: list[ObservationGroup] = Field(default_factory=list)
    include_figure: bool = True
    include_exports: bool = True
    export_kinds: list[str] = Field(default_factory=list)


class SpectrumRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    t_snapshots_input: str = "1e3, 1e4, 1e5, 1e6"
    nu_min: float = 1e8
    nu_max: float = 1e20
    num_nu: int = 200
    freq_unit: str = "Hz"
    show_nufnu: bool = False
    selected_instruments: list[str] = Field(default_factory=list)
    observation_groups: list[ObservationGroup] = Field(default_factory=list)
    include_figure: bool = True
    include_exports: bool = True
    export_kinds: list[str] = Field(default_factory=list)


class SkyMapRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    animate: bool = False
    t_obs: float = 1e6
    t_min: float = 1e4
    t_max: float = 1e7
    n_frames: int = 15
    nu_input: str = "1e9"
    fov: float = 500.0
    npixel: int = 256
    include_figure: bool = True
    include_exports: bool = True
    export_kinds: list[str] = Field(default_factory=list)
