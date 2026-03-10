from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


class SharedParams(BaseModel):
    d_L_mpc: float = Field(default=100.0, gt=0, le=1e6)
    theta_obs: float = Field(default=0.0, ge=0, le=1.57)
    jet_type: Literal["Top-hat", "Gaussian", "Power-law", "Two-component"] = "Top-hat"
    theta_c: float = Field(default=0.1, gt=0, le=1.57)
    E_iso: float = Field(default=1e52, ge=1e40, le=1e58)
    Gamma0: float = Field(default=300.0, ge=1.0, le=1e5)
    spreading: bool = False
    duration: float = Field(default=1.0, ge=0, le=1e8)
    k_e: float = Field(default=2.0, ge=0, le=20)
    k_g: float = Field(default=2.0, ge=0, le=20)
    theta_w: float = Field(default=0.3, gt=0, le=1.57)
    E_iso_w: float = Field(default=1e51, ge=1e40, le=1e58)
    Gamma0_w: float = Field(default=100.0, ge=1.0, le=1e5)
    medium_type: Literal["ISM", "Wind bubble", "Wind"] = "ISM"
    n_ism: float = Field(default=1.0, ge=1e-10, le=1e10)
    A_star: float = Field(default=0.1, ge=1e-10, le=1e10)
    k_m: float = Field(default=2.0, ge=0, le=4)
    eps_e: float = Field(default=0.1, gt=0, le=1)
    eps_B: float = Field(default=1e-3, gt=0, le=1)
    p: float = Field(default=2.3, gt=2.0, le=5.0)
    xi_e: float = Field(default=1.0, gt=0, le=1)
    ssc: bool = False
    kn: bool = False
    enable_rvs: bool = False
    eps_e_r: float = Field(default=0.1, gt=0, le=1)
    eps_B_r: float = Field(default=1e-3, gt=0, le=1)
    p_r: float = Field(default=2.3, gt=2.0, le=5.0)
    xi_e_r: float = Field(default=1.0, gt=0, le=1)
    rvs_ssc: bool = False
    rvs_kn: bool = False
    num_t: int = Field(default=100, ge=2, le=500)
    res_phi: float = Field(default=0.1, gt=0, le=1)
    res_theta: float = Field(default=0.25, gt=0, le=1)
    res_t: float = Field(default=10.0, gt=0, le=100)


MAX_FREQUENCIES = 30
MAX_SNAPSHOTS = 20


class LightCurveRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    frequencies_input: str = Field(default="1e9, R, 1keV", max_length=500)
    t_min: float = Field(default=1.0, gt=0)
    t_max: float = Field(default=1e8, gt=0)


class SpectrumRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    t_snapshots_input: str = Field(default="1e3, 1e4, 1e5, 1e6", max_length=500)
    nu_min: float = Field(default=1e8, gt=0)
    nu_max: float = Field(default=1e20, gt=0)
    num_nu: int = Field(default=200, ge=2, le=1000)


class SkyMapRequest(BaseModel):
    shared: SharedParams = Field(default_factory=SharedParams)
    animate: bool = False
    t_obs: float = Field(default=1e6, gt=0)
    t_min: float = Field(default=1e4, gt=0)
    t_max: float = Field(default=1e7, gt=0)
    n_frames: int = Field(default=15, ge=2, le=30)
    nu_input: str = Field(default="1e9", max_length=200)
    fov: float = Field(default=500.0, gt=0)
    npixel: int = Field(default=256, ge=2, le=1024)
