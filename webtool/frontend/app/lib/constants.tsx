import type { ReactNode } from "react";
import type { AxisName, Mode, SharedLinearSliderSpec, SharedLogSliderSpec, SharedParams, SharedSliderSpec } from "./types";

export const PLOT_CONFIG = { responsive: true, displaylogo: false } as const;
export const ALL_AXES: AxisName[] = ["xaxis", "yaxis", "yaxis2"];

export const FALLBACK_SHARED: SharedParams = {
  d_L_mpc: 100,
  z: 0.022,
  theta_obs: 0,
  flux_unit: "mJy",
  time_unit: "s",
  jet_type: "Top-hat",
  theta_c: 0.1,
  E_iso: 1e52,
  Gamma0: 300,
  spreading: false,
  duration: 1.0,
  k_e: 2.0,
  k_g: 2.0,
  theta_w: 0.3,
  E_iso_w: 1e51,
  Gamma0_w: 100,
  medium_type: "ISM",
  n_ism: 1.0,
  A_star: 0.1,
  k_m: 2.0,
  eps_e: 0.1,
  eps_B: 1e-3,
  p: 2.3,
  xi_e: 1.0,
  ssc: false,
  kn: false,
  enable_rvs: false,
  eps_e_r: 0.1,
  eps_B_r: 1e-3,
  p_r: 2.3,
  xi_e_r: 1.0,
  rvs_ssc: false,
  rvs_kn: false,
  num_t: 100,
  res_phi: 0.1,
  res_theta: 0.25,
  res_t: 10,
};

export const Y_UNIT_OPTIONS = ["mJy", "Jy", "uJy", "cgs", "AB mag", "erg/cm²/s"] as const;
export const PLOT_FLUX_UNIT_OPTIONS = ["mJy", "Jy", "uJy", "cgs", "AB mag"] as const;
export const TIME_UNIT_OPTIONS = ["s", "day", "hr", "min"] as const;
export const FREQ_UNIT_OPTIONS = ["Hz", "GHz", "keV", "MeV"] as const;

export const SKY_DEFAULT_N_FRAMES = 15;
export const SKY_MAX_N_FRAMES = 30;
export const SKY_MAX_PIXEL_ANIMATE = 512;
export const SKY_MAX_PIXEL_STATIC = 1024;
export const SKY_PIXEL_OPTIONS = [64, 128, 256, 512, 1024] as const;

export const ENABLE_INTERACTIVE_DOWNSAMPLE = true;
export const INTERACTIVE_LIGHTCURVE_NUM_T_MAX = 120;
export const INTERACTIVE_SPECTRUM_NUM_NU_MAX = 120;
export const INTERACTIVE_SKY_PIXEL_MAX_STATIC = 256;
export const INTERACTIVE_SKY_PIXEL_MAX_ANIMATE = 256;
export const INTERACTIVE_SKY_FRAMES_MAX = 8;

export const AUTO_RUN_DEBOUNCE_IDLE_MS = 10;
export const SLIDER_COMMIT_INTERVAL_MS = 20;
export const SPECTRUM_TEXT_COMMIT_DEBOUNCE_MS = 220;
export const API_STATUS_REFRESH_DIRECT_MS = 2000;
export const COLD_START_HINT_MS = 1200;

export const AXIS_EPS = 1e-6;
export const H0_KM_S_MPC = 67.4;
export const C_KM_S = 299792.458;

export const MODE_OPTIONS: { value: Mode; label: string }[] = [
  { value: "lightcurve", label: "Light Curve" },
  { value: "spectrum", label: "Spectrum" },
  { value: "skymap", label: "Sky Image" },
];

export const INSTRUMENT_TYPE_BY_NAME: Record<string, string> = {
  VLA: "Radio",
  ALMA: "Radio",
  MeerKAT: "Radio",
  ngVLA: "Radio",
  "Rubin/LSST": "Optical/IR",
  JWST: "Optical/IR",
  WFST: "Optical/IR",
  "SVOM/VT": "Optical/IR",
  "Swift/XRT": "X-ray",
  Chandra: "X-ray",
  "EP/WXT": "X-ray",
  "EP/FXT": "X-ray",
  "SVOM/MXT": "X-ray",
  "SVOM/ECLAIRs": "X-ray",
  "Swift/BAT": "Gamma-ray",
  "Fermi/GBM": "Gamma-ray",
  "Fermi/LAT": "Gamma-ray",
  CTA: "Gamma-ray",
};

export const INSTRUMENT_GROUP_ORDER = ["Radio", "Optical/IR", "X-ray", "Gamma-ray", "Other"] as const;

export const OBS_STORAGE_LC_KEY = "afterglow:webtool:obs:lightcurve";
export const OBS_STORAGE_SED_KEY = "afterglow:webtool:obs:spectrum";
export const BOOKMARKS_STORAGE_KEY = "afterglow:webtool:bookmarks:v1";
export const MAX_BOOKMARKS = 30;
export const URL_STATE_PARAM = "state";
export const URL_STATE_VERSION = 1;
export const URL_STATE_MAX_CHARS = 3500;

export const FREQ_HELP_TEXT = [
  "Accepts Hz values (1e9), unit suffixes (1GHz, 1keV), filter names, instrument bands, or custom ranges ([0.3keV,10keV]).",
  "",
  "Filters: U B V R I J H Ks | u b v uvw1 uvm2 uvw2 | g r i z | w VT_B VT_R |",
  "F225W F275W F336W F438W F475W F555W F606W F625W F775W F814W F850LP F105W F110W F125W F140W F160W",
  "",
  "Bands: XRT BAT FXT WXT MXT ECLAIRs LAT GBM",
  "",
  "Units: Hz kHz MHz GHz | eV keV MeV GeV",
].join("\n");

export const JET_TYPE_OPTIONS: SharedParams["jet_type"][] = ["Top-hat", "Gaussian", "Power-law", "Two-component"];
export const MEDIUM_TYPE_OPTIONS: { value: SharedParams["medium_type"]; label: string }[] = [
  { value: "ISM", label: "ISM" },
  { value: "Wind", label: "Wind" },
  { value: "Wind bubble", label: "Wind+ISM" },
];

export const DEFAULT_APP_VERSION = "2.0.1";
export const DOWNLOAD_TEXT_META = {
  csv: { ext: "csv", mime: "text/csv" },
  json: { ext: "json", mime: "application/json" },
} as const;

export const BIBTEX_TEXT = `@ARTICLE{2026JHEAp..5000490W,
       author = {{Wang}, Yihan and {Chen}, Connery and {Zhang}, Bing},
        title = "{VegasAfterglow: A high-performance framework for gamma-ray burst afterglows}",
      journal = {Journal of High Energy Astrophysics},
         year = 2026,
        month = feb,
       volume = {50},
          eid = {100490},
        pages = {100490},
          doi = {10.1016/j.jheap.2025.100490},
archivePrefix = {arXiv},
       eprint = {2507.10829},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2026JHEAp..5000490W},
}

@book{Zhang2018,
  author    = {Zhang, Bing},
  title     = {{The Physics of Gamma-Ray Bursts}},
  publisher = {Cambridge University Press},
  year      = {2018},
  doi       = {10.1017/9781139226530}
}`;

const logSpec = (
  key: keyof SharedParams,
  label: ReactNode,
  minExp: number,
  maxExp: number,
  step: number,
  defaultExp: number,
): SharedLogSliderSpec => ({ kind: "log", key, label, minExp, maxExp, step, defaultExp });

const linearSpec = (
  key: keyof SharedParams,
  label: ReactNode,
  min: number,
  max: number,
  step: number,
  decimals?: number,
): SharedLinearSliderSpec => ({ kind: "linear", key, label, min, max, step, decimals });

export const RADIATION_SLIDER_ROWS: SharedSliderSpec[][] = [
  [
    logSpec(
      "eps_e",
      <>
        log10(ε<sub>e</sub>)
      </>,
      -5,
      0,
      0.05,
      -1,
    ),
    logSpec(
      "eps_B",
      <>
        log10(ε<sub>B</sub>)
      </>,
      -6,
      0,
      0.05,
      -3,
    ),
  ],
  [
    linearSpec("p", "p", 2.01, 3.0, 0.01),
    logSpec(
      "xi_e",
      <>
        log10(ξ<sub>e</sub>)
      </>,
      -3,
      0,
      0.05,
      0,
    ),
  ],
];

export const REVERSE_SHOCK_DURATION_SPEC = logSpec("duration", "log10(Jet Duration (s))", -1, 4, 0.05, 0);

export const REVERSE_SHOCK_SLIDER_ROWS: SharedSliderSpec[][] = [
  [
    logSpec(
      "eps_e_r",
      <>
        log10(ε<sub>e,r</sub>)
      </>,
      -5,
      0,
      0.05,
      -1,
    ),
    logSpec(
      "eps_B_r",
      <>
        log10(ε<sub>B,r</sub>)
      </>,
      -6,
      0,
      0.05,
      -3,
    ),
  ],
  [
    linearSpec(
      "p_r",
      <>
        p<sub>r</sub>
      </>,
      2.01,
      3,
      0.01,
    ),
    logSpec(
      "xi_e_r",
      <>
        log10(ξ<sub>e,r</sub>)
      </>,
      -3,
      0,
      0.05,
      0,
    ),
  ],
];

export const RESOLUTION_SLIDERS: SharedLinearSliderSpec[] = [
  linearSpec("res_phi", "φ ppd", 0.05, 0.5, 0.05, 2),
  linearSpec("res_theta", "θ ppd", 0.1, 2, 0.05, 2),
  linearSpec("res_t", "t ppd", 1, 30, 0.5, 1),
];

export const JET_THETA_C_SPEC = linearSpec(
  "theta_c",
  <>
    θ<sub>c</sub>
  </>,
  0.01,
  1.0,
  0.01,
);

export const JET_ENERGY_ROW: SharedSliderSpec[] = [
  logSpec(
    "E_iso",
    <>
      log10(E<sub>iso</sub> (erg))
    </>,
    48,
    57,
    0.1,
    52,
  ),
  logSpec(
    "Gamma0",
    <>
      log10(Γ<sub>0</sub>)
    </>,
    1,
    3.5,
    0.05,
    Math.log10(300),
  ),
];

export const JET_POWERLAW_ROW: SharedSliderSpec[] = [
  linearSpec(
    "k_e",
    <>
      k<sub>e</sub>
    </>,
    0.5,
    10,
    0.1,
  ),
  linearSpec(
    "k_g",
    <>
      k<sub>g</sub>
    </>,
    0.5,
    10,
    0.1,
  ),
];

export const JET_TWOCOMP_ROW: SharedSliderSpec[] = [
  linearSpec(
    "theta_w",
    <>
      θ<sub>w</sub>
    </>,
    0.05,
    1.5,
    0.01,
  ),
  logSpec(
    "E_iso_w",
    <>
      log10(E<sub>iso,w</sub>)
    </>,
    48,
    55,
    0.1,
    51,
  ),
];

export const JET_TWOCOMP_GAMMA_SPEC = logSpec(
  "Gamma0_w",
  <>
    log10(Γ<sub>0,w</sub>)
  </>,
  1,
  3,
  0.05,
  2,
);

export const MEDIUM_ISM_SPEC = logSpec(
  "n_ism",
  <>
    log10(n<sub>ism</sub> (cm⁻³))
  </>,
  -5,
  5,
  0.1,
  0,
);

export const MEDIUM_WIND_SPEC = logSpec(
  "A_star",
  <>
    log10(A<sub>*</sub>)
  </>,
  -3,
  2,
  0.1,
  -1,
);

export const MEDIUM_KM_SPEC = linearSpec(
  "k_m",
  <>
    k<sub>m</sub>
  </>,
  0,
  4,
  0.1,
);

export const MEDIUM_WIND_FLOOR_SPEC = logSpec(
  "n_ism",
  <>
    log10(n<sub>floor</sub> (cm⁻³))
  </>,
  -8,
  5,
  0.1,
  0,
);
