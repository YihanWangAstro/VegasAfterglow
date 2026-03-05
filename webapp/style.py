"""CSS styling and BibTeX/clipboard helpers for the webapp."""

SIDEBAR_CSS = """<style>
    /* ── Color tokens ── */
    :root {
        --accent: #4ecdc4;
        --bg: #0e1117;
        --bg-raised: #161b22;
        --bg-input: #21262d;
        --text: #e6edf3;
        --text-muted: #8b949e;
        --border: #30363d;
    }

    /* ── Main area: white background, dark text ── */
    [data-testid="stAppViewContainer"] {
        background-color: #fff !important;
    }
    [data-testid="stMainBlockContainer"] * {
        color: #333 !important;
    }
    [data-testid="stHeader"] {
        background-color: #fff !important;
    }
    [data-testid="stHeader"] * {
        color: #333 !important;
    }

    /* ── Global sidebar: compact, small font ── */
    [data-testid="stSidebar"] *:not(.katex):not(.katex *) {
        font-size: 0.72rem !important;
    }
    [data-testid="stSidebar"] .katex {
        font-size: 0.72rem !important;
    }
    [data-testid="stSidebar"] * {
        min-height: 0 !important;
    }
    [data-testid="stSidebar"] [data-testid="stVerticalBlock"] { gap: 0.4rem; }

    /* ── Widgets: remove extra padding ── */
    [data-testid="stSidebar"] .stSlider,
    [data-testid="stSidebar"] .stCheckbox,
    [data-testid="stSidebar"] .stSelectbox,
    [data-testid="stSidebar"] .stTextInput,
    [data-testid="stSidebar"] .stNumberInput {
        padding-top: 0; padding-bottom: 0;
    }

    /* ── Slider tick labels ── */
    [data-testid="stSidebar"] .stSlider [data-testid="stTickBarMin"],
    [data-testid="stSidebar"] .stSlider [data-testid="stTickBarMax"] {
        font-size: 0.6rem !important;
        color: var(--text-muted) !important;
    }

    /* ── Select box: compact + dark ── */
    [data-testid="stSidebar"] .stSelectbox [data-baseweb="select"] {
        background-color: transparent !important;
    }
    [data-testid="stSidebar"] .stSelectbox [data-baseweb="select"],
    [data-testid="stSidebar"] .stSelectbox [data-baseweb="select"] * {
        max-height: 1.5rem !important;
        padding-top: 0 !important;
        padding-bottom: 0 !important;
        line-height: 1.5rem !important;
    }
    [data-testid="stSidebar"] .stSelectbox [data-baseweb="select"] > div {
        background-color: var(--bg-input) !important;
        border-color: var(--border) !important;
    }
    [data-testid="stSidebar"] .stSelectbox [data-baseweb="select"] svg {
        width: 0.6rem !important; height: 0.6rem !important;
        max-height: none !important;
    }

    /* ── Dropdown menu (portal): compact + dark ── */
    [data-baseweb="popover"],
    [data-baseweb="menu"],
    [role="listbox"] {
        background-color: var(--bg-input) !important;
    }
    [data-baseweb="popover"] *,
    [data-baseweb="menu"] *,
    [role="listbox"] * {
        font-size: 0.72rem !important;
        min-height: 0 !important;
    }
    [data-baseweb="popover"] li,
    [role="listbox"] li {
        padding-top: 0.12rem !important;
        padding-bottom: 0.12rem !important;
    }
    [data-baseweb="popover"] li:hover,
    [role="listbox"] li:hover {
        background-color: var(--bg-raised) !important;
    }

    /* ── Text input: compact + dark ── */
    [data-testid="stSidebar"] .stTextInput [data-baseweb="input"],
    [data-testid="stSidebar"] .stTextInput [data-baseweb="base-input"] {
        background-color: var(--bg-input) !important;
        border-color: var(--border) !important;
        height: 1.5rem !important;
        min-height: 0 !important;
        padding: 0 !important;
    }
    [data-testid="stSidebar"] .stTextInput input {
        padding: 0 0.4rem !important;
        height: 1.5rem !important;
        line-height: 1.5rem !important;
        background-color: transparent !important;
        color: var(--text) !important;
    }

    /* ── Number input: compact + dark, side-by-side +/- ── */
    [data-testid="stSidebar"] [data-testid="stNumberInputContainer"] {
        display: flex !important;
        flex-direction: row !important;
        height: 1.5rem !important;
    }
    [data-testid="stSidebar"] .stNumberInput [data-baseweb="input"],
    [data-testid="stSidebar"] .stNumberInput [data-baseweb="base-input"] {
        background-color: var(--bg-input) !important;
        border-color: var(--border) !important;
        height: 1.5rem !important;
        min-height: 0 !important;
        padding: 0 !important;
    }
    [data-testid="stSidebar"] .stNumberInput input {
        padding: 0 0.4rem !important;
        height: 1.5rem !important;
        line-height: 1.5rem !important;
        background-color: transparent !important;
        color: var(--text) !important;
    }
    /* Button container: side by side instead of stacked */
    [data-testid="stSidebar"] [data-testid="stNumberInputContainer"] > div:last-child {
        display: flex !important;
        flex-direction: row !important;
        height: 1.5rem !important;
    }
    [data-testid="stSidebar"] .stNumberInput button {
        background-color: var(--bg-input) !important;
        border-color: var(--border) !important;
        color: var(--text) !important;
        height: 1.5rem !important;
        width: 1.2rem !important;
        min-width: 0 !important;
        padding: 0 !important;
    }
    [data-testid="stSidebar"] .stNumberInput button:hover {
        background-color: var(--bg-raised) !important;
    }

    /* ── Expander ── */
    [data-testid="stSidebar"] .stExpander {
        margin-top: 0.15rem; margin-bottom: 0.15rem;
    }
    [data-testid="stSidebar"] .stExpander details summary {
        padding: 0.25rem 0.4rem !important;
    }

    /* ── Sidebar grouped containers ── */
    [data-testid="stSidebar"] [data-testid="stVerticalBlockBorderWrapper"] {
        border: 1px solid var(--border) !important;
        border-radius: 0.5rem !important;
        background-color: rgba(255,255,255,0.03) !important;
        padding: 0.1rem !important;
        margin-bottom: 0.3rem !important;
    }

    /* ── Help tooltips: compact font ── */
    [data-baseweb="tooltip"] * {
        font-size: 0.72rem !important;
    }

    /* ── Action buttons: accent outline ── */
    .stDownloadButton button,
    [data-testid="stMainBlockContainer"] .stButton button {
        background: transparent !important;
        border: 1px solid var(--accent) !important;
        color: var(--accent) !important;
        font-size: 0.8rem !important;
        padding: 0.3rem 1rem !important;
        border-radius: 0.4rem !important;
        transition: background 0.2s, color 0.2s;
    }
    .stDownloadButton button:hover,
    [data-testid="stMainBlockContainer"] .stButton button:hover {
        background: var(--accent) !important;
        color: var(--bg) !important;
    }

    /* ── Mobile: responsive layout ── */
    @media (max-width: 768px) {
        [data-testid="stMainBlockContainer"] {
            padding-left: 0.5rem !important;
            padding-right: 0.5rem !important;
        }
        .stPlotlyChart, .js-plotly-plot {
            max-width: 100vw !important;
            overflow-x: auto !important;
        }
        .stDownloadButton button,
        [data-testid="stMainBlockContainer"] .stButton button {
            font-size: 0.7rem !important;
            padding: 0.25rem 0.5rem !important;
        }
    }
    </style>"""

BIBTEX = r"""@ARTICLE{2026JHEAp..5000490W,
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
}"""

_bibtex_js = BIBTEX.strip().replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")
CLIPBOARD_JS = f"""<script>
(function() {{
  var text = `{_bibtex_js}`;
  try {{
    window.parent.navigator.clipboard.writeText(text);
  }} catch(e) {{
    var ta = document.createElement('textarea');
    ta.value = text;
    ta.style.position = 'fixed';
    ta.style.left = '-9999px';
    document.body.appendChild(ta);
    ta.select();
    document.execCommand('copy');
    document.body.removeChild(ta);
  }}
}})();
</script>"""
