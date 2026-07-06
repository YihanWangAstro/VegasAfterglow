"""Shared constants and helpers for the validation suite (no plotting deps)."""

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

try:
    import VegasAfterglow as va
    HAS_VA = True
except ImportError:
    va = None
    HAS_VA = False

ROOT_DIR = Path(__file__).parents[2]

# Convergence thresholds
MAX_ERROR_THRESHOLD, MEAN_ERROR_THRESHOLD = 0.15, 0.05


def masked_relative_errors(flux, ref, log_t):
    """Relative errors |flux/ref - 1| under the suite's convergence mask:
    positive finite flux in both curves, |dlogF/dlogt| <= 4 on the reference
    (steep de-beamed rises carry meaningless relative errors), and above 1e-6
    of the reference band peak. Returns None when nothing survives the mask."""
    flux = np.asarray(flux)
    ref = np.asarray(ref)
    valid = (ref > 0) & np.isfinite(ref) & (flux > 0) & np.isfinite(flux)
    if not valid.any():
        return None
    log_ref = np.full_like(log_t, np.nan)
    pos = ref > 0
    log_ref[pos] = np.log10(ref[pos])
    slope = np.gradient(log_ref, log_t)
    valid &= np.abs(slope) <= 4
    if not valid.any():
        return None
    valid &= ref > ref[valid].max() * 1e-6
    return np.abs(flux[valid] - ref[valid]) / ref[valid]

# Fiducial resolution settings (single source of truth, used by benchmark and report)
FIDUCIAL_RESOLUTION = (0.1, 0.25, 10)
DIM_INDEX = {"phi": 0, "theta": 1, "t": 2}
FIDUCIAL_VALUES = {dim: FIDUCIAL_RESOLUTION[idx] for dim, idx in DIM_INDEX.items()}


def find_fiducial_index(values, fiducial: float) -> int:
    if not values:
        return -1
    for i, v in enumerate(values):
        if abs(v - fiducial) < 1e-9:
            return i
    diffs = [abs(v - fiducial) for v in values]
    return diffs.index(min(diffs))


def worst_fiducial_error(errors_by_band: Dict, fid_idx: int) -> float:
    """Return the worst (max absolute) error across all bands at the fiducial index."""
    worst = 0.0
    for errs in errors_by_band.values():
        if errs and 0 <= fid_idx < len(errs):
            v = errs[fid_idx]
            if v is not None and v == v:
                worst = max(worst, abs(v))
    return worst


def get_git_commit() -> str:
    """Get the short git commit hash of HEAD."""
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def load_json_safe(path: Path) -> Tuple[Optional[Dict], Optional[str]]:
    """Load JSON with error handling. Returns (data, None) on success or (None, error_msg) on failure."""
    if not path.exists():
        return None, f"Results not found: {path}"
    try:
        return json.loads(path.read_text()), None
    except json.JSONDecodeError as e:
        return None, f"Invalid JSON: {e}"


def _parse_cmake_flags(cmake_path: Path) -> str:
    """Parse compiler flags from CMakeLists.txt based on current platform."""
    import platform
    import re

    if not cmake_path.exists():
        return "unknown"

    try:
        content = cmake_path.read_text()
        system = platform.system()

        if system == "Windows":
            # Look for MSVC flags
            match = re.search(r'elseif\s*\(\s*MSVC\s*\)\s*\n\s*add_compile_options\s*\(([^)]+)\)', content)
            if match:
                flags = match.group(1).strip()
                return flags.replace('\n', ' ').strip()
        else:
            # GNU/Clang flags (macOS/Linux)
            flags = []
            in_gnu_section = False
            for line in content.split('\n'):
                if 'CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"' in line:
                    in_gnu_section = True
                elif in_gnu_section:
                    if 'elseif' in line or 'endif' in line:
                        break
                    match = re.search(r'add_compile_options\s*\(([^)]+)\)', line)
                    if match:
                        flags.append(match.group(1).strip())

            if flags:
                return ' '.join(flags)
    except Exception:
        pass

    return "unknown"


def get_runtime_build_info() -> Dict[str, str]:
    """Get build/system information at runtime."""
    import platform
    info = {"Version": "unknown", "Python": sys.version.split()[0], "Platform": f"{platform.system()} {platform.machine()}",
            "Commit": "unknown", "Compiler": "unknown", "Flags": "unknown", "CPU": "unknown"}

    # Get version from VegasAfterglow package (trim commit suffix like +gb28e78a26.d20260204)
    if HAS_VA:
        full_version = getattr(va, "__version__", "unknown")
        # Remove the +commit.date suffix for cleaner display
        info["Version"] = full_version.split("+")[0] if "+" in full_version else full_version

    # Get CPU info
    try:
        cpu_brand = None
        # On macOS, use sysctl for detailed CPU info (platform.processor() only returns "arm" or "i386")
        if platform.system() == "Darwin":
            result = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                cpu_brand = result.stdout.strip()
        # On Linux, parse /proc/cpuinfo
        elif platform.system() == "Linux":
            try:
                with open("/proc/cpuinfo") as f:
                    for line in f:
                        if line.startswith("model name"):
                            cpu_brand = line.split(":")[1].strip()
                            break
            except (IOError, OSError):
                pass
        # Fallback to platform.processor()
        if not cpu_brand:
            cpu_brand = platform.processor()
        if cpu_brand and cpu_brand not in ("", "arm", "i386", "x86_64", "aarch64"):
            info["CPU"] = cpu_brand
        elif cpu_brand:
            # Use machine type as fallback with arch
            info["CPU"] = f"{platform.machine()} ({cpu_brand})"
    except Exception:
        pass

    info["Commit"] = get_git_commit()

    # Get compiler info from compile_commands.json (if available)
    compile_commands = ROOT_DIR / "compile_commands.json"
    if compile_commands.exists():
        try:
            data = json.loads(compile_commands.read_text())
            if data:
                parts = data[0].get("command", "").split()
                if parts:
                    compiler_name = Path(parts[0]).name
                    try:
                        result = subprocess.run([compiler_name, "--version"], capture_output=True, text=True, check=True)
                        info["Compiler"] = result.stdout.split("\n")[0].strip()
                    except (subprocess.CalledProcessError, FileNotFoundError):
                        info["Compiler"] = compiler_name
        except Exception:
            pass

    # Get compiler flags from CMakeLists.txt (more reliable than compile_commands.json)
    cmake_flags = _parse_cmake_flags(ROOT_DIR / "CMakeLists.txt")
    if cmake_flags != "unknown":
        info["Flags"] = cmake_flags

    return info
