"""vegasgen CLI end-to-end: argument parsing, model run, CSV/JSON output."""
import csv
import json
import sys

import pytest

from VegasAfterglow import cli


def _run(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["vegasgen"] + argv)
    cli.main()


def test_csv_output(tmp_path, monkeypatch):
    """CSV output has a column header plus exactly num_t data rows, and the first data row has at least two columns with positive time and flux values."""
    out = tmp_path / "lc.csv"
    _run(monkeypatch, ["--nu", "1e17", "--num_t", "16", "--t_max", "1e6",
                       "--output", str(out), "--format", "csv"])
    assert out.exists()
    with open(out) as f:
        rows = [r for r in csv.reader(f) if r and not r[0].startswith("#")]
    assert len(rows) == 17  # column header + 16 time points
    header, first = rows[0], rows[1]
    assert len(header) == len(first) >= 2
    assert float(first[0]) > 0  # time column
    assert float(first[1]) > 0  # flux column


def test_json_output(tmp_path, monkeypatch):
    """JSON output parses to a dict containing light-curve content (numeric values, a flux field, or a time key)."""
    out = tmp_path / "lc.json"
    _run(monkeypatch, ["--nu", "1e17", "--num_t", "12", "--t_max", "1e6",
                       "--output", str(out), "--format", "json"])
    data = json.loads(out.read_text())
    assert isinstance(data, dict)
    payload = json.dumps(data)
    assert "1e" in payload or "flux" in payload.lower() or "t" in data


def test_named_band_frequency(tmp_path, monkeypatch):
    """A named band ("XRT") is accepted for --nu and the run still writes the output file."""
    out = tmp_path / "lc_band.csv"
    _run(monkeypatch, ["--nu", "XRT", "--num_t", "10", "--t_max", "1e6",
                       "--output", str(out), "--format", "csv"])
    assert out.exists()


def test_bad_frequency_rejected(monkeypatch):
    """An unparseable --nu value aborts with SystemExit or ArgumentTypeError rather than running the model."""
    import argparse

    # Raised as ArgumentTypeError from parse_frequencies (post-parse_args), not
    # converted to a clean SystemExit by argparse; accept either.
    with pytest.raises((SystemExit, argparse.ArgumentTypeError)):
        _run(monkeypatch, ["--nu", "not-a-frequency", "--num_t", "8"])
