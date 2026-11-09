"""Characterization (golden master) tests for the whole `debreak` pipeline.

They pin the CURRENT behaviour, not necessarily the correct one, so that a
refactor can be checked for "no output changed". If a change is intended,
regenerate the snapshots with:  pytest --update-golden
"""
import os
import shutil
import subprocess
import sys

import pytest

from conftest import REPO

LAUNCHER = os.path.join(REPO, "tests", "run_debreak.py")
GOLDEN = os.path.join(REPO, "tests", "golden")
BAM = "testdata/test_read.bam"  # relative on purpose: it ends up in the VCF header
REF = "testdata/test_ref.fa"

# Header lines that legitimately differ between runs.
VOLATILE_HEADERS = ("##fileDate=", "##CommandLine=")


def _normalize_vcf(text):
    return "".join(
        line + "\n"
        for line in text.splitlines()
        if not line.startswith(VOLATILE_HEADERS)
    )


def _sorted_lines(text):
    return "".join(line + "\n" for line in sorted(text.splitlines()))


def _run(outdir, *extra):
    cmd = [sys.executable, LAUNCHER, "--bam", BAM, "-o", str(outdir) + "/", "-t", "4", *extra]
    proc = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc


def _collect(outdir):
    """Outputs worth pinning, normalized to be run-independent."""
    out = {}
    with open(os.path.join(outdir, "debreak.vcf")) as f:
        out["debreak.vcf"] = _normalize_vcf(f.read())
    with open(os.path.join(outdir, "debreak-allsv-merged-final")) as f:
        out["debreak-allsv-merged-final"] = _sorted_lines(f.read())
    raw_dir = os.path.join(outdir, "sv_raw_calls")
    raw = "".join(open(os.path.join(raw_dir, n)).read() for n in sorted(os.listdir(raw_dir)))
    out["sv_raw_calls.txt"] = _sorted_lines(raw)
    with open(os.path.join(outdir, "log.txt")) as f:
        out["log.txt"] = f.read()
    return out


def _check_golden(name, outputs, update, required=True):
    gdir = os.path.join(GOLDEN, name)
    if update:
        shutil.rmtree(gdir, ignore_errors=True)
        os.makedirs(gdir)
        for fname, content in outputs.items():
            with open(os.path.join(gdir, fname), "w") as f:
                f.write(content)
        pytest.skip("golden files updated: " + gdir)
    if not os.path.isdir(gdir) and not required:
        pytest.skip(f"no golden files for '{name}' yet: generate them on Linux with "
                    "the external tools (CI: run workflow with update_golden=true)")
    assert os.path.isdir(gdir), "no golden files, run: pytest --update-golden"
    for fname, content in outputs.items():
        with open(os.path.join(gdir, fname)) as f:
            expected = f.read()
        assert content == expected, f"{name}/{fname} differs from golden"


SCENARIOS = {
    "default": [],
    "no_genotype": ["--no_genotype"],
    "tumor": ["--tumor"],
    "min_support_10": ["-m", "10"],
    "depth_30_maxcov_80": ["-d", "30", "--maxcov", "80"],
    "min_size_500": ["--min_size", "500"],
}


@pytest.mark.parametrize("name", SCENARIOS)
def test_pipeline_matches_golden(name, tmp_path, update_golden):
    _run(tmp_path / "out", *SCENARIOS[name])
    _check_golden(name, _collect(tmp_path / "out"), update_golden)


def test_default_run_finds_all_six_simulated_svs(tmp_path):
    """Sanity check independent of snapshots: the test BAM has 6 embedded SVs."""
    _run(tmp_path / "out")
    with open(tmp_path / "out" / "debreak.vcf") as f:
        records = [l.split("\t") for l in f if not l.startswith("#")]
    svtypes = sorted(r[7].split("SVTYPE=")[1].split(";")[0] for r in records)
    assert svtypes == ["DEL", "DEL", "DUP", "INS", "INS", "INV"]


def test_missing_input_exits_with_error(tmp_path):
    proc = subprocess.run(
        [sys.executable, LAUNCHER, "-o", str(tmp_path / "out")],
        cwd=REPO, capture_output=True, text=True,
    )
    assert proc.returncode == 1
    assert "No input file given" in proc.stdout


def test_rescue_dup_without_ref_logs_warning(tmp_path):
    _run(tmp_path / "out", "--rescue_dup")
    log = (tmp_path / "out" / "log.txt").read_text()
    assert "Reference genome is required for duplication rescue" in log


# --- modules needing external binaries (Linux / bioconda env) --------------

def _have(*tools):
    return all(shutil.which(t) for t in tools)


@pytest.mark.skipif(not _have("minimap2"), reason="needs minimap2")
def test_rescue_dup_matches_golden(tmp_path, update_golden):
    _run(tmp_path / "out", "--rescue_dup", "-r", REF)
    _check_golden("rescue_dup", _collect(tmp_path / "out"), update_golden, required=False)


@pytest.mark.skipif(not _have("minimap2", "wtdbg2", "wtpoa-cns", "bsalign"),
                    reason="needs minimap2, wtdbg2, bsalign")
def test_full_pipeline_matches_golden(tmp_path, update_golden):
    _run(tmp_path / "out", "--poa", "--rescue_large_ins", "--rescue_dup", "-r", REF)
    _check_golden("full", _collect(tmp_path / "out"), update_golden, required=False)
