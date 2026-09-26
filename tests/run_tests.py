"""Regression tests for Ambuild's Poreblazer fork.

    python3 tests/run_tests.py --fork src/poreblazer.exe --upstream UPSTREAM_EXE [--threads 1 4]
    python3 tests/run_tests.py --fork src/poreblazer.exe --write-exact-reference

Runs every case in CASES (upstream's example frameworks and cells built by Ambuild)
and checks that:

  * with the default labelling, the fork's output matches upstream's exactly at every
    thread count: its log (less the line naming the labelling) and every file it writes,
    including nitrogen_network.grd, which holds every grid cube's pore radius;
  * with exact labelling, the output does not depend on the thread count and matches
    tests/reference_exact.json.

Only the Python standard library is needed. --write-exact-reference rewrites the
reference from the fork's exact labelling; review the changes before committing.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CASES_DIR = os.path.join(HERE, "cases")

# (name, directory, cubelet size). The large frameworks use a coarser grid to keep the
# upstream runs short; the output must match at any grid spacing. MOF180 is hexagonal
# (gamma = 120), so it exercises the non-orthorhombic code paths.
CASES = [
    ("HKUST1", os.path.join(ROOT, "HKUST1"), 0.2),
    ("IRMOF1", os.path.join(ROOT, "IRMOF1"), 0.2),
    ("MIL47V", os.path.join(ROOT, "MIL47V"), 0.4),
    ("MOF180", os.path.join(ROOT, "MOF180"), 0.4),
    ("ambuild_20A", os.path.join(CASES_DIR, "ambuild_20A"), 0.2),
    ("ambuild_30A_dense", os.path.join(CASES_DIR, "ambuild_30A_dense"), 0.2),
    ("ambuild_40A_sparse", os.path.join(CASES_DIR, "ambuild_40A_sparse"), 0.2),
]
EXACT_REFERENCE = os.path.join(HERE, "reference_exact.json")
LABELLING_LINE = "Percolation labelling:"


def defaults_dat(case_dir, cubelet, labelling):
    """defaults.dat for a case: its first six values (or the shared ones), the grid
    spacing, and the visualisation line (grid file on, plus the labelling if exact)"""
    source = os.path.join(case_dir, "defaults.dat")
    if not os.path.isfile(source):
        source = os.path.join(CASES_DIR, "defaults.dat")
    with open(source) as f:
        lines = [line.rstrip() for line in f]
    values = lines[:6]
    values[3] = "%g" % cubelet
    values.append("2, 1" if labelling == "exact" else "2")
    return "\n".join(values) + "\n\n"


def run(exe, case, threads, labelling):
    """Run exe on a case; return the log (less the labelling line), the sha256 of every
    file written, and the wall time"""
    name, case_dir, cubelet = case
    rundir = tempfile.mkdtemp(prefix="pb_%s_" % name)
    try:
        # only the inputs: input.dat, the structure it names and the atom types (the
        # example directories also hold outputs of old runs)
        with open(os.path.join(case_dir, "input.dat")) as f:
            xyz = f.readline().split()[0]
        atoms = os.path.join(case_dir, "UFF.atoms")
        if not os.path.isfile(atoms):
            atoms = os.path.join(CASES_DIR, "UFF.atoms")
        for path in (os.path.join(case_dir, "input.dat"), os.path.join(case_dir, xyz), atoms):
            shutil.copy(path, rundir)
        inputs = {"input.dat", xyz, "UFF.atoms"}
        with open(os.path.join(rundir, "defaults.dat"), "w") as f:
            f.write(defaults_dat(case_dir, cubelet, labelling))
        inputs.add("defaults.dat")
        with open(os.path.join(rundir, "input.dat")) as f:
            stdin = f.read()
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        start = time.time()
        proc = subprocess.run([os.path.abspath(exe)], input=stdin, cwd=rundir, env=env,
                              capture_output=True, text=True)
        wall = time.time() - start
        if proc.returncode != 0:
            raise RuntimeError("%s failed on %s (exit %d):\n%s" % (exe, name, proc.returncode,
                                                                   proc.stderr[-2000:]))
        log = [line for line in proc.stdout.splitlines() if LABELLING_LINE not in line]
        files = {}
        for entry in sorted(os.listdir(rundir)):
            if entry in inputs:
                continue
            with open(os.path.join(rundir, entry), "rb") as f:
                files[entry] = hashlib.sha256(f.read()).hexdigest()
        return {"log": log, "files": files, "wall": wall,
                "labelling": [l.split(":", 1)[1].strip() for l in proc.stdout.splitlines()
                              if LABELLING_LINE in l]}
    finally:
        shutil.rmtree(rundir, ignore_errors=True)


def differences(a, b):
    """What differs between two runs' outputs"""
    diff = []
    if set(a["files"]) != set(b["files"]):
        diff.append("files written: %s vs %s" % (sorted(a["files"]), sorted(b["files"])))
    diff += ["%s differs" % f for f in sorted(set(a["files"]) & set(b["files"])) if a["files"][f] != b["files"][f]]
    if a["log"] != b["log"]:
        for i, (x, y) in enumerate(zip(a["log"], b["log"])):
            if x != y:
                diff.append("log line %d: %r vs %r" % (i + 1, x.strip(), y.strip()))
                break
        else:
            diff.append("log lengths %d vs %d" % (len(a["log"]), len(b["log"])))
    return diff


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fork", required=True, help="the fork's poreblazer.exe")
    parser.add_argument("--upstream", help="upstream 3.0.5's poreblazer.exe, to compare the default labelling with")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--cases", nargs="+", help="run only these cases")
    parser.add_argument("--write-exact-reference", action="store_true")
    args = parser.parse_args()
    cases = [c for c in CASES if not args.cases or c[0] in args.cases]

    failures = []
    if args.write_exact_reference:
        reference = {}
        for case in cases:
            r = run(args.fork, case, args.threads[0], "exact")
            reference[case[0]] = {"log": r["log"], "files": r["files"]}
            print("%-20s exact  %6.1f s  %d files" % (case[0], r["wall"], len(r["files"])), flush=True)
        with open(EXACT_REFERENCE, "w") as f:
            json.dump(reference, f, indent=1, sort_keys=True)
            f.write("\n")
        print("wrote", EXACT_REFERENCE)
        return 0

    with open(EXACT_REFERENCE) as f:
        exact_reference = json.load(f)
    for case in cases:
        name = case[0]
        upstream = run(args.upstream, case, 1, "poreblazer") if args.upstream else None
        if upstream:
            print("%-20s upstream           %6.1f s" % (name, upstream["wall"]), flush=True)
        for threads in args.threads:
            fork = run(args.fork, case, threads, "poreblazer")
            status = "ran"
            if fork["labelling"] != ["poreblazer"]:
                failures.append("%s: default run reported labelling %s" % (name, fork["labelling"]))
            if upstream:
                diff = differences(fork, upstream)
                status = "matches upstream" if not diff else "DIFFERS: " + "; ".join(diff)
                if diff:
                    failures.append("%s default, %d threads: %s" % (name, threads, "; ".join(diff)))
            print("%-20s default %2d threads %6.1f s  %s" % (name, threads, fork["wall"], status), flush=True)
        for threads in args.threads:
            exact = run(args.fork, case, threads, "exact")
            if exact["labelling"] != ["exact"]:
                failures.append("%s: exact run reported labelling %s" % (name, exact["labelling"]))
            ref = exact_reference.get(name)
            diff = differences(exact, ref) if ref else ["no reference"]
            if diff:
                failures.append("%s exact, %d threads: %s" % (name, threads, "; ".join(diff)))
            print("%-20s exact   %2d threads %6.1f s  %s" % (
                name, threads, exact["wall"], "matches reference" if not diff else "DIFFERS: " + "; ".join(diff)),
                flush=True)

    if failures:
        print("\n%d failure(s):" % len(failures))
        for f in failures:
            print("  " + f)
        return 1
    print("\nAll cases passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
