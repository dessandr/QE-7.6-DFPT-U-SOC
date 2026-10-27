"""Compare two soc_u_port_regress.sh output directories.

Dynamical-matrix files are compared number by number. For *.out files only
lines with physical results are compared (energies, frequencies, dielectric
tensor, effective charges, Hubbard parameters, occupations), so timings and
dates are ignored. Usage: python3 -I soc_u_port_compare.py REF_DIR NEW_DIR
Exit status 1 if any compared value differs or a file is missing.
"""
import pathlib
import re
import sys

NUM = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')
RESULT_LINE = re.compile(
    r'^!|total energy|Fermi energy|highest occupied|freq|omega\(|'
    r'dielectric|Effective charges|^\s+\(\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\)|'
    r'Px\s|Py\s|Pz\s|Ex\s|Ey\s|Ez\s|atom\s+\d+\s+[A-Za-z]+\s+Mixed|'
    r'Hubbard_U|U\(\s*\d|chi0|chi\s|Tr\[ns|eigenvalues|Atom:\s*\d|'
    r'total magnetization|absolute magnetization|magnetization\s+=',
    re.IGNORECASE)


def numbers(text, results_only):
    vals = []
    for line in text.splitlines():
        if results_only and not RESULT_LINE.search(line):
            continue
        vals += [float(x.replace('D', 'E').replace('d', 'e')) for x in NUM.findall(line)]
    return vals


def compare(ref, new, results_only):
    a = numbers(ref.read_text(errors='replace'), results_only)
    b = numbers(new.read_text(errors='replace'), results_only)
    if len(a) != len(b):
        return False, f'count {len(a)} vs {len(b)}'
    if not a:
        return False, 'no result values found (failed run?)'
    diff = max((abs(x - y) for x, y in zip(a, b)), default=0.0)
    scale = max((abs(x) for x in a), default=0.0)
    return diff == 0.0, f'n={len(a):6d} max|ref|={scale:.3e} max|diff|={diff:.3e}'


def main():
    ref_root, new_root = map(pathlib.Path, sys.argv[1:3])
    ok_all = True
    for ref in sorted(p for p in ref_root.rglob('*') if p.is_file()):
        rel = ref.relative_to(ref_root)
        name = rel.name
        if name == 'BINARIES.sha256' or '/out/' in f'/{rel}':
            continue
        is_dyn = '.dyn' in name or name.endswith('.xml') and 'dyn' in name
        if not (is_dyn or name.endswith('.out')):
            continue
        new = new_root / rel
        if not new.exists():
            print(f'MISSING   {rel}')
            ok_all = False
            continue
        same, info = compare(ref, new, results_only=not is_dyn)
        ok_all &= same
        print(f"{'IDENTICAL' if same else 'DIFFERENT'} {str(rel):45s} {info}")
    print('ALL IDENTICAL' if ok_all else 'DIFFERENCES FOUND')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
