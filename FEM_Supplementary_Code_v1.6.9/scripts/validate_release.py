#!/usr/bin/env python3
"""Check package integrity and available syntax tools; do not run analyses."""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quick_suite_problems(root):
    """Check the same suite inventory consumed by run_MS_validation.R."""
    inventory = root / 'scripts/quick_check_suites.csv'
    if not inventory.is_file():
        return ['Missing scripts/quick_check_suites.csv']
    with inventory.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ['suite', 'script']:
            return ['Quick-check inventory must have suite,script columns']
        rows = list(reader)
    if not rows:
        return ['Quick-check inventory is empty']
    problems = []
    names = set()
    for row in rows:
        name, relative = row['suite'], row['script']
        if not name or not re.fullmatch(r'[A-Za-z0-9_]+', name) or name in names:
            problems.append('Invalid or duplicate suite name: ' + str(name))
        names.add(name)
        path = (root / (relative or '')).resolve()
        if not relative or not path.is_relative_to(root.resolve()) or not path.is_file():
            problems.append('Missing or invalid quick-check script: ' + str(relative))
    return problems


def validate(root, require_r=False):
    checks = []

    def add(name, passed, detail):
        checks.append({'check': name, 'status': 'PASS' if passed else 'FAIL', 'detail': detail})

    paths = {str(p.relative_to(root)): p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    expected = {}
    checksum_file = root / 'SHA256SUMS.txt'
    if checksum_file.exists():
        for line in checksum_file.read_text(encoding='utf-8').splitlines():
            value, name = line.split(maxsplit=1)
            name = name.removeprefix('*').removeprefix('./')
            if name in expected:
                add('checksum duplicate', False, name)
            expected[name] = value
        mismatch = [n for n, h in expected.items() if n not in paths or digest(paths[n]) != h]
        extra = sorted(set(paths) - set(expected) - {'SHA256SUMS.txt'})
        add('SHA256 integrity', not mismatch and not extra, {'entries': len(expected), 'mismatch_or_missing': mismatch, 'unlisted': extra})
    else:
        add('SHA256 integrity', False, 'SHA256SUMS.txt missing')

    manifest = root / 'FILE_MANIFEST.csv'
    if manifest.exists():
        with manifest.open(newline='', encoding='utf-8') as handle:
            rows = list(csv.DictReader(handle))
        names = [row['path'] for row in rows]
        bad = [row['path'] for row in rows if row['path'] not in paths or paths[row['path']].stat().st_size != int(row['size_bytes'])]
        wanted = set(paths) - {'FILE_MANIFEST.csv', 'SHA256SUMS.txt'}
        add('file manifest', not bad and len(names) == len(set(names)) and set(names) == wanted, {'entries': len(rows), 'bad_size_or_missing': bad, 'difference': sorted(set(names) ^ wanted)})
    else:
        add('file manifest', False, 'FILE_MANIFEST.csv missing')

    provenance = root / '06_biodiversity_potential_QRF' / 'SOURCE_PROVENANCE.csv'
    if provenance.exists():
        with provenance.open(newline='', encoding='utf-8') as handle:
            rows = list(csv.DictReader(handle))
        bad = [row['packaged_path'] for row in rows if row['packaged_path'] not in paths or digest(paths[row['packaged_path']]) != row['packaged_sha256']]
        names = [row['packaged_path'] for row in rows]
        add('imported source provenance', bool(rows) and len(names) == len(set(names)) and not bad, {'sources': len(rows), 'mismatch_or_missing': bad})
    else:
        add('imported source provenance', False, 'SOURCE_PROVENANCE.csv missing')
    break_provenance = root / '07_break_year_analysis' / 'SOURCE_PROVENANCE.csv'
    if break_provenance.exists():
        with break_provenance.open(newline='', encoding='utf-8') as handle:
            rows = list(csv.DictReader(handle))
        bad = [r['packaged_path'] for r in rows if r['packaged_path'] not in paths or digest(paths[r['packaged_path']]) != r['packaged_sha256']]
        engine = root / '07_break_year_analysis/upstream/break_year_v5_6_2'
        actual = {str(p.relative_to(root)) for p in engine.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
        recorded = {r['packaged_path'] for r in rows}
        add('break-year upstream provenance', bool(rows) and not bad and actual == recorded and len(rows) == len(recorded), {'files': len(rows), 'mismatch_or_missing': bad, 'inventory_difference': sorted(actual ^ recorded)})
    else:
        add('break-year upstream provenance', False, 'SOURCE_PROVENANCE.csv missing')
    version = (root / 'VERSION').read_text().strip()
    info = json.loads((root / 'PACKAGE_INFO.json').read_text())
    specification = json.loads((root / 'SOURCE_SPECIFICATION.json').read_text())
    add('version', bool(re.fullmatch(r'\d+\.\d+\.\d+', version))
        and f'version: "{version}"' in (root / 'CITATION.cff').read_text()
        and info.get('software_version') == version and specification.get('release') == version,
        'VERSION, CITATION.cff, PACKAGE_INFO.json and SOURCE_SPECIFICATION.json must agree')

    broken_links = []
    for name, path in paths.items():
        if path.suffix.lower() != '.md':
            continue
        for target in re.findall(r'\[[^\]]*\]\(([^\s)]+)\)', path.read_text(encoding='utf-8')):
            if re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:', target) or target.startswith('#'):
                continue
            target = target.split('#', 1)[0]
            if target and not (path.parent / target).exists():
                broken_links.append({'file': name, 'target': target})
    add('local documentation links', not broken_links, broken_links)
    suite_problems = quick_suite_problems(root)
    add('quick-check entry points', not suite_problems, suite_problems)
    required_docs = [info.get('main_guide', 'README.md'), info.get('run_guide', 'README.md'),
                     'DEPENDENCIES.md', 'IMPLEMENTATION_CHOICES.md']
    missing_docs = [name for name in required_docs if name and name not in paths]
    add('documented release interfaces', not missing_docs, missing_docs)
    clutter = [str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and
               (p.suffix in {'.pyc', '.patch', '.zip'} or p.name in {'.DS_Store', '.Rhistory', '.RData'}
                or '__pycache__' in p.parts or 'legacy_pre_v1.6' in p.parts)]
    add('submission layout', not clutter, {'development_artifacts': clutter})

    for suffix, tool, prefix in [('.py', None, None), ('.sh', 'bash', ['-n']), ('.js', 'node', ['--check'])]:
        files = sorted(p for p in paths.values() if p.suffix == suffix)
        failures = []
        executable = shutil.which(tool) if tool else None
        if tool and not executable:
            checks.append({'check': suffix + ' syntax', 'status': 'SKIP', 'detail': tool + ' unavailable'})
            continue
        for file in files:
            try:
                if suffix == '.py':
                    ast.parse(file.read_text(encoding='utf-8'), filename=str(file))
                else:
                    result = subprocess.run([executable, *prefix, str(file)], capture_output=True, text=True)
                    if result.returncode:
                        failures.append({'file': str(file.relative_to(root)), 'error': result.stderr.strip()})
            except Exception as exc:
                failures.append({'file': str(file.relative_to(root)), 'error': str(exc)})
        add(suffix + ' syntax', not failures, {'files': len(files), 'failures': failures})

    rscript = shutil.which('Rscript')
    if rscript:
        result = subprocess.run([rscript, '--vanilla', str(root / 'scripts/check_R_syntax.R'), str(root)], capture_output=True, text=True)
        add('R syntax', result.returncode == 0, (result.stdout + result.stderr).strip())
    else:
        checks.append({'check': 'R syntax', 'status': 'FAIL' if require_r else 'SKIP', 'detail': 'Rscript unavailable; R parsing and execution are not verified'})
    checks.append({'check': 'scientific reproduction', 'status': 'NOT_RUN', 'detail': 'Integrity and syntax do not reproduce scientific results; original data and the analysis environment are required.'})
    return {'version': version, 'ok': not any(x['status'] == 'FAIL' for x in checks), 'checks': checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--require-r', action='store_true', help='Fail rather than skip if Rscript is unavailable.')
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()
    report = validate(args.root.resolve(), args.require_r)
    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        for check in report['checks']:
            print(f"{check['status']}: {check['check']} — {check['detail']}")
    return 0 if report['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
