#!/usr/bin/env python2.7
"""Capture reference outputs of the Python 2 biofrills code base.

This script must run under Python 2.7 with Biopython 1.76, with biofrills/
and scripts/ unmodified since the commit recorded in py27/versions.json; see
README.md beside it for the environment.  It writes JSON files under
tests/data/py27/:

* library/<input>.json -- library functions applied to each aligned input,
* library/tie-batch.json -- consensus of many small tie-heavy alignments,
* library/chisqprob.json -- stats.chisq.chisqprob over a grid of (x, df),
* scripts/<case>.json -- stdout, stderr, exit status and written files of
  each command-line script case.

Every function call records either {"value": ...} or {"error": {"type": ...,
"message": ...}}, so that failures of the original code are part of the
reference.  Floats are written with Python 2.7's repr, the shortest decimal
string that reads back to the identical double.

seq-purge and seq-diff import biocma, whose current version is Python 3
only.  The last Python 2 commit of biocma, BIOCMA_COMMIT, is extracted from
the biocma repository beside this one (or from $BIOCMA_REPO) into a scratch
directory and put on the scripts' module path.

Usage, from the repository root:

    <python2.7> tests/data/capture_py27.py
"""

import io
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
INPUTS = os.path.join(HERE, 'inputs')
OUT = os.path.join(HERE, 'py27')

sys.path.insert(0, REPO)

import Bio
import networkx
import numpy
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

from biofrills import alnutils, consensus, pairutils
from biofrills.stats import chisq

ALIGNED = [
    'cdc2-domain.aln.fasta',
    'scttlre-domain.aln.fasta',
    'fikk.aln.fasta',
    'synth-edge.aln.fasta',
    'synth-ties.aln.fasta',
    'synth-mixedcase.aln.fasta',
    'synth-single.aln.fasta',
    'synth-x.aln.fasta',
]

SCALINGS = ['none', 'sum1', 'max1', 'avg1', 'andy']

# Sequence weights passed to aa_counts, aa_frequencies and consensus.  The
# names are recorded in the fixtures; a consumer of the fixtures rebuilds the
# values from them.
#   None            -- the argument's default
#   True            -- aa_counts computes sequence_weights(aln) itself
#   "seqweights"    -- sequence_weights(aln, 'none'), as cladecompare passes
#   "unit"          -- [1] * len(aln), cladecompare's unweighted mode
WEIGHT_SPECS = [None, True, 'seqweights', 'unit']

BIOCMA_COMMIT = 'eac0c57eb83a9498e53ccdeb9cbc3fe21a5826a7'


def build_weights(aln, spec):
    """The sequence weights that a WEIGHT_SPECS name stands for."""
    if spec == 'seqweights':
        return alnutils.sequence_weights(aln, 'none')
    if spec == 'unit':
        return [1] * len(aln)
    return spec


def outcome(func, *args, **kwargs):
    """Call func; return {"value": result} or {"error": {...}}."""
    try:
        return OrderedDict([('value', func(*args, **kwargs))])
    except Exception as exc:
        return OrderedDict([('error', OrderedDict([
            ('type', type(exc).__name__),
            ('message', str(exc)),
        ]))])


def records_table(records):
    """SeqRecords as [id, description, sequence] rows."""
    return [[rec.id, rec.description, str(rec.seq)] for rec in records]


def plain_dict(mapping):
    """A Counter or defaultdict as a key-sorted dict, for JSON."""
    return OrderedDict(sorted(mapping.items()))


def capture_alignment(fname):
    """Library functions applied to one aligned input file."""
    aln = AlignIO.read(os.path.join(INPUTS, fname), 'fasta')
    result = OrderedDict()
    result['input'] = 'inputs/' + fname
    result['ids'] = [rec.id for rec in aln]

    result['sequence_weights'] = [
        OrderedDict([('kwargs', {'scaling': scaling}),
                     ('result', outcome(alnutils.sequence_weights, aln,
                                        scaling))])
        for scaling in SCALINGS]

    for func_name in ('aa_counts', 'aa_frequencies'):
        func = getattr(alnutils, func_name)
        cases = []
        for spec in WEIGHT_SPECS:
            for gap_chars in ('-.', '-.X'):
                def call(spec=spec, gap_chars=gap_chars):
                    weights = build_weights(aln, spec)
                    return plain_dict(func(aln, weights=weights,
                                           gap_chars=gap_chars))
                cases.append(OrderedDict([
                    ('weights', spec),
                    ('kwargs', {'gap_chars': gap_chars}),
                    ('result', outcome(call)),
                ]))
        result[func_name] = cases

    cases = []
    for simple in (False, True):
        for spec in ([None] if simple else [None, 'seqweights', 'unit']):
            for trim_ends in (True, False):
                for gap_threshold in (0.5, 0.6, 0.8):
                    def call(spec=spec, simple=simple, trim_ends=trim_ends,
                             gap_threshold=gap_threshold):
                        return consensus.consensus(
                            aln, weights=build_weights(aln, spec),
                            gap_threshold=gap_threshold, simple=simple,
                            trim_ends=trim_ends)
                    cases.append(OrderedDict([
                        ('weights', spec),
                        ('kwargs', OrderedDict([
                            ('simple', simple),
                            ('trim_ends', trim_ends),
                            ('gap_threshold', gap_threshold)])),
                        ('result', outcome(call)),
                    ]))
    result['consensus'] = cases

    cases = []
    for threshold in (0.4, 0.5, 0.9):
        for weights in (None, False):
            def call(threshold=threshold, weights=weights):
                return records_table(alnutils.blocks(aln, threshold, weights))
            cases.append(OrderedDict([
                ('weights', weights),
                ('kwargs', {'threshold': threshold}),
                ('result', outcome(call)),
            ]))
    result['blocks'] = cases

    result['remove_empty_cols'] = outcome(
        lambda: records_table(alnutils.remove_empty_cols(aln)))

    # All pairs i < j, as rows [i, j, outcome].
    seqs = [str(rec.seq) for rec in aln]
    for func_name in ('score_pairwise', 'identity_abs', 'identity_avg'):
        func = getattr(pairutils, func_name)
        result[func_name] = [
            [i, j, outcome(func, seqs[i], seqs[j])]
            for i in range(len(seqs))
            for j in range(i + 1, len(seqs))]
    return result


def capture_tie_batch():
    """Consensus of each small alignment in tie-batch.fasta."""
    groups = OrderedDict()
    for rec in SeqIO.parse(os.path.join(INPUTS, 'tie-batch.fasta'), 'fasta'):
        groups.setdefault(rec.id.split('_')[0], []).append(rec)
    result = OrderedDict()
    result['input'] = 'inputs/tie-batch.fasta'
    result['sequence_weights_kwargs'] = {'scaling': 'avg1'}
    result['consensus_kwargs'] = [OrderedDict([('simple', simple),
                                               ('trim_ends', trim_ends)])
                                  for simple in (False, True)
                                  for trim_ends in (True, False)]
    alignments = []
    for key, recs in groups.items():
        aln = MultipleSeqAlignment(recs)
        alignments.append(OrderedDict([
            ('id', key),
            ('sequence_weights',
             outcome(alnutils.sequence_weights, aln, 'avg1')),
            ('consensus', [outcome(consensus.consensus, aln, **kw)
                           for kw in result['consensus_kwargs']]),
        ]))
    result['alignments'] = alignments
    return result


def capture_chisqprob():
    """chisqprob over a grid of x and df, as [x, df, outcome] rows."""
    xs = [-1.0, 0.0, 1e-6, 0.001, 0.1, 0.5, 1.0, 1.5, 2.0, 3.84, 5.0, 7.5,
          10.0, 18.0, 19.0, 20.0, 30.0, 50.0, 100.0, 300.0, 1000.0]
    dfs = [0, 1, 2, 3, 4, 5, 10, 19, 20, 50, 100, 0.5, 2.5]
    return OrderedDict([
        ('function', 'biofrills.stats.chisq.chisqprob'),
        ('rows', [[x, df, outcome(chisq.chisqprob, x, df)]
                  for x in xs for df in dfs]),
    ])


# Script cases: (case name, script, argv, stdin input file or None).  A case
# is named <script>.<input>[.<options>], options abbreviated as in argv.
# Arguments name files in inputs/, which are copied into a scratch working
# directory, so recorded output refers to them by bare file name.
SCRIPT_CASES = [
    ('seq-blocks.cdc2', 'seq-blocks', ['cdc2-domain.aln.fasta'], None),
    ('seq-blocks.cdc2.t09', 'seq-blocks',
     ['-t', '0.9', 'cdc2-domain.aln.fasta'], None),
    ('seq-blocks.cdc2.w.t09', 'seq-blocks',
     ['-w', '-t', '0.9', 'cdc2-domain.aln.fasta'], None),
    ('seq-blocks.edge', 'seq-blocks', ['synth-edge.aln.fasta'], None),
    ('seq-blocks.edge.w.t08', 'seq-blocks',
     ['-w', '-t', '0.8', 'synth-edge.aln.fasta'], None),

    ('seq-split.s10', 'seq-split', ['-s', '10', 'cdc2-domain.fasta'], None),
    ('seq-split.s1e1', 'seq-split', ['-s', '1e1', 'cdc2-domain.fasta'], None),
    ('seq-split.b3000', 'seq-split', ['-b', '3000', 'cdc2-domain.fasta'],
     None),
    ('seq-split.ns3', 'seq-split', ['-ns', '3', 'cdc2-domain.fasta'], None),
    ('seq-split.nb3', 'seq-split', ['-nb', '3', 'cdc2-domain.fasta'], None),

    ('seq-purge.h', 'seq-purge', ['-h', 'dupes.fasta'], None),
    ('seq-purge.i', 'seq-purge', ['-i', 'dupes.fasta'], None),
    ('seq-purge.s', 'seq-purge', ['-s', 'dupes.fasta'], None),
    ('seq-purge.a', 'seq-purge', ['-a', 'dupes.fasta'], None),
    ('seq-purge.b', 'seq-purge', ['-b', 'dupes.fasta'], None),
    ('seq-purge.c', 'seq-purge', ['-c', 'dupes.fasta'], None),
    ('seq-purge.s.stdin', 'seq-purge', ['-s'], 'dupes.fasta'),
    ('seq-purge.fikk-cma.s', 'seq-purge', ['-s', 'fikk.cma'], None),
    ('seq-purge.fikk-cma.c', 'seq-purge', ['-c', 'fikk.cma'], None),

    ('seq-diff.cdc2.i', 'seq-diff',
     ['-i', 'cdc2-domain.fasta', 'diff-query.fasta'], None),
    ('seq-diff.cdc2.h', 'seq-diff',
     ['-h', 'cdc2-domain.fasta', 'diff-query.fasta'], None),
    ('seq-diff.cdc2.b', 'seq-diff',
     ['-b', 'cdc2-domain.fasta', 'diff-query.fasta'], None),
    ('seq-diff.cdc2.a', 'seq-diff',
     ['-a', 'cdc2-domain.fasta', 'diff-query.fasta'], None),
    ('seq-diff.cdc2.s', 'seq-diff',
     ['-s', 'cdc2-domain.fasta', 'diff-query.fasta'], None),
    ('seq-diff.query.i', 'seq-diff',
     ['-i', 'diff-query.fasta', 'cdc2-domain.fasta'], None),
    ('seq-diff.cdc2.i.stdin', 'seq-diff', ['-i', 'diff-query.fasta'],
     'cdc2-domain.fasta'),
    ('seq-diff.fikk-cma.i', 'seq-diff',
     ['-i', 'fikk.cma', 'fikk-subset.cma'], None),

    ('seq-identity.cdc2', 'seq-identity', ['cdc2-domain.aln.fasta'], None),
    ('seq-identity.edge', 'seq-identity', ['synth-edge.aln.fasta'], None),

    ('seq-extract.cdc2', 'seq-extract',
     ['cdc2-domain.aln.fasta', '1-5,10-20,34,35,240-'], None),
    ('seq-extract.edge.stdin', 'seq-extract', ['-3,7'],
     'synth-edge.aln.fasta'),

    ('seq-orfs.default', 'seq-orfs', ['dna.fasta'], None),
    ('seq-orfs.m20', 'seq-orfs', ['-m', '20', 'dna.fasta'], None),

    ('seq-sort.default', 'seq-sort', ['cdc2-domain.fasta'], None),
    ('seq-sort.a', 'seq-sort', ['-a', 'cdc2-domain.fasta'], None),
    ('seq-sort.s', 'seq-sort', ['-s', 'cdc2-domain.fasta'], None),
    ('seq-sort.l', 'seq-sort', ['-l', 'cdc2-domain.fasta'], None),

    ('seq-convert.fasta-tab', 'seq-convert', ['fasta', 'tab'],
     'cdc2-domain.fasta'),
    ('seq-convert.fasta-clustal', 'seq-convert', ['fasta', 'clustal'],
     'synth-edge.aln.fasta'),

    ('seq-composition.cdc2', 'seq-composition', ['cdc2-domain.fasta'], None),
    ('seq-composition.mixedcase', 'seq-composition',
     ['synth-mixedcase.aln.fasta'], None),

    ('seq-dupes.dupes', 'seq-dupes', ['dupes.fasta'], None),

    ('seq-intersect.plain', 'seq-intersect',
     ['intersect-query.fasta', 'cdc2-domain.fasta'], None),
    ('seq-intersect.r', 'seq-intersect',
     ['-r', 'intersect-query.fasta', 'cdc2-domain.fasta'], None),
    ('seq-intersect.stdin', 'seq-intersect', ['cdc2-domain.fasta'],
     'intersect-query.fasta'),

    ('seq-pick.f5', 'seq-pick', ['-f', '5', 'cdc2-domain.fasta'], None),
    ('seq-pick.l5', 'seq-pick', ['-l', '5', 'cdc2-domain.fasta'], None),
    ('seq-pick.s5', 'seq-pick', ['-s', '5', 'cdc2-domain.fasta'], None),
    ('seq-pick.c5', 'seq-pick', ['-c', '5', 'cdc2-domain.aln.fasta'], None),
    ('seq-pick.k3', 'seq-pick', ['-k', '3', 'cdc2-domain.aln.fasta'], None),
]


def run_script_case(script, argv, stdin_name, biocma_dir):
    """Run one script case in a scratch copy of inputs/; record the results."""
    workdir = tempfile.mkdtemp(prefix='biofrills-py27-')
    try:
        for fname in os.listdir(INPUTS):
            shutil.copy(os.path.join(INPUTS, fname), workdir)
        before = set(os.listdir(workdir))
        env = dict(os.environ, PYTHONPATH=REPO + os.pathsep + biocma_dir,
                   PYTHONHASHSEED='0')
        stdin = (open(os.path.join(workdir, stdin_name), 'rb')
                 if stdin_name else open(os.devnull, 'rb'))
        try:
            proc = subprocess.Popen(
                [sys.executable, os.path.join(REPO, 'scripts', script)]
                + argv,
                cwd=workdir, env=env, stdin=stdin,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
        finally:
            stdin.close()
        written = OrderedDict()
        for fname in sorted(set(os.listdir(workdir)) - before):
            with open(os.path.join(workdir, fname), 'rb') as handle:
                written[fname] = handle.read().decode('utf-8')
        # Tracebacks name absolute paths; make them machine-independent.
        stderr = (stderr.decode('utf-8')
                  .replace(os.path.realpath(workdir), '<workdir>')
                  .replace(workdir, '<workdir>')
                  .replace(REPO, '<repo>')
                  .replace(biocma_dir, '<biocma>')
                  .replace(sys.prefix, '<python-prefix>'))
        return OrderedDict([
            ('script', script),
            ('argv', argv),
            ('stdin', stdin_name),
            ('returncode', proc.returncode),
            ('stdout', stdout.decode('utf-8')),
            ('stderr', stderr),
            ('written_files', written),
        ])
    finally:
        shutil.rmtree(workdir)


def write_json(relpath, data):
    """Write data under OUT as indented ASCII JSON."""
    path = os.path.join(OUT, relpath)
    if not os.path.isdir(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    text = json.dumps(data, indent=1, separators=(',', ': '),
                      ensure_ascii=True)
    with io.open(path, 'w', encoding='ascii', newline='\n') as handle:
        handle.write(text.decode('ascii') + u'\n')


def extract_biocma(destdir):
    """Extract the biocma package at BIOCMA_COMMIT into destdir."""
    biocma_repo = os.environ.get(
        'BIOCMA_REPO', os.path.join(os.path.dirname(REPO), 'biocma'))
    tarpath = os.path.join(destdir, 'biocma.tar')
    subprocess.check_call(['git', '-C', biocma_repo, 'archive', '-o', tarpath,
                           BIOCMA_COMMIT, 'biocma'])
    with tarfile.open(tarpath) as archive:
        archive.extractall(destdir)
    os.remove(tarpath)


def main():
    """Check the environment, then rewrite every fixture under OUT."""
    assert sys.version_info[:2] == (2, 7), sys.version
    assert not sys.flags.hash_randomization, "unset PYTHONHASHSEED or -R"
    assert pairutils.score_pairwise.__module__ == 'biofrills.pairutils', \
        "compiled cpairutils is active; capture the pure-Python path"
    assert chisq.chisqprob.__module__ == 'biofrills.stats.chisq', \
        "scipy's chisqprob is active; capture the pure-Python path"

    # Extract biocma first, so that a missing repository fails before the
    # existing fixtures are deleted.
    biocma_dir = tempfile.mkdtemp(prefix='biocma-py2-')
    try:
        extract_biocma(biocma_dir)
        if os.path.isdir(OUT):
            shutil.rmtree(OUT)
        for fname in ALIGNED:
            write_json(os.path.join('library',
                                    fname[:-len('.fasta')] + '.json'),
                       capture_alignment(fname))
        write_json(os.path.join('library', 'tie-batch.json'),
                   capture_tie_batch())
        write_json(os.path.join('library', 'chisqprob.json'),
                   capture_chisqprob())
        for case, script, argv, stdin_name in SCRIPT_CASES:
            write_json(os.path.join('scripts', case + '.json'),
                       run_script_case(script, argv, stdin_name, biocma_dir))
    finally:
        shutil.rmtree(biocma_dir)

    sources = ['biofrills', 'scripts', 'setup.py']
    # The last commit that changed the code under test, rather than HEAD,
    # so that regenerating on a later commit reproduces this file too.
    sources_commit = subprocess.check_output(
        ['git', 'log', '-1', '--format=%H', '--'] + sources,
        cwd=REPO).decode('ascii').strip()
    sources_dirty = subprocess.check_output(
        ['git', 'status', '--porcelain', '--'] + sources,
        cwd=REPO).decode('ascii').strip()
    write_json('versions.json', OrderedDict([
        ('python', sys.version.split()[0]),
        ('implementation', platform.python_implementation()),
        ('platform', platform.platform()),
        ('machine', platform.machine()),
        ('biopython', Bio.__version__),
        ('numpy', numpy.__version__),
        ('networkx', networkx.__version__),
        ('biocma_commit', BIOCMA_COMMIT),
        ('biofrills_sources_commit', sources_commit),
        ('biofrills_sources_modified', bool(sources_dirty)),
    ]))


if __name__ == '__main__':
    main()
