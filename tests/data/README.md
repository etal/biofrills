# Reference outputs from the Python 2 code base

The files under `py27/` record what biofrills computed on Python 2.7, before the port to Python 3.
They are the reference against which the port is checked: a ported function should reproduce these values, or the difference should be explained by a documented, deliberate change.

## Producing versions

`py27/versions.json` is written by the capture script and is the authoritative record.
At capture time it read:

| Component | Version |
|---|---|
| Python | 2.7.15 (CPython, conda-forge build) |
| Platform | Darwin 25.6.0, x86_64 under Rosetta 2 on Apple silicon |
| Biopython | 1.76 |
| NumPy | 1.16.5 |
| NetworkX | 2.2 (used only by `seq-pick`) |
| biocma | commit `eac0c57`, the last Python 2 commit (used only by `seq-purge` and `seq-diff`) |
| biofrills | sources (`biofrills/`, `scripts/`, `setup.py`) as of commit `36684bb`, unmodified |

The compiled `cpairutils` extension was not built and SciPy was not installed, so `pairutils` and `stats.chisq.chisqprob` ran their pure-Python implementations.
The capture script asserts both conditions.
The Cython `score_pairwise` declares its accumulator and return value as single-precision `float`, but every score recorded here is a sum of integers and half-integers small enough to be exact in single precision, so the compiled extension would have returned the same values.

## Regenerating

Later commits port and change the biofrills sources, so the capture runs in a checkout whose sources are unchanged since commit `36684bb`, such as commit `615d27c`, which added this directory.
A worktree beside this repository keeps the default locations of the `biocma` and `cladecompare` repositories valid:

```sh
git worktree add ../biofrills-py27 615d27c
python3 ../biofrills-py27/tests/data/make_inputs.py   # optional: rewrites the worktree's inputs/
```

The `py27-oracle` pixi environment defined in `pyproject.toml` pins the producing versions listed above.
It is built for osx-64, as the capture was, and runs under Rosetta on Apple silicon; it is not available on Linux.
From this repository, the following command runs the worktree's copy of `capture_py27.py`:

```sh
pixi run -e py27-oracle capture
```

Pixi warns that the environment requires an x86_64 processor; with the committed `pixi.lock` the installed packages run under Rosetta regardless.

Without pixi, or on Linux, conda can build an equivalent environment; conda-forge publishes no Python 2.7 build for osx-arm64, but its osx-64 build runs under Rosetta:

```sh
CONDA_SUBDIR=osx-64 conda create -p /tmp/bf-py27 --override-channels -c conda-forge \
    python=2.7 biopython=1.76 numpy networkx
```

On Linux, the same command without `CONDA_SUBDIR` gives a native linux-64 environment.
Then run the worktree's copy of the capture with that environment's interpreter:

```sh
/tmp/bf-py27/bin/python ../biofrills-py27/tests/data/capture_py27.py
```

In the worktree, `git status --short tests/data/py27` then lists every fixture file the capture changed; it lists none when the capture reproduces the committed files.
A capture on another platform or with another Python 2.7 release also changes `versions.json`, which records both.

`make_inputs.py` reads the test files of the `cladecompare` and `biocma` repositories beside this one (or as given by `--cladecompare` and `--biocma`) and writes the committed files under `inputs/`; rerunning it reproduces them byte for byte.
`capture_py27.py` extracts biocma at commit `eac0c57` from the `biocma` repository beside this one (or from `$BIOCMA_REPO`) into a scratch directory, then deletes and rewrites `py27/`.
It requires hash randomization to be off, which is the Python 2 default; leave `PYTHONHASHSEED` unset or set it to `0`.
Two consecutive captures produce identical files.

## Inputs

| File | Source | Used for |
|---|---|---|
| `cdc2-domain.aln.fasta` | cladecompare `test/cdc2-domain.fasta_aln.cma` | library functions, scripts |
| `scttlre-domain.aln.fasta` | cladecompare `test/scttlre-domain.fasta_aln.cma` | library functions |
| `fikk.aln.fasta` | biocma `test/fikk.cma` | library functions |
| `cdc2-domain.fasta` | cladecompare `test/cdc2-domain.fasta` | scripts on unaligned sequences |
| `fikk.cma` | biocma `test/fikk.cma` | the CMA branches of `seq-purge` and `seq-diff` |
| `fikk-subset.cma` | the first 15 sequences of `fikk.cma` | `seq-diff` |
| `synth-*.aln.fasta` | written by `make_inputs.py` | edge cases, described in that script |
| `tie-batch.fasta` | random, seeded | consensus tie-breaking |
| `dna.fasta` | random, seeded | `seq-orfs` |
| `dupes.fasta` | built from scttlre-domain sequences | `seq-purge`, `seq-dupes` |
| `intersect-query.fasta` | identifiers from `cdc2-domain.fasta` | `seq-intersect` |
| `diff-query.fasta` | records from `cdc2-domain.fasta`, some altered, and one new record | `seq-diff` |

The `.aln.fasta` files are the MAPGAPS alignments of the cladecompare and biocma test sets, reduced to their consensus columns by the rule of `cladecomparelib.core.cma_blocks`, which drops lowercase insert residues.
That is the form in which cladecompare passes alignments to biofrills.
The files `cdc2-domain.fasta` and `fikk.cma` are copied unchanged.
The FASTA files `cdc2-domain.fasta`, `scttlre-domain.fasta` and `fikk-full.fasta` in those repositories are unaligned.

`tie-batch.fasta` holds 200 small alignments over a five-letter alphabet with gaps.
Record identifiers have the form `t<alignment>_s<row>`, and the capture groups records into alignments by the part before `_`.

## Layout of `py27/`

Every function call is recorded as `{"value": ...}` or, where the original code raised an exception, `{"error": {"type": ..., "message": ...}}`.
Exceptions are part of the reference.
For example, the entropy consensus raises `KeyError` on `synth-mixedcase.aln.fasta`, where a residue type occurs only in lowercase, and on the single-sequence alignment `sequence_weights` divides by zero, as does every call that computes weights from it.
Floats are written with Python 2.7's `repr`, the shortest decimal string that reads back as the identical double, so a JSON reader recovers the computed values exactly.
Integers and floats remain distinct: a count of `3` and a weighted count of `3.0` are different values.

`library/<input>.json`, one per aligned input, holds the input path (`input`), the record identifiers in alignment order (`ids`), and one section per function:

- `sequence_weights`: one case per scaling (`none`, `sum1`, `max1`, `avg1`, `andy`).
- `aa_counts`, `aa_frequencies`: every combination of a weights specification and `gap_chars` of `-.` or `-.X`.
- `consensus`: the entropy algorithm for three weights specifications, and the simple algorithm, each with `trim_ends` true and false and `gap_threshold` 0.5, 0.6 and 0.8.
  These cover the calls made by fammer (`trim_ends=False, gap_threshold=0.6`) and cladecompare (`gap_threshold=0.8` with its sequence weights).
- `blocks`: `threshold` 0.4, 0.5 and 0.9, with `weights` `None` (sequence weights computed internally) and `False` (unweighted); records as `[id, description, sequence]`.
- `remove_empty_cols`: records as `[id, description, sequence]`.
- `score_pairwise`, `identity_abs`, `identity_avg`: every pair of rows, as `[i, j, outcome]` with `i < j`.

In every section except `remove_empty_cols` and the pairwise ones, each case gives the keyword arguments of the call in `kwargs`, so that it can be replayed as `func(aln, **kwargs)` once the weights are supplied.
The `weights` field of an `aa_counts`, `aa_frequencies` or `consensus` case names the sequence weights passed to the function:

| Name | Value |
|---|---|
| `null` | `None`, the default |
| `true` | `True`; `aa_counts` then computes `sequence_weights(aln)` itself |
| `"seqweights"` | `sequence_weights(aln, 'none')`, as cladecompare passes |
| `"unit"` | `[1] * len(aln)`, cladecompare's unweighted mode |

`sequence_weights(aln)` defaults to the `none` scaling, so `true` and `"seqweights"` produce the same weights, and their recorded results agree; they are two routes to one value, not independent evidence.
Where computing the `"seqweights"` weights itself fails, as on the single-sequence alignment, the error is recorded as the outcome of the case even though the function under test was never called.

`library/tie-batch.json` records, for each small alignment, the sequence weights under the scaling in `sequence_weights_kwargs` and the consensus under each combination of `simple` and `trim_ends` in `consensus_kwargs`.

`library/chisqprob.json` names the function in `function` and records `chisqprob(x, df)` in `rows` as `[x, df, outcome]`, over a grid of `x` from -1 to 1000 and `df` from 0 to 100, including two non-integer values of `df`.
The test for `x <= 0` precedes the test for `df <= 0`, so `x = 0` returns `1.0` for every `df`, including `df = 0`.

`scripts/<case>.json` records one command-line run: the script, its arguments, the input file given on standard input (if any), the exit status, standard output, standard error, and the contents of any files the script wrote.
A case is named `<script>.<input>[.<options>]`, with options abbreviated as on the command line; for example, `seq-blocks.cdc2.w.t09` runs `seq-blocks -w -t 0.9 cdc2-domain.aln.fasta`.
A final `.stdin` labels a case whose input arrives on standard input, and `.default` one with no options; the labels are not exhaustive, and the `argv` and `stdin` fields of each file record exactly how the script was run.
Each run takes place in a scratch directory holding copies of `inputs/`, so output names input files without a directory.
Standard error has the scratch directory, the repository root, the biocma extraction directory and the Python installation prefix replaced by `<workdir>`, `<repo>`, `<biocma>` and `<python-prefix>`.

## Interpreting a mismatch

- **Consensus ties.**
  Where several residues tie, the original code takes whichever comes first in Python 2's dictionary iteration order, which follows string hashes rather than any documented rule.
  The port replaces this with an explicit tie-breaking rule, so consensus strings may differ where ties occur; `tie-batch.json` exists to measure that difference.
  The same dictionary-order dependence affects the ordering of tied clustering coefficients in `seq-pick -c` and the medoids that `seq-pick -k` derives from them.
- **Summation order of weighted frequencies.**
  `aa_frequencies` normalizes by `sum(counts.values())`, a floating-point sum taken in dictionary order.
  With non-integer weights (`true` and `"seqweights"`), summing the same counts in another order changes the normalizer in its last bits, so a port that iterates in insertion order differs from these fixtures by about one unit in the last place in many frequencies.
  For example, the `"seqweights"` normalizer for `cdc2-domain.aln.fasta` is `1451.1196865831323` in hash order and `1451.1196865831325` summed in sorted key order.
  Compare these values with a relative tolerance on the order of `1e-15`, not with equality; with integer weights the sums are exact and equality holds.
- **Exception messages.**
  The exception type is the reference; the message text is the interpreter's.
  For example, Python 2 reports integer division by zero as `integer division or modulo by zero` and Python 3 as `division by zero`.
- **Log line numbers.**
  The `seq-intersect` cases log through `cmdutils.log_config`, whose format includes the source line of the logging call (`sequtils [@71]`), so editing `sequtils.py` above that line changes standard error without changing behavior.
- **Recorded defects.**
  The fixtures record the original behavior, defects included.
  For example, `seq-sort -l` is accepted but ignored, so `seq-sort.l` and `seq-sort.default` are identical, and `seq-diff -s` stops with a `ValueError` traceback because `cdc2-domain.fasta` contains one sequence twice.
  A port that fixes such a defect differs from the fixture by design, and the test should state the expected difference.
- **Biopython output.**
  `seq-convert`, `seq-extract`, `seq-purge`, `seq-sort`, `seq-intersect` and `seq-diff` write records through Biopython's writers, and `seq-dupes` through `SeqRecord.format`, so their output also depends on the Biopython version.
  The BiopythonWarning about partial codons in the `seq-orfs` cases is another example.
  The CMA output of `seq-purge` and `seq-diff` is written by biocma.
- **Not captured.**
  `seq-pick -r` and `stats.uniqueof20` sample at random.
  `logoutils` requires WebLogo, whose computations are its own rather than biofrills'.
  The on-disk index cache of `seq-intersect -c` is a side effect rather than an output.
