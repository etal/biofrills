# BioFrills

Two filenames reach this one file: OMP and Codex read it as `AGENTS.md`, while Claude Code reads it as `CLAUDE.md`, a symlink recorded in git.
Where symlinks do not survive checkout (Git for Windows outside Developer Mode, or a DrvFS mount such as `/mnt/c`), `CLAUDE.md` arrives as a stub naming this file; read `AGENTS.md` directly rather than pasting content into the stub.

Utilities for biomolecular sequence analysis, built on Biopython: a library of alignment and sequence-set functions, and 13 `seq-*` command-line programs.
The 2013 Python 2 code base has been ported to Python 3.12+ and current Biopython; `README.rst` is the user and developer documentation.

## Consumers

BioFrills exists mainly to support [CladeCompare](https://github.com/etal/cladecompare) and [Fammer](https://github.com/etal/fammer), which pin it to a git commit.
The names they call are the compatibility contract: `alnutils.sequence_weights`, `alnutils.aa_frequencies`, `alnutils.blocks`, `alnutils.remove_empty_cols`, `consensus.consensus` (with `weights`, `gap_threshold` and `trim_ends`) and `logoutils.aln2logodata`.
A difference in the last place can decide a gap column whose weighted gap fraction lies exactly at `gap_threshold`.
A gap decision in `consensus` becomes a match column of a Fammer MAPGAPS profile, or a column that CladeCompare assigns p = 1 instead of testing, and a change in `aa_frequencies` moves CladeCompare's pseudocounts and with them its p-values.
Treat any change to these functions as a change of output, and check it against the reference data below.

## Layout

- `biofrills/`: the library modules; `biofrills/cli/` holds one module per console script, registered in `[project.scripts]`, and in `__init__.py` the argument helpers they share.
- `tests/`: the pytest suite; `tests/alignments.py` holds shared helpers, importable because `pythonpath = ["tests"]`.
- `tests/data/`: the Python 2.7 reference outputs and the tools around them, documented in `tests/data/README.md`.

## Commands

```sh
pixi run test              # newest Python in the lock file
pixi run -e py312 test     # Python 3.12, the oldest supported version
pixi run lint              # ruff check
pixi run fmt               # ruff format
uv run pytest              # the same suite with the PyPI dependencies of uv.lock
```

The command-line programs are tested as installed, so after a change to `[project.scripts]` run `pixi install` or `uv sync` before the tests.

## Reference data

- `tests/data/py27/` is captured output of the pre-port sources; never edit it by hand.
  The capture runs from a git worktree of the pre-port commit `615d27c` beside this repository, never from the working tree, in a Python 2.7 environment such as the osx-64 `py27-oracle` pixi environment; `tests/data/README.md` gives the steps.
- `tests/data/capture_py27.py` is frozen Python 2 code; ruff excludes it, and it must not be modernized.
- `tests/reference_cases.py` enumerates every recorded library case of a function that still exists, and both `tests/test_reference.py` and `tests/data/replay.py` use it.
  Floating-point values are compared with the relative tolerances defined there, each of which `tests/data/README.md` justifies; do not widen one to make a test pass.
- A change that deliberately moves a recorded library outcome is recorded by regenerating `tests/data/differences.json` with `pixi run python tests/data/replay.py`, and each new entry is explained under "Differences in library output" in `tests/data/README.md`.
  `git diff tests/data/differences.json` then shows exactly which outcomes moved.
- `tests/test_cli.py` runs each installed console script as a subprocess, with the default warning filters and `PYTHONHASHSEED=0`, as the recorded runs had; it sets `COLUMNS=80`, the width argparse assumed when those runs wrote to a pipe.
  It normalizes the log line numbers (`[@N]`) that the log format prints and reduces a traceback to its first and last lines.
  A change that deliberately moves a recorded command-line run is listed in `EXPECTED` there, as a function deriving the expected run from the recorded one, and explained in `tests/data/README.md` under "Command-line differences" or under the section that names its cause, such as "Consensus ties".
  Its CMA cases are skipped unless biocma can be imported; with a checkout of biocma beside this repository, `PYTHONPATH=$PWD/../biocma uv run pytest tests/test_cli.py` runs them, so run it after any change to the CMA paths of `seq-diff` or `seq-purge`.

## Numerical behavior

- `consensus` decides weighted gap columns, and `alnutils.blocks` weighted non-gap columns, with `math.fsum` and a `1e-12` margin below the threshold (`alnutils._GAP_TOLERANCE`), which makes the decision independent of row order and of the rounding of the sequence weights.
  The weighted counts of `aa_counts` still accumulate in row order, the `sum1`, `avg1` and `andy` scalings of `sequence_weights` still divide by a built-in `sum`, and the `aa_frequencies` normalizer still sums in dictionary order; changing any of these moves recorded output and needs its own decision.
- Ties between residue types in `consensus` are broken by an explicit rule ending in alphabetical order; never let dictionary or set iteration order decide an output.

## Letter case

- `blocks`, `remove_empty_cols` and `consensus` follow the alignment convention of HMMER and CHAIN: an uppercase residue is aligned, a lowercase residue is an insertion, `-` is a deletion and `.` a gap in an insert column.
  The rule is defined once, in `alnutils._GAP_CHARS`, `alnutils._case_marks_inserts` and `alnutils._has_aligned_residue`; change it there rather than in one of its callers.
- An alignment without any uppercase residue marks nothing by case and is read without regard to it; MAFFT writes nucleotide alignments in lowercase, so this is a common input rather than a corner case.
- `blocks` gives a column whose residues are all lowercase a proportion of zero, and `remove_empty_cols` keeps it, because that column holds residues of its sequences and the callers of `remove_empty_cols` write its result out as an alignment: Fammer its group profile alignments and its structural seed alignments, and CladeCompare a combined alignment that it splits back into its two parts.
  Fammer passes it only protein alignments, from MAFFT, which writes protein residues in uppercase, and from TM-align, whose recorded output holds no lowercase residue, so in Fammer the rule changes no output.
  The asymmetry is deliberate.
- `aa_counts`, `aa_frequencies` and `sequence_weights` count residues as written, so `a` and `A` are distinct types there, whereas the background frequencies of `consensus` fold case.
  Changing that would move CladeCompare's pseudocounts and p-values and needs its own decision.
- `pairutils.score_pairwise`, `identity_abs` and `identity_avg` compare residues without regard to case, and `seq-identity` and `seq-pick -c` and `-k` inherit that.
  They do not apply the insert convention: a function of two rows cannot tell whether case marks inserts in the alignment they came from, and `seq-pick -k` charges a medoid `pairutils.ungap_len`, its number of residues with lowercase included, which must equal the identity of the sequence with itself.
- Lowercase does not always mark an insertion: BLAST and PSI-BLAST use it for masking when `-lcase_masking` is given, which the convention would misread.

## Pitfalls

- pytest turns warnings into errors.
- ruff's automatic fixes are not behavior-neutral on this code base; the unsafe B006 fix to `make_simple_col_consensus` once raised the number of differing recorded `consensus` outcomes by 284.
  Apply such fixes by hand and replay the reference data.
- Library modules that log do so through `logging.getLogger(__name__)`; do not add a `NullHandler`, which would hide warnings from callers that configure no logging.
- networkx (`graph` extra), weblogo (`logo` extra) and biocma (CMA input) are optional: `alnutils.to_graph` and `seq-pick` import networkx inside the function, and `seq-diff` and `seq-purge` import biocma through `biofrills.cli.import_cma`, which exits with a message when CMA input is given without it.
  Keep these imports out of module level, except in `logoutils`, which exists to wrap weblogo.
- Runtime dependency floors default to the versions that Ubuntu 26.04 ships; raise one only for a concrete problem, such as a security advisory or a conflict.
  After a change to a dependency, update both `uv.lock` and `pixi.lock`.

## Conventions

- Format with `ruff format`, and list each formatting-only commit in `.git-blame-ignore-revs`.
- Keep modernization and behavior changes in separate commits, so that a numerical difference never hides in a formatting-only commit.
- Commit messages: an imperative subject line, then a short paragraph on what changed and why; American English.
- In Markdown and reStructuredText files, put each sentence on its own line.
