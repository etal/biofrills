# Reference outputs from the Python 2 code base

The files under `py27/` record what biofrills computed on Python 2.7, before the port to Python 3.
They are the reference against which the port is checked: a ported function should reproduce these values, or the difference should be explained in this file.

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
The current code returns values in both cases; see "Differences in library output" below.
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

In the sections that call a function on the whole alignment, each case gives the keyword arguments of the call in `kwargs`, so that it can be replayed as `func(aln, **kwargs)` once the weights are supplied.
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

`tests/test_cli.py` replays every script case against the installed console scripts in the same way: as a separate process in a scratch copy of `inputs/`, with `PYTHONHASHSEED` set to `0`.
It compares the exit status, standard output and written files exactly, and standard error after the normalizations described below under "Log line numbers" and "Recorded defects".
After those normalizations no recorded line contains one of the placeholders above, so the replay substitutes no paths.
The cases that read CMA files run only where biocma can be imported; with a Python 3 checkout of biocma beside this repository, `PYTHONPATH=$PWD/../biocma uv run pytest tests/test_cli.py` runs them.
The path must be absolute, because each command runs in a scratch directory.
The cases whose current output differs, explained below under "Consensus ties", "Ambiguous residues in `seq-purge -c`" and "Command-line differences", are listed in `EXPECTED` in `tests/test_cli.py`.

## Interpreting a mismatch

- **Consensus ties.**
  Where several residues tie, the original code takes whichever comes first in Python 2's dictionary iteration order, which follows string hashes rather than any documented rule.
  The current code breaks such ties by an explicit rule instead, so consensus strings may differ where ties occur; the cases that changed are listed under "Differences in library output" below, and `tie-batch.json` exists to measure that difference.
  The same order dependence affected `seq-pick -c`, which sorts records by clustering coefficient, and `seq-pick -k`, which starts from those records.
  In `cdc2-domain.aln.fasta`, the records `CDC2-1_Sm` and `CDC2-2_Sm` have the same aligned sequence, so their weighted clustering coefficients are equal in exact arithmetic.
  As computed by networkx they differ by up to one unit in the last place, depending on the order in which networkx iterates over sets, that is, on the hash seed.
  `seq-pick -c` therefore compares coefficients rounded to 12 decimal places and orders equal ones by record identifier.
  With hash seeds 0 to 7, the coefficient of each record of that alignment varies by at most `3.9e-16`, while the closest two distinct coefficients differ by `8.9e-6`, so the rounding removes the variation and merges no distinct values.
  `seq-pick.c5` then matches the fixture with each hash seed tried, from 0 to 15.
  In `seq-pick.k3` the k-medoids cost of the two records ties exactly, and the current code chooses `CDC2-2_Sm`, the first of the two in the file, where the fixture has `CDC2-1_Sm`; this choice does not depend on the hash seed.
- **Summation order of weighted frequencies.**
  `aa_frequencies` normalizes by `sum(counts.values())`, a floating-point sum taken in dictionary order, which is string-hash order on Python 2 and insertion order on Python 3.
  With non-integer weights (`true` and `"seqweights"`), summing the same counts in another order changes the normalizer in its last place, and every frequency with it.
  For example, the `"seqweights"` normalizer for `cdc2-domain.aln.fasta` is `1451.1196865831323` in hash order and `1451.1196865831325` in insertion order, which is also the correctly rounded sum.
  Under Python 3 this moves 80 weighted frequencies of `cdc2-domain.aln.json` and 32 of `synth-x.aln.json`, by at most 2 units in the last place; the weighted frequencies of the other inputs, and all counts, are unchanged.
  Compare these values with a relative tolerance on the order of `1e-15`, not with equality; with integer weights the sums are exact and equality holds.
- **Compensated summation.**
  Since Python 3.12, the built-in `sum` of floating-point numbers compensates for rounding error (Neumaier's algorithm), where Python 2 added the terms one at a time.
  `sequence_weights` divides by `sum(seq_weights)` under the scalings `sum1`, `avg1` and `andy`, so these weights differ from the fixtures by at most 3 units in the last place: 35 values in `scttlre-domain.aln.json`, 8 in `synth-x.aln.json`, and 185 in 39 of the 200 alignments of `tie-batch.json`.
  The largest relative difference is `4.2e-16`, so the tolerance given above for weighted frequencies covers these values too.
  The scalings `none` and `max1` take no such sum and match the fixtures exactly.
  `consensus` divides by a sum of these weights to compute the weighted gap fraction of a column, and there a difference in the last place can decide a column whose exact fraction equals the threshold.
  Column 5 of `tie-batch.json` alignment `t101` holds two gaps, a `C` and an `A`, and in exact arithmetic the `avg1` weights of its four rows are 20/21, 22/21, 20/21 and 22/21, so the weighted gap fraction of the column is exactly one half.
  As stored, the two weights of 20/21 differ by one unit in the last place: Python 2's term-by-term total, `3.999999999999999`, gives a fraction of exactly `0.5`, and the correctly rounded total, `3.9999999999999996`, gives `0.49999999999999994`.
  `consensus` therefore sums the weights with `math.fsum`, which is correctly rounded and so does not depend on the order of the rows, and counts a weighted gap fraction at most `1e-12` below `gap_threshold` as reaching it; the column is a gap, as recorded, so the case produces no entry in `differences.json`.
  The measurements behind that margin cover the fixture inputs and 2,100 random alignments: 600 of 4 to 40 rows and 5 to 60 columns, and 1,500 of 3 to 8 rows and 10 to 60 columns, with gap probabilities from 0.1 to 0.6.
  There, the largest rounding error of the fraction was `5e-16`, as it was in synthetic alignments of up to 1,000 rows or 20,000 columns, and the smallest distance from the threshold of a fraction not exactly equal to it was `9.5e-6`.
  The margin is thus about 2,000 times the first and about seven orders of magnitude below the second.
- **Platform math library.**
  `chisqprob` calls `math.lgamma`, whose result depends on the platform's C math library.
  The fixtures were captured on x86_64 macOS under Rosetta 2, where Python 3.12 reproduces every `chisqprob` value exactly.
  On arm64 macOS, `math.lgamma(5.0)` is `3.1780538303479444` rather than `3.178053830347945`, and 10 of the 273 recorded values, all with `df` 2.5, 10 or 100, differ from the fixture by up to 246 units in the last place, a relative difference of at most `2.9e-14`.
  These values need a relative tolerance of about `1e-13`, larger than the one given above for summation.
- **Exception messages.**
  The exception type is the reference; the message text is the interpreter's.
  For example, Python 2 reports integer division by zero as `integer division or modulo by zero` and Python 3 as `division by zero`.
- **Log line numbers.**
  The `seq-intersect` cases log through `cmdutils.log_config`, whose format includes the source line of the logging call (`sequtils [@71]` in the fixtures), so editing `sequtils.py` above that call changes standard error without changing behavior.
  `tests/test_cli.py` replaces that number with `N` on both sides.
- **Recorded defects.**
  The fixtures record the original behavior, defects included.
  For example, `seq-sort -s` was documented to sort by similarity but sorted by increasing length, and `seq-diff -s` stops with a `ValueError` traceback because `cdc2-domain.fasta` contains one sequence twice.
  The frames of that traceback name files and lines of the interpreter, of Biopython and of the script, so `tests/test_cli.py` compares only its first and last lines.
  A port that fixes such a defect differs from the fixture by design, and the test should state the expected difference.
  For the library cases, the expected differences are recorded in `differences.json`, described under "Differences in library output" below.
- **Ambiguous residues in `seq-purge -c`.**
  The recorded rule split a sequence at each `X` and purged the record when every segment occurred in some kept sequence, in any order and possibly in different kept sequences; a leading, trailing or doubled `X` yields an empty segment, which occurs in every sequence.
  The current rule, which fammer's `update_fasta.ungap_and_unique` also applies, purges a record containing `X` only where the whole sequence, with each `X` matching any single character, occurs within one kept sequence; an `X` in a kept sequence matches only an `X`.
  In `seq-purge.c` the record `x_fragment_of_1` is therefore kept, because its two segments both occur in the first record of `dupes.fasta` but 20 residues apart.
  Standard error reads `Purged 5 sequences` instead of `Purged 6 sequences`, and standard output gains that record between `psu|PFF0750w` and `x_novel`.
  The other `seq-purge` cases are unchanged, including `seq-purge.fikk-cma.c`, where both rules keep the one record that contains `X`.
  `replay.py` covers only the library cases, so this difference is recorded here and in `tests/test_cli.py` rather than in `differences.json`.
- **Command-line differences.**
  Five script cases differ from their fixtures because the port corrects the command-line programs.
  In `seq-orfs.default` and `seq-orfs.m20`, standard error no longer holds Biopython's warning about a partial codon, because `seq-orfs` cuts each reading frame to a whole number of codons before translating it.
  Biopython already ignored the trailing partial codon, so standard output is unchanged.
  In `seq-identity.cdc2` and `seq-identity.edge`, the rows of the grid no longer end in a tab; the values are unchanged.
  `seq-sort.s` now exits with status 2 and a usage error, because `seq-sort -s` has been removed rather than implemented.
  It sorted by increasing length, not by similarity, and every similarity measure in biofrills requires aligned input, while `seq-sort` sorts unaligned sequences; `seq-pick -c` selects the most central sequences of an alignment.
  `seq-sort -l` now selects the default order explicitly, so `seq-sort.l` still matches its fixture.
- **Biopython output.**
  `seq-convert`, `seq-extract`, `seq-purge`, `seq-sort`, `seq-intersect` and `seq-diff` write records through Biopython's writers, and `seq-dupes` through `SeqRecord.format`, so their output also depends on the Biopython version.
  The CMA output of `seq-purge` and `seq-diff` is written by biocma.
- **Not captured.**
  `seq-pick -r` and `stats.uniqueof20` sample at random.
  `logoutils` requires WebLogo, whose computations are its own rather than biofrills'.
  The on-disk index cache of `seq-intersect -c` is a side effect rather than an output.

## Differences in library output

`replay.py` calls every recorded library case again with the current sources of this repository and writes each case whose outcome differs from the recorded one to `differences.json`.
It compares outcomes as `tests/test_reference.py` does, through the shared module `tests/reference_cases.py`, so the floating-point differences described above under "Summation order of weighted frequencies", "Compensated summation" and "Platform math library" are compared with a tolerance and produce no entries.
An entry identifies the case as the fixture does: the fixture file relative to `py27/`, the section, the alignment for a `tie-batch.json` case, and the index of the case within its section, which is `null` for a section that records a single outcome.
It then repeats the fields that describe the call (`weights` and `kwargs` where the fixture records them, `rows` for a pairwise case, `x` and `df` for a `chisqprob` case), and gives the recorded outcome, the current outcome, and, where both are strings of equal length, the 1-based columns at which they differ.
Every entry is explained below, and all four kinds are deliberate changes of behavior.
In the third and fourth kinds the recorded outcome is an exception, so they change no outcome that the original code completed.
A test of the current code expects the fixture value except where `differences.json` records another.
After a change that moves a recorded output, regenerate the file from the repository root and explain the new entries here:

```sh
pixi run python tests/data/replay.py
```

Before the port, the unmodified sources at commit `36684bb` produced an empty list under Python 2.7.

- **Column order with `trim_ends=False`** (10 cases, all in `synth-edge.aln.json`).
  A column with no residue to report, either all gaps or all lowercase insert residues, yielded its gap immediately, ahead of any C-terminal residues that `consensus` was still holding back, so the string left column order.
  In `synth-edge.aln.json`, column 14 holds a single `Q` and column 15 is all gaps; the recorded `...YRGH-Q` becomes `...YRGHQ-`, which keeps each character at its own column, as callers that pair the consensus with alignment columns assume.
  The other `trim_ends=False` cases of that file already ended in `Q-`, because the `Q` column was not a majority-gap column at their threshold.
  Cases with `trim_ends=True` are unchanged, because consensus holds nothing back in that mode.
  The three MAPGAPS alignments are unchanged as well: in them no residue-less column follows held-back C-terminal residues.
- **Tie-breaking rule** (23 cases).
  The entropy algorithm chooses the residue with the largest relative-entropy term and, among residues whose terms are exactly equal, the alphabetically first.
  The simple algorithm keeps its two heuristics, preferring the residue that most often follows the previous consensus residue and then the residue least frequent in the whole alignment, and adds alphabetical order as the final key.
  Both rules order tied residues totally, so the choice among them no longer depends on dictionary iteration order: under Python 2.7 with `PYTHONHASHSEED` set to 1, 2 and 3, the previous code changed 5, 86 and 5 recorded cases relative to the default ordering, and the current code changes none.
  The weighted background frequencies of the entropy algorithm are still normalized by a sum taken in dictionary order (see "Summation order of weighted frequencies" above); that can move a term by about one unit in the last place, but no recorded case has two distinct terms that close.
  Every exact tie of the entropy term in the fixtures is between residues with equal column and background frequencies, so no intermediate key, such as background frequency, could separate them.
  No case of the three MAPGAPS alignments changes; the 23 that do are these:
  - `synth-ties.aln.json`, all 18 entropy cases: column 8 holds `W`, `W`, `Y`, `Y`, and the consensus there changes from `Y` to `W`.
  - `tie-batch.json` alignment `t045`, both entropy cases: `D` and `E` tie in column 8, which changes from `E` to `D`.
  - `tie-batch.json` alignment `t157`, both simple cases: four residues tie in column 6, which changes from `E` to `D`; column 7 then changes from `D` to `A`, because the first heuristic follows the new residue in column 6.
  - `tie-batch.json` alignment `t081`, simple algorithm with `trim_ends=False`: three residues tie in column 4, which changes from `E` to `D`.
    With `trim_ends=True` the majority-gap column 3 does not reach the algorithm, so the first heuristic conditions on a different preceding column and settles the tie as before.
- **Sequence weights of an alignment without informative columns** (28 cases, all in `synth-single.aln.json`).
  `sequence_weights` takes its weights from the columns in which at least two residues, and at least a fifth of the rows, are not gaps, and a single sequence has no such column, so every raw weight was zero and every scaling divided by zero, as did every call that computes weights from it.
  Where nothing in the alignment tells the sequences apart, the current code weights them equally: `1.0` each under every scaling except `sum1`, which gives `1/N`.
  Up to rounding, these are the weights that every scaling already gives an alignment of identical sequences without gaps, and under the `none` scaling they sum to the number of sequences, the effective number of independent sequences when no redundancy can be measured.
  The 5 `sequence_weights` cases now give `[1.0]`.
  The other 23 cases now give the outcomes recorded for the same call with unit weights or without weights: the 8 weighted `aa_counts` and `aa_frequencies` cases match the `"unit"` cases, the 12 entropy `consensus` cases with weights `null` and `"seqweights"` give `MKV-LAYRGH`, the input sequence, and the 3 `blocks` cases with weights `null` give `MKVLAYRGH`, as with `false`.
- **Case-insensitive background frequencies** (18 cases, all in `synth-mixedcase.aln.json`).
  `consensus` uppercases a column that mixes lowercase and uppercase residues, but counted the background frequencies of the alignment as written, so an uppercased residue type that occurs elsewhere only in lowercase had no background frequency, and the entropy algorithm raised `KeyError`.
  The background frequencies of both algorithms now fold case as well; `alnutils.aa_frequencies` itself still counts residues as written.
  All 18 entropy cases of that file now complete: they give `MK--LWGH` with `gap_threshold` 0.5 and `MK-VLWGH` with 0.6 and 0.8, as the simple algorithm already did.
  No other input contains a lowercase residue, so no other recorded case changes; on an alignment that does, background frequencies, and with them a consensus residue, can differ from those of the original code.
