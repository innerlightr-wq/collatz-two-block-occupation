# Least Realizers and Occupation in an Aligned Two-Block Collatz Family

**Elias De Jesús** (ORCID 0009-0007-0190-9143)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24384828.svg)](https://doi.org/10.5281/zenodo.24384828)

Technical note and its full verification package.
**Paper: [`manuscript/occupation-note.pdf`](manuscript/occupation-note.pdf)** (13 pages),
archived at **[doi:10.5281/zenodo.24384828](https://doi.org/10.5281/zenodo.24384828)**.

## Citation

> De Jesús, Elias (2026). *Least Realizers and Occupation in an Aligned
> Two-Block Collatz Family.* Zenodo. https://doi.org/10.5281/zenodo.24384828

The concept DOI [10.5281/zenodo.24384827](https://doi.org/10.5281/zenodo.24384827)
always resolves to the latest version; the version DOI above pins the
deposited file.

The archived PDF is the reviewed release. The copy in this repository is the
same text with two additions to Appendix A: the link to this repository, and
runtimes replaced by measurements from a fresh clone. No mathematical content
differs.

Among the 64 residue classes mod 64, exactly two carry the aligned six-step
parity blocks `A = 111100` and `B = 111010`, namely 15 and 7, and on those
classes six steps of the Collatz shortcut map act as `F_A(x) = (81x+65)/64` and
`F_B(x) = (81x+73)/64`. The note asks whether one fixed positive integer can
follow these two blocks forever.

**The open problem is fixed-seed exit:** no positive integer is known to be
excluded from following `A` and `B` indefinitely, and proving
`D_m(0) ≤ 84m` for arbitrarily large `m` — or merely an unbounded deficit
below `5m·log_{81/64}64 = 88.274237898…m` — would settle it. Nothing in this
repository proves that, and nothing here bears on the Collatz conjecture
itself.

## Quick check (30 seconds, low memory)

```bash
python3 verification/scripts/rev2_checks.py
```

Expected final line:

```
47/47 checks passed
```

## Full verification

All scripts are standard-library Python 3, exact integer or rational
arithmetic, no network access, deterministic output. Each exits 0 on success
and nonzero on any failed check. Reference logs for every run are in
`verification/logs/`; timings will differ from machine to machine, the
mathematical output should not. All timings quoted here were measured on one
idle machine from a fresh clone of this repository.

### Light runs (well under 1 GiB each)

```bash
cd verification/scripts
python3 rev2_checks.py        #     47 checks    1 s   claims new in revisions 2-3
python3 note_checks.py 6      #  37912 checks    1 s   every numerical claim in the note
python3 windowcount.py 40     #    110 checks    7 s   B_t census, valuation ceiling
python3 prefix6m.py 6         #  34314 checks    1 s   D_m(0), mean-one lemma, comparison family
python3 averaging.py 40       #   1046 checks   10 s   Theorem 4.6 ingredients
python3 seedcount.py 6        #   4204 checks   17 s   Propositions 5.1-5.2 and Table 3
python3 m7.py                 #   prints D_7(0) = 2     6 s
python3 depth48.py            #     38 checks   <1 s   depth-48 audit and the density reduction
```

### Heavy runs — up to **1.3 GiB resident each**

Each stores 2^24 multi-precision integers; the largest measured peak is
**1.29 GiB**. **Run them one at a time**, not in parallel:

```bash
cd verification/scripts
for s in min49.py verify49.py note_verifier.py independent48.py; do
    echo "=== $s"; python3 "$s"; echo "exit=$?"
done
```

| script | what it establishes | time |
|---|---|---|
| `min49.py` | the exact least 49-block realizer, by cyclic successor search over all 2^49 words | 105 s |
| `verify49.py` | the same minimality, independently: an exact cyclic-interval count below the claimed value, using the complementary split of the letters | 95 s |
| `note_verifier.py` | the depth-48 control: exactly one 49-block-free start below 2^241 | 34 s |
| `independent48.py` | the same count by a structurally different implementation (branching on the first six letters, full 288-bit residues, cyclic interval test) | 4 min |

### Witness versus minimality, for the depth-49 result

These are different claims and are checked separately.

* That the 75-digit integer of Theorem 3.1 **is** a realizer is a *witness*
  property: it is replayed block by block and, independently, through its 294
  individual shortcut steps with parities recorded before each step, both
  routes returning the same word and the same terminal state 31 mod 64. This
  replay is instant and lives in `note_checks.py` (§3).
* That it is the **least** such integer is an *exhaustive* claim, established
  by `min49.py` and re-established by `verify49.py` with a different algorithm
  and a different split. Neither is implied by the replay.

## Contents

```
manuscript/occupation-note.tex    LaTeX source (pdflatex/latexmk, no external deps)
manuscript/occupation-note.pdf    the 13-page note
verification/scripts/             every verifier named in Appendix A of the note
verification/logs/                a reference log per script, each ending EXIT=0
CITATION_CHECKLIST.md             bibliography provenance and outstanding items
```

Rebuild the paper with:

```bash
cd manuscript && latexmk -pdf occupation-note.tex
```

## What the note contains

* an **occupation theorem**: if `b` positive integers below `H` complete `t`
  blocks, then `r_min(t+b) ≥ (H + 73/17)ρ^(-b) - 73/17`, with `ρ = 81/64`;
* the **exact least realizer at depth 49**, a 75-digit integer whose terminal
  state is 31 mod 64, so it completes exactly 49 blocks and then exits;
* the depth-48 control: every positive integer below `2^241 ≈ 3.53×10^72`
  exits before completing block 49;
* a **depth-averaging theorem**: for every fixed `ℓ`, `D_m(0) ≤ 2^(5m-ℓ)`
  outside a set of upper density at most `3·2^(-2ℓ)`, hence for every `ε > 0`
  the set `{m : D_m(0) > ε·2^(5m)}` has density zero — a vanishing *fraction*,
  which still permits an absolute count far above `84m`;
* a half-word reduction expressing the same count as a correlation between two
  explicit lists, with the bin bound `D_m(0) ≤ 2 + Σ e_j`.

## License

The MIT license in this repository covers the verification code. The technical
note itself is archived under **CC BY 4.0**, as recorded in the Zenodo deposit
[10.5281/zenodo.24384828](https://doi.org/10.5281/zenodo.24384828).
