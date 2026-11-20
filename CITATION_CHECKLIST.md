# Citation checklist — revision 3 (released and archived)

**This note is archived: doi:10.5281/zenodo.24384828** (version),
concept DOI 10.5281/zenodo.24384827, deposited 20 November 2026 under
CC BY 4.0. Verified first-hand against the Zenodo API: the record's title,
sole creator ("De Jesus, Elias"), publication date, licence and the single file
`occupation-note.pdf` all match, and `https://doi.org/10.5281/zenodo.24384828`
resolves to `https://zenodo.org/records/24384828`.

**Deposited file versus this repository.** The archived PDF is textually
identical to the reviewed release (`two-block-occupation-note-rev3`). The copy
in this repository differs in Appendix A only: it adds the sentence linking to
the public verification package, and it replaces the estimated runtimes with
measurements taken from a fresh clone. The mathematical text, tables and
bibliography are unchanged. If you want the archive to carry the link and the
measured timings, upload the current PDF as a new version under the concept
DOI; that is a deposit decision, not made here.


Every entry below names the claim it supports and the status of its
bibliographic record. Nothing here blocks the build; items marked
**DOI to confirm — Elias** are the only ones needing your attention.

## Verified records (DOI printed in the bibliography)

| key | role in the paper | where | record |
|---|---|---|---|
| `Terras1976` | parity-vector cylinder background for Lemma 1.3 | §1.4 | Acta Arith. **30** (1976) 241–252, doi:10.4064/aa-30-3-241-252 |
| `Lagarias1985` | general survey background | §1.4 | Amer. Math. Monthly **92** (1985) 3–23, doi:10.2307/2322189 |
| `MonksYazinski2004` | the earlier lower-density inequality liminf(h/ℓ) ≥ ln2/ln3, **new citation in this revision** | §1.4, Problem 1 | Discrete Math. **275** (2004) 219–236, doi:10.1016/S0012-365X(03)00125-0 |
| `DeJesusSturmian2026` | §1.4: the Sturmian background — the irrational-seed statement and its independent proofs, and the quantitative/transcendence results | §1.4 | E. De Jesús, *Transcendence of the 3x+1 conjugacy map on Sturmian words*, version 2, 5 Nov 2026. **Verification status:** theorem numbers and section numbers checked **first-hand against the local Version 2 build** `~/transcendence-paper-v2/main.pdf` — Theorem 1.2 ("Quantitative form"), Corollary 1.4 ("Unconditional range", γ < γ* = 0.8258977…, which contains log₃2), §1.2 "Prior irrationality results", §9 "Scope and open problems", §10 "Related work". **Status of the identifiers, as of this revision:** the v2 posting date of 5 November 2026 is confirmed on Hexagon, and the DOI `10.5281/zenodo.24015456` is corroborated by an independent November preprint that cites it with that exact identifier; the Zenodo record itself has not been inspected directly. Neither URL could be retrieved in this environment (no network access). The DOI is printed in the bibliography; the "author-supplied" qualification lives here rather than in print. The Zenodo record for that paper is still uninspected. |
| `LopezStoll2009` | the open question itself: their abstract records that it is unknown whether any aperiodic parity vector has an eventually periodic (equivalently rational) Φ-value — now cited at that sentence in §1.3, verified first-hand from the paper | §1.3 | Integers **9** (2009) #A13, 141–162, doi:10.1515/INTEG.2009.014 |

## Preprint (no DOI claimed)

| key | role | status |
|---|---|---|
| `LopezStoll2021` | Theorem 1, the general density claim that the restricted consequence of this paper would follow from | arXiv:2101.12747**v1** [math.NT], 29 Jan 2021. **doi:10.48550/arXiv.2101.12747** (arXiv-issued via DataCite), now printed in the bibliography; the `v1` identifier is retained because that is the version used. No journal record located, no v2. The paper does not rely on its proof. |

## Accompanying computations (no bibliography entry)

The depth-averaging and half-word results are proved in this paper, so the two
working reports behind them are **not** cited: there is no separate manuscript
under either title, and neither entry was referenced in the body. The
calculations ship with the paper and are credited in Appendix A as
`averaging.py` (Theorem 4.6) and `seedcount.py` (Propositions 5.1–5.2 and
Table 3), each with its log. Nothing here needs a DOI.

## Checked and deliberately **not** cited

* *Finite-Window Survival Counts and Residue Refinement for the Odd Collatz
  Map* (E. De Jesús, 1 Nov 2026), Propositions 16 and 19 — inspected for an
  occupancy or refinement comparison. **No result of this revision uses it**,
  so it is not cited. (Unpublished 15-page draft; antecedent record
  doi:10.5281/zenodo.22255293 relates to the earlier window note, not to this
  draft. **DOI to confirm — Elias** if you decide to cite it after all.)
* Earlier conditional-lift and EOC notes — checked; the prefix lift of
  Theorem 4.6 is derived here from Lemma 1.3 and Lemma 4.3 of this paper, not
  imported, so no citation is owed.
* rStation Papers 103 and 104, and J. Allikvere, *Periodic shadows, drift
  pressure, and transcendence of the 3x+1 conjugacy map* — these bear on the
  Thue–Morse and Fourier investigations, which remain in their own reports and
  supply **no lemma or comparison to this manuscript**. If any of that material
  is later imported, cite the exact version: rStation Paper 104 (Lyapunov
  Weaver, 5 Nov 2026, Theorem 1) and Paper 103, both web-hosted at
  rstation.ai — **DOI to confirm — Elias**, none known to have been assigned;
  Allikvere, preprint, LaTeX build 13 Aug 2026 — **DOI to confirm — Elias**.

## Section numbering: local build versus posted version

The numbering was checked first-hand against the local Version 2 build
`~/transcendence-paper-v2/main.pdf`, in which the back matter is §9 "Scope and
open problems" and §10 "Related work". The **posted** v2 has those as **§8 and
§9**: the local build is one section ahead of the public file. This does not
affect the citation, which is to **§1.2 "Prior irrationality results"**, a
section present and identically numbered in both; the posted §1.2 credits four
independent routes to Sturmian irrationality (Pham, Allikvere, Cassidy and
De Jesús), where the local build lists three. Theorem 1.2 ("Quantitative
form") and Corollary 1.4 ("Unconditional range") agree in both.

## Correction made to the suggested wording

The §1.4 paragraph was supplied with "De Jesús (2026, §9)" for the tracing of
the independent irrationality proofs. **That section number is wrong.** In the
Version 2 build, §9 is "Scope and open problems" and §10 is "Related work";
the three strands (Monks–Yazinski Theorem 2.7(b) below β; López–Stoll §4 for
the *real* value Φ_ℝ, a different number; and De Jesús together with Cassidy,
independently and contemporaneously, Cassidy's note posted first) are traced in
**§1.2, "Prior irrationality results."** The paragraph as inserted cites
§1.2 and keeps the attribution collective.

## Errata applied in the final pass

* Removed the two self-referential bibliography entries and their checklist
  items; the scripts remain credited in Appendix A.
* Added the arXiv-issued DOI to the 2021 preprint, retaining `v1`.
* §3: "full 294-bit residues" corrected to "residues modulo $2^{288}$" — the
  294-bit modulus belongs to depth 49, and the same paragraph already states
  $64^{48}=2^{288}$.
* The 2009 Sturmian paper now supports a specific sentence in §1.3 and is cited
  there; verified first-hand from the paper's abstract.
* Rebuild check: all six bibliography entries are cited, and every in-text
  citation resolves (no uncited entries, no dangling keys).

## Disposition of the two formerly "unpublished report" entries

* `Averaging2026` ("Depth averaging in the two-block count") — **no distinct
  manuscript exists**; the title was a label for an internal calculation. Entry
  removed in revision 2; the result is Theorem 4.6 with Corollary 4.7 of this
  note, and the calculation is credited in Appendix A as `averaging.py`.
* `SeedCount2026` ("The affine-carry seed count for aligned two-block words")
  — **no distinct manuscript exists**, same disposition. The results are
  Propositions 5.1–5.2 and Table 3 of this note, calculation credited in
  Appendix A as `seedcount.py`.

## Rule applied

No DOI was inferred from a local filename or a title match. Where a record
could not be verified from the material at hand, the entry gives author,
title, date or version, and status, and omits the DOI field rather than
guessing.
