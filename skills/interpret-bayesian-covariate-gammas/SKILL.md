---
name: "interpret-bayesian-covariate-gammas"
description: "Use when interpreting posterior estimates of a per-group covariate coefficient (e.g., a \"gamma\" effect fit separately per species, site, or other subgroup) from a Bayesian hierarchical model — especially when there are multiple covariate variants/models to compare and the user wants a structured, convergence-aware read of the results rather than just a table dump."
---

## Interpreting per-group Bayesian covariate effects ("gamma" results)

Use this structure whenever asked to interpret posterior summaries of a covariate
coefficient that was fit separately for many groups (species, sites, cohorts, etc.),
especially when several covariate variants exist and need comparing. The goal is a
narrative interpretation, not a table dump — name specific groups, don't just report
percentages.

Do the analysis with a script (Python/R/bash), not by eyeballing the CSV. Never modify
the user's files, models, or manuscript as part of "interpret" — that's a separate,
later step the user asks for explicitly.

### 0. Check what data you actually have

Before interpreting anything, check whether the file has just a point estimate (e.g.
posterior mean) or the full posterior summary (mean, median, sd, q5/q95, rhat,
ess_bulk, ess_tail). If only a point estimate is available:
- State this limitation explicitly and up front, not buried in a footnote.
- Report sign/direction patterns as directional hints only — never say "credible" or
  "significant" without an interval to back it up.
- If the user has a way to regenerate the file with full posterior summaries (e.g., a
  known upstream script that already computes q5/q95 for other parameters), say so
  and offer it as a next step.

If multiple files exist covering the same groups/models (e.g., a lightweight
per-run "diagnostics" file vs. a full posterior-summary file), check whether their
convergence columns measure the same thing. A "whole-model" max-rhat/min-ess (worst
value across every parameter in that group's fit — hundreds of nuisance parameters)
is a different, usually much scarier, number than a diagnostic computed for the one
parameter you care about. Reconcile any such discrepancy explicitly instead of
letting the scarier number silently override the more relevant one, and explain what
each file's diagnostic actually represents.

### 1. Convergence and estimation quality, first, before any values

Before interpreting a single coefficient, check:
- **Rhat**: flag > 1.01 as mild, > 1.05 as a real problem. Name the specific
  group/model combinations that are flagged, sorted worst-first — don't just give a
  count.
- **Effective sample size** (bulk and tail, or whatever the pipeline reports):
  compare against whatever threshold the project already uses (e.g., stated in a
  Methods section or fitting script) and flag anything below it.
- **Mean vs. median**: report the largest `|mean - median|` gap across all rows, and
  the typical gap. If it's small everywhere, say plainly that mean and median are
  interchangeable here and move on — don't hedge unnecessarily. If large for
  specific groups, flag those as likely skewed/unstable posteriors.
- Any group with a genuinely bad convergence flag (e.g., rhat > 1.1, ESS in the
  tens) should be named explicitly and described as unreliable — recommend excluding
  it from reported results or refitting it, not just noting the number and moving on.

### 2. Go through each covariate model one at a time

Never lump multiple covariate variants into one summary. For each one, report in its
own paragraph:
- N groups with an estimate.
- Sign split on point estimates (n positive / n negative).
- If CIs are available: N credible (interval excludes zero), split into credibly
  positive vs. credibly negative, out of the total tested.
- Mean and median of the point estimates, and the range.
- **Name the specific groups** with credible (or, if no CI, the largest-magnitude)
  effects, sorted by estimate — this is what lets patterns become visible. A bare
  percentage hides the interesting part.

### 3. Cross-covariate consistency

If two or more covariates conceptually overlap (e.g., a "current level of X" and a
"rate of change into X" covariate), explicitly:
- Correlate their point estimates across groups.
- Count how many groups agree vs. disagree in sign.
- Flag any group where two covariates are *both* credible but point in opposite
  directions. Don't treat this as an error — explain the plausible interpretation
  (e.g., a static level and its recent rate of change are genuinely different
  quantities and can diverge for a specific group).

### 4. Call out unusual groups

For each covariate model, explicitly name:
- Any credible group whose direction is opposite the dominant direction in that
  model. Offer a plausible domain explanation if one is inferable from what's known
  about that group (e.g., a species known to prefer human-modified habitat showing
  up as the exception to a "natural habitat is good" pattern).
- Any group with unusually large effect magnitude relative to the rest.

### 5. Synthesize, then give a clear bottom line

Close with:
- Which covariate(s) have the strongest, most one-directional, most credible signal
  vs. which are weak or inconsistent — this is usually the single most decision-useful
  sentence for the user.
- If comparing across multiple runs/populations (e.g., the same covariate structure
  fit to two different species groups), explicitly compare the patterns rather than
  just reporting each in isolation — note where they agree and where they diverge.
- Any remaining caveats (missing CIs, unresolved convergence issues, small sample
  sizes for particular groups) and a concrete next step to resolve them if relevant.

### 6. Save the write-up

Once the interpretation has been given in chat, save it as a plain-text file named
after the group/dataset being analyzed (e.g. `<group>.txt` — species group,
site name, cohort, whatever the natural grouping label is), placed alongside the
source CSVs it was derived from. Use plain-text section headers (no markdown bold
markers), and keep the same structure as the chat response: mean-vs-median,
convergence, one section per covariate model, cross-covariate nuance, and a bottom
line. This makes the interpretation a durable project artifact rather than
something that only exists in the conversation. If a comparison to a
previously-interpreted group exists (e.g. a prior `<other_group>.txt`), include
that comparison in the bottom line the same way it was discussed in chat.

### Style

Write in flowing paragraphs, one per covariate model plus a synthesis paragraph —
not nested bullet lists. Keep numbers precise and inline (means, credible counts,
specific CI bounds for called-out groups). Don't editorialize about what the user
should conclude scientifically — describe the pattern and let them decide what it
means for their work, except where a plausible domain explanation directly explains
an outlier.
