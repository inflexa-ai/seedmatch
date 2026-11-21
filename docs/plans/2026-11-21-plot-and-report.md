# Plots and the report HTML: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port each plot of SeedMatchR to `seedmatch.plot`, and the `gt` table of `SeedMatchReport` to `seedmatch.report_html`. The corrected rules of `D-063` to `D-070` apply. These are the parts:

- the ECDF plot of `deseq_fc_ecdf` and `.plot_ecdf`, with the theme and the palette of `R/themes.R`
- the plots `dens.plot`, `cdf.plot`, and `diff.plot` of `wass_dist`
- the histogram of `kmer_hist`
- the figure of `plotAlignment` with `output.type = "ggplot"`
- the range figures of `plotIRanges` and `plotRangesDF`, as one function
- the PR and ROC figures that `plot` of PRROC 1.4 draws for the result of `benchmark_preds`
- the `gt` table of `SeedMatchReport` (`.make_table`)

**Architecture:** The data of a figure comes from a core function. The figure comes from a function of `seedmatch.plot`. These are the new core functions:

- `seedmatch.stats.wasserstein.wasserstein_difference`: the grid and the ECDF difference of `diff.plot`.
- `seedmatch.stats.density`: the `bw.nrd0` bandwidth and the exact Gaussian density of `plot_wasserstein`. `seedmatch` does not export these names.
- `seedmatch.kmers.kmer_histogram`: the bins of the k-mer histogram.
- `seedmatch.report.html.report_html`: the HTML text of the report table. It uses only the `html` module of the standard library.

`seedmatch.plot` is a package with these modules:

- `__init__.py`: the exports, and an `ImportError` that names the extra `plot` if matplotlib is absent.
- `_style.py`: the palette, the line widths, the text of an R number, the size, limit, and color checks, and the two styles of Axes.
- `ecdf.py`, `wasserstein.py`, `kmers.py`, `alignment.py`, `ranges.py`, and `benchmark.py`: one public function each.

A plot function makes a `matplotlib.figure.Figure` with the object API. It does not import `matplotlib.pyplot`, and it does not change `matplotlib.rcParams`. It returns the figure and writes no file. The caller writes a file with `Figure.savefig`.

**Tech Stack:** Python 3.12, numpy, pandas, and matplotlib 3.11.1 of the extra `plot`. The standard library `html` module for the report table. R 4.5.3 with ggplot2 4.0.2, gt 1.4.0, scales 1.4.0, PRROC 1.4, and SeedMatchR 2.0.0 for the goldens.

**Spec:** `docs/design.md`, sections "Parity contract" ("Not identical by design"), "Architecture" ("Packages" and "Rules"), "Divergences (defect fixes)", "Tests" ("Coverage"), "Performance and bounds", and "Delivery order" (step 9).

**Later plans (not in this plan):** none.

**Closed items:**

- The PR and ROC plots. `docs/plans/2026-11-20-small-analyses.md:19` names them for delivery step 9, and Task 7 ports them. Task 9 closes the line of the small-analyses plan.
- The plot of a network. The network plan names it for delivery step 9 (`docs/plans/2026-11-20-network.md:27`). SeedMatchR has no function that draws a network: `grep` of `R/` and `vignettes/` finds no plot call in the network code. Thus there is no port, and Task 9 closes the line of the network plan.

## Global Constraints

- Runtime dependencies: `numpy`, `pandas`, `numba`, and `scipy`. The extra `plot` adds `matplotlib`. No other dependency. The port does not use jinja2, because pandas 2.3.3 rejects jinja2 3.0.3 of the store lock (see "Facts from the planner probes").
- `import seedmatch` does not import matplotlib or `seedmatch.plot`. The core functions of this plan are in the core, and `seedmatch` exports them. Only `seedmatch.plot` exports the plot functions.
- A plot function uses `matplotlib.figure.Figure` and its Axes. It does not import `matplotlib.pyplot`, and it does not change `matplotlib.rcParams`. It writes no file.
- The tests of a figure read the data of its artists: the points of each line, the bars, the texts, and the colors. They do not compare pixels. Pixel identity with ggplot2 is a non-goal of the design. The coverage omits `*/seedmatch/plot/*`, but each plot module has tests.
- Errors derive from `seedmatch.errors.SeedMatchError`. Invalid input raises `InputError`, and each message names the bad value and the expected value. A warning has the category `SeedMatchWarning`, and it goes through `seedmatch.errors.warn`.
- The message templates of this plan are normative. Each error test matches a stable part of its template. A template that names a helper of `seedmatch._checks` is the message of that helper.
- Each change of behavior from R gets an entry `## D-NNN: <title>` in `docs/divergences.md`. It also gets a row in the table of `docs/design.md`, and a test whose name contains the ID.
- This plan writes the entries of the new IDs `D-063` to `D-070`, after the entry of `D-062`. The table rows come after the row of `D-062`, in ID order. The IDs come in task order, thus each task adds its entries and its rows at the end. Task 4 also changes the text of `D-010`.
- A rename of an argument or a value, or a change of a type, is not a divergence. `docs/r-map.md` records it. The style of a figure (fonts, sizes, and margins) is not a divergence.
- The R source of truth is `../SeedMatchR` at commit `2fc0c4f`. The R scripts run with `tools/golden/run.sh`. The tests never run R and never use the internet.
- Licenses (`docs/design.md`, "Licenses"): ggplot2 and gt have the MIT license. The C source of `pretty` and `density` in R has the GPL license. The port does not copy it. `kmer_histogram` uses its own rule (`D-066`), and the density is the exact Gaussian sum.
- PRROC has the GPL-3 license. `plot_benchmark` does not translate `plot.PRROC`. Its rules come from the output of the planner probes.
- Each function obeys "Performance and bounds" of the design. Its docstring gives the time and the memory, and its task adds a line to that section of `docs/design.md`.
- Task 2 runs `uv lock` without `--offline`. This command reads the package index and changes only `uv.lock`. After it, `uv sync --locked --offline` must pass.
- Documents use ASD-STE100. Code comments obey the admission test of the project rules. Documents, names, and comments obey the banned-word rule of the project.
- Git: an implementer does no `git add` and no commit. Each task ends with a checkpoint. The controller commits.

## Facts from the planner probes

Each probe script and its output is in `tmp/claude/plan-plot/` of the planner. The R probes ran with `tools/golden/run.sh`. The Python probes ran with `uv run --with matplotlib==3.11.1` and `uv run --with jinja2==3.0.3`. Machine: Apple M2 Max, Python 3.12.

Versions and the lock:

- R 4.5.3, ggplot2 4.0.2, gt 1.4.0, scales 1.4.0, and SeedMatchR 2.0.0 (`probe_ecdf.out`). ggplot2 and scales are in the system library of R. gt is in `tools/golden/.rlib`.
- `../inflexa/images/package-store/lock.amd64.json` pins `matplotlib==3.11.1`, `contourpy==1.4.0`, `cycler==0.12.1`, `fonttools==4.65.0`, `kiwisolver==1.5.1`, and `pillow==12.3.0`. It also pins `pyparsing==3.3.2`, `python-dateutil==2.9.0.post0`, `six==1.17.0`, and `packaging==26.3`. `lock.arm64.json` gives the same versions.
- `manifest.yaml:436-438` makes matplotlib a root, with the warm script `warm/matplotlib.py`. This script fills the font cache of the sandbox.
- The lock pins `jinja2==3.0.3` and `markupsafe==3.0.3`. `manifest.yaml:592` names `Jinja2` as a root.
- `uv.lock` of this repository has python-dateutil 4.0.0, but the store lock has 2.9.0.post0. Each lock has packaging 26.3.

pandas `Styler` with the jinja2 of the lock (`probe_styler_303.out`, `probe_styler_direct.out`, `probe_styler_316.out`):

- pandas 2.3.3 with jinja2 3.0.3 and markupsafe 3.0.3: `DataFrame.style` raises `AttributeError: The '.style' accessor requires jinja2`.
- The import of the `Styler` module raises `ImportError: Pandas requires version '3.1.2' or newer of 'jinja2' (version '3.0.3' currently installed).`
- With jinja2 3.1.6, `DataFrame.style` gives a `Styler`, and its HTML has the element `<style type="text/css">`.
- `DataFrame.to_html` works without jinja2. `Styler.background_gradient` takes its colors from matplotlib, thus it also works only with the extra.
- Decision of the controller: `report_html` writes the HTML with the `html` module of the standard library. The other choice was jinja2 3.1.2 or later in the inflexa lock, a change of another repository.

matplotlib without pyplot (`proto_pyplot2.out`, `proto_pyplot3.out`, `proto_rc.out`, `proto_rc2.out`, `proto2.out`):

- These calls load no `matplotlib.pyplot`: `Figure(layout="constrained")`, `add_subplot`, `subplots`, `plot`, `legend`, `suptitle`, `fill_between`, `stairs`, `axvline`, `bar`, `Figure.text`, `FigureCanvasAgg(fig).get_renderer()`, and `savefig` to PNG and to PDF.
- `matplotlib.get_backend()`, `dict(matplotlib.rcParams)`, and `matplotlib.rcParams == other` resolve the backend, and that loads pyplot.
- `{k: v for k, v in dict.items(matplotlib.rcParams) if k != "backend"}` loads no pyplot, and a comparison of two such snapshots shows a change of a parameter.
- `matplotlib.colors.is_color_like("blu")` is `False`, and `is_color_like("#3498db")` is `True`. `to_hex` of the fill `#F8766D` with alpha 0.25 gives `#f8766d40`.
- `Line2D.get_xydata()` gives the points of the call, not the points of the steps. A legend with `loc="upper left"` has the code 2.

`.plot_ecdf` and `deseq_fc_ecdf` (`probe_ecdf.out`, `probe_ecdf_legend.out`). The test data has 179 genes: 18 genes with `mer7m8 >= 1` and 161 genes with `mer7m8 == 0`.

- The p-value labels of `.plot_ecdf` (`R/deseq_fc_ecdf.R:312-324`):

| p-value | label |
| --- | --- |
| 0.05 | `0.05` |
| 0.0149 | `0.01` |
| 0.0105 | `0.01` |
| 0.01 | error `object 'pvalue_label' not found` |
| 0.0099999 | `0.01` |
| 0.00999 | `0.01` |
| 0.005 | `0.005` |
| 0.0025 | `0.002` |
| 0.001 | `0.001` |
| 0.000999 | `p<0.001` |
| 1e-10 | `p<0.001` |
| 0.996 | `1` |
| 1 | `1` |
| NaN | error `missing value where TRUE/FALSE needed` |
| NA | error `missing value where TRUE/FALSE needed` |
| NULL | `Undefined` |

- `20 / 2000 == 0.01` is `TRUE` in R. Thus a DTS p-value of 20 of 2000 permutations is exactly 0.01.
- The layer is `GeomStep` with `StatEcdf` and `linewidth` 1.5. The colors are `black` and `#3498db`.
- The group `Background` has 163 layer rows: `-Inf`, 161 distinct values, and `Inf`. The group `mer7m8` has 20 rows. The rows come in the order of the data, not of x.
- The title is `ECDF`. The subtitle is `Statistic: 0.51\np-value: p<0.001`. The axis labels are `Log2(Fold Change)` and `Cumulative Fraction`.
- The coordinate limits are -1 and 1. The panel ranges are -1.1 to 1.1 for x and -0.05 to 1.05 for y.
- R gives the KS statistic 0.51035196687370576 and the p-value 0.00010377312440832125. The port `ecdf_test` gives 0.510351966873706 and 0.00010377312440832123.
- `get_guide_data(pl, "colour")` gives the legend `Background: 161` (black) and `mer7m8: 18` (`#3498db`), with the title `Group`. `get_labels()` of the scale gives no labels, thus the golden script uses `get_guide_data`.
- A row of `res` in no gene list gives a third line in `grey50` with the legend key `NA`. Its layer rows are `-Inf`, 5, and `Inf`, with y 0, 1, and 1.
- `factor_order = c("mer7m8", "Background")` gives the legend `mer7m8: 161` (black) and `Background: 18`. The counts come from the order of `gene_lists` (`R/deseq_fc_ecdf.R:112`, `:121`).
- A `Background` list with 2 absent genes and 1 repeated gene gives `Background: 164`.
- An NA fold change of a seed gene gives the build warning ``Removed 1 row containing non-finite outside the scale range (`stat_ecdf()`).`` The seed group then has 19 layer rows, but the legend stays `mer7m8: 18`.
- `theme = ggplot2::theme(legend.position = "bottom")` leaves the legend position `inside`. `.plot_ecdf` adds `SeedMatchR.theme`, not its argument `theme` (`R/deseq_fc_ecdf.R:349`).
- 10 groups with the 9 colors of the palette stop the build with `Insufficient values in manual scale. 10 needed but only 9 provided.`
- `output_filename` with the name `ecdf.png` gives a file that starts with `%PDF-`, because `ggsave` gets `device = "pdf"` (`R/deseq_fc_ecdf.R:357`).
- 12,000 background genes give the legend `Background: 12,000`.
- The subtitles of the other tests: Wilcoxen `Statistic: 748 | p-value: 1`, Wass `Statistic: 2.196 | p-value: Undefined`, and DTS after `set.seed(1)` `Statistic: 74.912 | p-value: p<0.001`. The DTS p-value is 0.00025.
- `paste` of `round(x, 3)` gives `0.51`, `1`, `0`, `74.912`, and `12345.679` (15 significant digits).
- `SeedMatchR.theme` (`R/themes.R:2-14`): text `#34495e` of size 12, titles of size 14, and major grid lines `#6b8cae` with the line type 2 and the width 0.2. The legend is inside at (0.01, 0.99) with the justification (0, 1), and its title is bold of size 12. The palette is at `R/themes.R:17-18`.

`wass_dist` (`probe_wass.out`, `probe_dump.out`). The seed values are the 18 fold changes with `mer7m8 >= 1`, and the non-seed values are the 161 others.

- The areas are 2.196, 2.196, and 0 with `n.interp = 10000`, and 2.197, 2.197, and 0 with 100.
- `cdf.plot`: `GeomLine` with `StatEcdf` and `linewidth` 1.25. The colors are `black` (non-seed) and `red` (seed). The labels are `log2(FC)` and `eCDF`, and there is no title.
- The `cdf.plot` rows of the non-seed group start with `-Inf`, -3.44283360609722, and -2.29233729598998. Their y values are 0, 0.0062111801242236, and 0.0124223602484472. Thus `geom_line` joins the left edge of the panel at 0 to the first value at its height.
- `dens.plot`: the fills are `#F8766D` (non-seed) and `#00BFC4` (seed), with alpha 0.25 and a black outline. The labels are `log2(FC)` and `Density`.
- Each density group has 512 rows from -8.1646145214532293 to 3.4227253013710501, the pooled range of the two samples. The maxima are 0.38874830264622323 (non-seed) and 0.80556303947789587 (seed).
- `bw.nrd0` gives 0.50955898190682714 (non-seed) and 0.20292326433035393 (seed). For the non-seed sample, the IQR is 3.012701092342454 and the standard deviation is 1.5642894884738532.
- A prototype with `r_quantile` and `np.std(ddof=1)` gives the seed bandwidth exactly. It gives 0.5095589819068272 for the non-seed bandwidth, 1 unit in the last place from R.
- `density()` of R uses an FFT approximation. On the same grid, the exact Gaussian sum differs from the R layer by a fraction of the maximum. The fraction is at most 4.0e-4 (non-seed) and 1.34e-3 (seed). `np.linspace(low, high, 512)` is not bit-equal to the grid of R.
- `diff.plot`: 10000 rows, `GeomLine` in black, the labels `log2(FC)` and `Seed - NonSeed eCDF`. The title is `Wass(- censored) = 2.196   Wass(+ censored) = 0   Wass (total) = 2.196`.
- A larger case gives the title `Wass(- censored) = 1000.01   Wass(+ censored) = 0   Wass (total) = 1000.01`. Thus `paste` gives 15 significant digits.
- The private `_grid` and `_interpolate` of `seedmatch.stats.wasserstein` give the R difference at `n.interp = 100` within 1.1e-16. The x values differ from R in the last bit. `np.linspace(low, high, 100)` equals `_grid(low, high, 100, 0, 100)`.

`kmer_hist` (`probe_kmer_aln.out`, `probe_titles.out`, `kmer_hist.pdf`). The input is `Test_seqs.fa` with `width = 7`.

- `kmer_hist` returns `integer(0)` and draws 2 pages. Page 1 has the title `K-mer gene frequencies` and the x label `Number of genes with k-mer`. Page 2 has `Histogram of gene.counts` and `gene.counts`.
- `abline` draws the target line on page 2 only. `lw = 3` matches `lwd` as a partial name, and the PDF draws the line with the width 2.25 points.
- There are 16384 k-mers, with gene counts from 0 to 87. The target `CTCTATA` has 18 genes.
- `hist(gene.counts, breaks = 50)` gives the breaks 0, 2, ..., 88 (44 bins), with bins closed on the right. The counts are `3394 1773 1308 1276 1465 1430 1411 1160 944 656 460 361 235 138 119 65 43 32 21 21 16 9 7 11 6 2 5 4 5 1 2 0 0 0 3 0 0 0 0 0 0 0 0 1`.
- The counts 0, 1, and 2 occur 915, 1227, and 1252 times. The first bin of R holds the three values, 3394 = 915 + 1227 + 1252.
- The targets `CUCUAUA`, `ctctata`, `CTCTATAA`, and `NOTAKMER` give `NA`, and R draws no line and gives no message.
- `pretty(c(1, 7), 50, min.n = 1)` gives the step 0.1, thus 60 bins of width 0.1. `pretty(c(12, 9876), 50, min.n = 1)` gives the step 200 from 0 to 10000.
- A prototype of the rule of `kmer_histogram` (`proto_core.out`) on the same counts:
  - `max_bins=50`: the width 2 and 44 bins. The first counts are 2142, 2242, 1457, 1289, and 1354. The last count is 1, and the sum is 16384.
  - `max_bins=10`: the width 10 and 9 bins. The counts start with 8484, 6021, 1570, 225, and 59, and end with 4, 1, and 1.
  - `max_bins=3`: the width 50, with the counts 16359 and 25.

`plotAlignment` with `output.type = "ggplot"` and the patch of `D-045` (`probe_kmer_aln.out`). The match is the first `mer7m8` match, `ENSRNOG00000016275` from 73 to 79.

- The label of `geom_text` has these lines:

```text
miRNA       3'-UUUUGUCACAAGAACGAGAUAUU-5'   
                ||||||||||||||||||||||      
target   5'-...GAAACAGUGUUCUUGCUCUAUAA...-3'
```

- The text has the family `mono`, the face `bold`, the size 3.866058 (the default of ggplot2, 11 points), and `hjust` and `vjust` 0.5. The theme is `theme_void`.
- The row names of the data.frame are `siRNA `, six spaces, and `target` (`R/plotAlignment.R:135-137`).
- Without the patch, the `ggplot` output stops with `nchar(s1) == nchar(s2) is not TRUE`.

`plotIRanges` and `plotRangesDF` (`probe_ranges.out`, `proto_ranges.out`). The table `eight` has the starts 10, 30, ..., 150 and the ends start + 6.

- `NAMESPACE` does not export the functions. The formals are `ir, xlim, main, col, ...` and `ir, xlim, main, col, bar_height, ...`.
- `plotRangesDF(eight)` gives 6 rectangles, because `head()` keeps 6 rows (`R/plotIRanges.R:20`). xmin is 10, 30, ..., 110, and ymin is 0.875, 1.875, ..., 5.875.
- The fill is `#595959FF`, and the plot has no title. The panel x range is 4.7 to 121.3.
- `xlim = c(0, 1000)`, `main = "My title"`, and `col = "red"` give the same layer, no title, and no coordinate limits. Thus `plotRangesDF` ignores these options.
- `bar_height = 0.8` gives ymin 0.6 and ymax 1.4 for row 1. Thus this option has an effect.
- The y value of a row is `as.integer(rownames(ir))` (`:17-18`). The subset `eight[c(3, 7, 8), ]` gives ymin 2.875, 6.875, and 7.875.
- The row names `a`, `b`, and `c` give the warning `NAs introduced by coercion`, and ymin and ymax are NA.
- `data.frame(start = c(5, 20), end = c(10, 20))` gives xmin 20 and xmax 20 for row 2. Thus the bar of a range of one base has the width 0.
- `plotIRanges(IRanges(start = c(5, 20, 40), end = c(10, 20, 60)))` draws blue arrows from 5 to 10 at y 1 and from 40 to 60 at y 3.
- For the range from 20 to 20, R gives the warning `zero-length arrow is of indeterminate angle and so skipped`. R draws no arrow for it.
- The titles of `plotIRanges` are `IRanges Plot`, `Position`, and `Ranges`. The y range is 0 to n + 1. In `plotIRanges`, `xlim`, `main`, and `col` have an effect.
- An empty `IRanges()` stops `plotIRanges` with `need finite 'xlim' values`. A data.frame stops it with `invalid time series parameters specified (1)`.
- A prototype with one `PolyCollection` loads no pyplot. For `eight`, it gives 8 paths. The first path has the vertices (9.5, 0.875), (16.5, 0.875), (16.5, 1.125), and (9.5, 1.125).
- The limits of the prototype are (2.15, 163.85) and (0, 9). The face color is `#0000ff`.

`plot` of PRROC 1.4 for the result of `benchmark_preds` (`probe_prroc_src.out`, `probe_prroc.out`, `probe_prroc2.out`, `proto_benchmark.out`):

- `R/benchmark.R:11-13` documents the result as "Two figures". But the function returns the curve objects of `pr.curve` and `roc.curve` (`:32-40`), and it draws nothing.
- `plot(res$PRC)` and `plot(res$ROC)` call `plot.PRROC`. Its defaults are `xlim = c(0, 1)`, `ylim = c(0, 1)`, `auc.main = TRUE`, `auc.type = "integral"`, `color = TRUE`, `lwd = 3`, and `legend = 4`.
- The patched `benchmark_preds` with `padj = 0.15` and `l2fc = 0` gives 63 positives in 179 rows. The PR curve has 106 rows, and the ROC curve has 4 rows. The thresholds go from 0 to 2.
- The titles are `PR curve\nAUC = 0.4723177` and `ROC curve\nAUC = 0.591954`. `format` gives 7 significant digits.
- The axis labels are `Recall` and `Precision` for PR, and `FPR` and `Sensitivity` for ROC. The plot region is -0.04 to 1.04 for x and for y.
- `segments` draws the curve from each row to the next row, with `lwd = 3`. The PR curve has 105 segments, and the ROC curve has 3.
- The color of a segment comes from the threshold t of its first row: `scale.color[round(1 + 99 * (t - min) / (max - min))]`. min and max are the extremes of the thresholds.
- The scale is `hsv(h = seq(0, 1, length = 100) * 0.8, s = 1, v = 1)`. Color 1 is `#FF0000`, color 50 is `#00FF60`, color 53 is `#00FF85`, and color 100 is `#CC00FF`.
- The PR segments have the colors `#FF0000`, `#00FF60`, and `#CC00FF`. The first PR segment goes from (1, 0.3519553) to the same point.
- The ROC segments go from (0, 0) to (0.01724138, 0) in `#CC00FF`, then to (0.03448276, 0.2222222) in `#00FF60`, and then to (1, 1) in `#FF0000`.
- `layout` puts a color bar of 2.5 cm to the right of each plot. `image` draws the 100 colors from the lowest to the highest threshold, and `axis(4)` labels the bar.
- The unpatched `benchmark_preds` gives the ROC title `ROC curve\nAUC = 0.3762679`, because its ROC curve gets the labels as `scores.class1` (`D-006`).
- For one distinct score s, the color range goes from 0.9 s to s + 0.09 s. For s = 2, the range is 1.8 to 2.18, and each segment is `#00FF85`.
- For s = 0, each color is NA. For s = 0 and s = -3, `image` stops R with `increasing 'x' and 'y' values expected`.
- With an infinite score, the color of the threshold `Inf` is NA. Then R stops with `'to' must be a finite number`.
- `round(1 + 99 * 0.5)` is 50 in R, because R rounds half to even. Thus the index from 0 is `np.rint(1 + position) - 1`. `np.rint(99 * 0.5)` gives 50, the index of the next color.
- A prototype with one `LineCollection` for each curve, and `Figure.colorbar` of a `ScalarMappable`, loads no pyplot and keeps the rcParams. It gives the segments and the colors of R.
- In the prototype, `fig.axes` has 4 Axes, because each color bar is an Axes. The PR Axes and the ROC Axes come first.

The `gt` table of `SeedMatchReport` on the test data (`probe_gt.out`, `probe_titles.out`, `gt_report.html`):

- The columns come in this order: the stub `Category`, then `Full-D0` to `Full-D4`, `18mer-D0` to `18mer-D4`, and `15mer-D0` to `15mer-D4`. Then `8mer`, `7mer-m8`, `7mer-A1`, `6mer`, and `Total`. `tab_spanner` gathers the `15mer` columns from D0 to D4, but `report` has them from D4 to D0.
- The labels of the window columns are `D0` to `D4`. The spanners are `Full Guide Strand (g2:g23)`, `18-mer (g2:g19)`, and `15-mer (g2:g19)`.
- The title is `In-silico siRNA Binding Prediction`, and the subtitle is `Identifying siRNA hits in the transcriptome`. Each has `colspan="21"`. The row group `SeedMatchReport` has a title row.
- The count rows have 0 decimals, for example `125` and `42`. `fmt_number` of 12345 gives `12,345`. The percent row gives `45.83%`, `100.00%`, `0.00%`, and `33.60%` (`Total`). `fmt_percent` of 0.125 gives `0.12%`.
- `data_color` colors only the category columns of the row `% off-target`, not `Total`:
  - 100% gives `#FF0000` with white text, and 75% gives `#FF6846` with white text.
  - 45.83% gives `#FFA68B`, 33.33% gives `#FFBFAA`, 28.77% gives `#FFC8B5`, and 22.22% gives `#FFD5C6`.
  - 0% gives `#FFFFFF`.
- The table element has `id="mygt"`, and the CSS selects `#mygt` in each rule. The rows 2 and 4 have the class `gt_striped`.
- The port `seed_match_report` on the test data gives these values:
  - `In silico predictions`: `Full-D1` 1, `Full-D4` 1, `18mer-D3` 3, `18mer-D4` 24, `15mer-D4` 73, `15mer-D3` 18, `8mer` 1, `6mer` 4, and `Total` 125. Each other category is 0.
  - `Off-target predictions`: `Full-D1` 1, `Full-D4` 1, `18mer-D3` 1, `18mer-D4` 11, `15mer-D4` 21, `15mer-D3` 4, `6mer` 3, and `Total` 42.
  - `% off-target`: 100 for `Full-D1` and `Full-D4`, 33.333333 (`18mer-D3`), 45.833333 (`18mer-D4`), 28.767123 (`15mer-D4`), 22.222222 (`15mer-D3`), 75 (`6mer`), and 33.6 (`Total`).

The README data (`tests/golden/stats/ecdf.tsv`, input `readme`): 17 target genes and 121 background genes. KS with `greater` gives 0.45989304812834275 and 0.0010735381780312581. `r_round` gives the labels `0.46` and `0.001`. The Wass statistic, `neg.cens.auc`, is 2.033.

Measurements of the planner (`proto.out`, `proto2.out`):

- The ECDF figure of 22,000 values with `savefig` to PNG at 100 dpi takes 0.022 seconds. A line of 10^6 points with `savefig` takes 0.056 seconds.
- The density of 20,000 values at 512 points takes 0.048 seconds with blocks of 4096 values, and the peak is 48 MiB. Blocks of 1024 values give a work array of 512 · 1024 · 8 bytes, 4 MiB.
- The Wasserstein figure of 20,000 and 2,000 values with `n_interp=10000` and blocks of 1024 takes 0.072 seconds. With `savefig` to PNG at 100 dpi, it takes 0.463 seconds. The peak of `tracemalloc` is 8.8 MiB.
- The range prototype draws each bar in one `PolyCollection`. With 10^5 rows, it takes 0.705 seconds, and 4.228 seconds with `savefig` to PNG at 100 dpi.
- The peak of `tracemalloc` for 10^5 rows is 43.0 MiB. Thus each row costs 42 µs and 451 bytes with `savefig` (`proto_benchmark.out`). With 10^4 rows, a row costs 51 µs and 475 bytes.
- The benchmark prototype with 10^5 distinct scores has 200,003 curve rows. It takes 1.895 seconds, and 8.974 seconds with `savefig`, and the peak is 71.4 MiB.
- The benchmark prototype with 10^6 scores of 8 values has 262,604 curve rows, because the PR curve has inner rows for about 300,000 positives. It takes 2.571 seconds, 11.864 seconds with `savefig`, and 98.4 MiB.
- The alignment text with `left_pad=20` and `right_pad=20` has lines of 84 characters. In a figure of 0.09 · 84 + 0.4 by 0.9 inches with the font size 10, the text is inside the figure.

## Shared interfaces

A task makes each of these names, and a later task uses it:

```python
# seedmatch/plot/_style.py (Task 2)
PALETTE: Final = (
    "black", "#3498db", "#2ecc71", "#f1c40f", "#e74c3c", "#9b59b6", "#1abc9c", "#f39c12", "#d35400",
)
TEXT_COLOR: Final = "#34495e"

def linewidth_points(linewidth: float) -> float: ...
def pasted(value: float) -> str: ...
def expanded(low: float, high: float) -> tuple[float, float]: ...
def as_size(size: object) -> tuple[float, float]: ...
def as_limits(name: str, value: object) -> tuple[float, float]: ...
def as_color(name: str, value: object) -> str: ...
def new_figure(size: tuple[float, float]) -> Figure: ...
def style_theme_axes(axes: Axes) -> None: ...
def style_gray_axes(axes: Axes) -> None: ...
def as_written[T: Text](text: T) -> T: ...

# seedmatch/stats/wasserstein.py (Task 3): renamed from _as_wasserstein_sample, with the same messages.
def as_wasserstein_sample(name: str, values: ArrayLike) -> npt.NDArray[np.float64]: ...

# seedmatch/stats/wasserstein.py (Task 3): renamed from _Ecdf and _ecdf, with the same fields.
class Ecdf(NamedTuple): ...  # knots, heights, widths, steps
def ecdf(sample: npt.NDArray[np.float64]) -> Ecdf: ...

# seedmatch/stats/density.py (Task 3)
DENSITY_POINTS: Final = 512

def nrd0_bandwidth(values: npt.NDArray[np.float64]) -> float: ...
def gaussian_density(
    values: npt.NDArray[np.float64], grid: npt.NDArray[np.float64], bandwidth: float
) -> npt.NDArray[np.float64]: ...
```

Rules of `_style.py`:

- `PALETTE` is `SeedMatchR.palette` (`R/themes.R:17-18`).
- `linewidth_points` gives the width in points of a ggplot2 `linewidth`: `linewidth * 72.27 / 25.4 * 72 / 96`. Thus 1.5 gives 3.2009, 1.25 gives 2.6674, and 0.5 gives 1.0670. A base graphics `lwd` of 1 is 0.75 points (`kmer_hist.pdf`).
- `pasted` gives the text of `paste` of R for a number: `format(value, ".15g")`, with `-0` as `0`. Thus 0.51 gives `0.51`, 748.0 gives `748`, 1e-05 gives `1e-05`, and 1000.01 gives `1000.01`. R gives `1e+05` for 100000, but `pasted` gives `100000`, as `docs/r-map.md` records.
- `expanded` gives the limits of ggplot2 with the default expansion: `(low - 0.05 * (high - low), high + 0.05 * (high - low))`.
- `as_size` returns a pair of finite numbers above 0 as floats. Another value raises `size (0, 5.0) is not valid; expected (width, height) in inches, two finite numbers above 0`. The message gives the value as `repr` of the input, thus `(0, 5)` gives `size (0, 5)`. A bool is not a number.
- `as_limits` returns a pair of finite numbers with low < high as floats. Another value raises `x_limits (1.0, -1.0) is not valid; expected (low, high), two finite numbers with low < high`. The message starts with `name`, and a bool is not a number.
- `as_color` returns the value if it is a `str` that `matplotlib.colors.is_color_like` accepts. A value that is not a `str` raises the message of `as_text`.
- Another text raises `palette[1] 'blu' is not a matplotlib color; expected a color name or a hex code, for example '#3498db'`. The message starts with `name`.
- `new_figure` returns `Figure(figsize=size, layout="constrained", facecolor="white")`.
- `style_theme_axes` applies the values of `SeedMatchR.theme` to one Axes:
  - the axis labels in `TEXT_COLOR` with the size 12, and the tick labels with the size 9.6
  - a white panel with an edge of `#2c3e50`
  - major grid lines of `#6b8cae`, dashed, with the width `linewidth_points(0.2)`
- `style_gray_axes` applies the panel of `theme_grey` of ggplot2: the face color `#EBEBEB`, white grid lines with the width `linewidth_points(0.5)`, and no edge lines.
- `as_written` sets `parse_math=False` on a Text and returns it. Each plot function calls it on each title or label that comes from an argument. Thus a `$` in a user text draws as written.
- The sizes are near the sizes of R, not equal. These functions change only the given Axes, never `rcParams`.

Rules of `seedmatch.stats.density`:

- `nrd0_bandwidth` gives `0.9 * min(sd, iqr / 1.34) * n ** -0.2`, the `bw.nrd0` rule of R. `sd` is `np.std(values, ddof=1)`, and `iqr` is `r_quantile(values, 0.75) - r_quantile(values, 0.25)`. If `iqr` is 0, the rule uses `sd`. A sample of `as_wasserstein_sample` has two distinct values, thus `sd` is above 0.
- `gaussian_density` gives, for each grid point g, the sum over the values v of `exp(-0.5 * ((g - v) / h) ** 2)`, divided by `n * h * sqrt(2 * pi)`. It reads the values in blocks of the module constant `_BLOCK = 1024`, and it reads `_BLOCK` at each call. Thus the work memory is `grid.size * 1024 * 8` bytes, and the time is O(`grid.size` · n).

## Review Focus

1. The HTML of the report without `Styler`. The controller chose the `html` module of the standard library, because pandas 2.3.3 rejects jinja2 3.0.3 of the store lock. Owner: Task 8.
2. No `matplotlib.pyplot`, and no change of `rcParams`, after each plot function and `savefig`. The snapshot must use `dict.items` without the key `backend`, because the other forms load pyplot. Owner: Task 2.
3. The bins of the k-mer histogram are not the bins of R (`D-066`). The first bar is 2142, and the first bar of R is 3394. Owner: Task 4.
4. The density agrees with R within 0.5% of its maximum, because R uses an FFT approximation. The planner measured 0.04% and 0.134%. Owner: Task 3.
5. The x limits and the ends of the lines. `plot_ecdf` and `plot_wasserstein` use the 5% expansion of ggplot2. The ECDF lines end at the limits, where R draws to the edge of the panel. Owner: Task 2 and Task 3.
6. The colors of the percent row are a linear blend from white to red, with black text. `gt` blends in another color space, and it gives white text on a dark cell. 0% and 100% give the colors of `gt`. Owner: Task 8.
7. The diff of `uv.lock`. It must add matplotlib and its dependencies at the versions of the store lock, and change python-dateutil from 4.0.0 to 2.9.0.post0. Then `uv sync --locked --offline` must pass. Owner: Task 2.
8. The times of each figure, and the width of the alignment figure with the font size 11. Owner: each task.
9. The cost of `plot_ranges` for a large table. The controller ruled that the function has no row limit, because a limit drops rows without a message, as `head()` does. One `PolyCollection` holds each bar. The planner measured 42 µs and 451 bytes for each row with `savefig`. Owner: Task 6.
10. The rows of `plot_ranges` (`D-068`). The port draws each row at its position in the table, from `start - 0.5` to `end + 0.5`. R draws 6 rows at the positions of the row names. Owner: Task 6.
11. The cost of `plot_benchmark`. It draws each row of the two curves, and the PR curve has up to 2K + 2 + max(100, P) rows. The planner measured 11.9 seconds and 98 MiB for 262,604 rows with `savefig`. Owner: Task 7.
12. The colors of `plot_benchmark` (`D-069`). The index of a color must use the rounding of R, half to even, after the addition of 1. A curve without two distinct finite thresholds gets the middle color. Owner: Task 7.

---

### Task 1: R goldens of the plots and of the report table

**Files:**
- Add: `docs/plans/2026-11-21-plot-and-report.md` (this plan)
- Make: `tools/golden/R/plots.R`, and the files of `tests/golden/plots/` (written by the script)
- Modify: `tools/golden/R/check_env.R`, `tools/golden/README.md`

**Interfaces:**
- Produces: the goldens below. Each later task reads them with `tests/helpers.py::golden("plots/<file>")` and `golden_note("plots/<file>", pattern)`.

- [ ] **Step 1: Write `tools/golden/R/plots.R`.** The script sources `common.R`, and it writes each golden with `write_golden` into `tests/golden/plots/`. The notes of each golden give the versions of ggplot2, gt, scales, and PRROC. These are the rules of the script:
  - The first statement after `source` is `grDevices::pdf(NULL)`. `ggplot_build` and `hist` open a device, and without this statement R writes `Rplots.pdf` into the root of the repository.
  - `library(ggplot2)` comes before the plots. `load_test_data()` gives `res`, `seqs`, and `guide`. The seed genes are `res$gene_id[res$mer7m8 >= 1]`, and the background genes are `res$gene_id[res$mer7m8 == 0]`.
  - A note has the form `<key>: <value>`. A line end in a value becomes ` | `.
- The sections of the script:
  - **ECDF lines.** `ecdf_lines.tsv` holds the columns `case`, `group`, `x`, and `y`. These are the rows of `ggplot_build(pl)$data[[1]]`, sorted by x in each group, as the step line goes. `group` is the legend key of the row: `Background`, `mer7m8`, or `NA`. The cases:
    - `test_data`: `deseq_fc_ecdf(res, list(Background = background, mer7m8 = seed_genes))`.
    - `no_group_r`: the same call, with one more row of `res`. The row is a copy of row 1, with the `gene_id` `ENSNOGROUP` and the `log2FoldChange` 5.
  - Notes of the ECDF section, in `ecdf_lines.tsv`. Each legend comes from `get_guide_data(pl, "colour")` as `<label> [<colour>]`, with ` / ` between the keys:
    - `legend test_data`, `legend no_group_r`, `legend factor_order_reversed` (`factor_order = c("mer7m8", "Background")`), and `legend absent_and_repeated` (two absent genes and the first background gene again in `Background`).
    - `legend na_fold_change`: the `log2FoldChange` of the first seed gene is `NA`. Also `na_fold_change rows of mer7m8` and the build warning.
    - `legend big_counts`: a synthetic `res` of 12,000 background genes and 345 target genes.
    - `subtitle KS`, `subtitle Wilcoxen`, `subtitle Wass`, and `subtitle DTS`. The DTS call comes after `set.seed(1)`.
    - `labels` (the title and the axis labels), `layer` (the geom, the stat, the line width, and the colors), `coord limits`, and `panel range` (x and y of `layout$panel_params`).
    - `theme argument`: the legend position with `theme = ggplot2::theme(legend.position = "bottom")`.
    - `ten groups`: the error of the build with 10 gene lists of `res$gene_id`.
    - `ggsave to a png name`: the first 5 bytes of the file of `output_filename = file.path(tempdir(), "ecdf.png")`.
    - `theme`: the values of `SeedMatchR.theme` in "Facts from the planner probes".
  - **ECDF labels.** `ecdf_labels.tsv` holds the columns `p_value` and `label`. For each p-value of the table in "Facts from the planner probes", call `.plot_ecdf` with `test_results = list(statistic = c(D = 0.123456), p.value = p)`. `label` is the text after `p-value: ` of the subtitle, or `ERROR: <message>`. `p_value` is the text of `double_to_string`, or `NaN`, `NA`, or `NULL`.
  - **Wasserstein.** The cases are `test_data` (the seed and the non-seed fold changes) and `wass_40x300_ties` (x is the seed sample, and y is the non-seed sample, of `tests/golden/stats/samples.tsv`). For each case, `wass_dist(seed, non.seed)` gives the plots:
    - `wasserstein_ecdf.tsv`: `case`, `group`, `x`, and `y`, the rows of `cdf.plot` sorted by x in each group. `group` is `non-seed` or `seed`.
    - `wasserstein_density.tsv`: `case`, `group`, `x`, and `density`, the 512 rows of each group of `dens.plot`.
    - `wasserstein_difference.tsv`: `case`, `n_interp`, `x`, and `difference`, the rows of `diff.plot` of `wass_dist(..., n.interp = 100)`.
  - Notes of the Wasserstein section, in `wasserstein_density.tsv`: `bandwidth <case> <group>` (`format(bw.nrd0(x), digits = 17)`), and `title <case> n_interp 100` and `title <case> n_interp 10000`. Also `panel x range <case>`, and the geoms, the colors, the fills, the alpha, the line widths, and the labels of each plot.
  - **K-mer histogram.** The script computes `gene.counts` with the calls of `R/kmer_hist.R:19-28` on `seqs` with `width = 7`. `kmer_hist.tsv` holds the columns `gene_count` and `kmers` of `as.data.frame(table(gene.counts))`, with `gene_count` as an integer. The notes:
    - `hist breaks` and `hist counts` of `hist(gene.counts, breaks = 50, plot = FALSE)`, with `,` between the values.
    - `target CTCTATA` and the targets `CUCUAUA`, `ctctata`, `CTCTATAA`, and `NOTAKMER`: `gene.counts[target][[1]]`.
    - `return value`, `pages`, and `title calls` of `kmer_hist(seqs, "CTCTATA")`. The script draws to `pdf(file.path(tempdir(), "kmer.pdf"), compress = FALSE)`, and it traces `graphics::title` as the planner probe does.
    - `abline width`: the PDF holds the text `2.25 w`.
  - **Alignment.** `alignment_plot.tsv` holds the columns `case`, `line`, and `text`. The case `mer7m8_1` is the match of the case `mer7m8_1` of `tests/golden/analyses/alignment.tsv`, `ENSRNOG00000016275` from 73 to 79. Call `plotAlignment(..., output.type = "ggplot")` inside `with_seedmatchr_patch("tools/golden/patches/plotAlignment.R", ...)`. The lines are `strsplit` of the label of the layer data.
  - Notes of the alignment section: `text` (the family, the face, the size, `hjust`, and `vjust`), `theme`, `data.frame row names` in brackets, and `unpatched` (the error without the patch).
  - **Ranges.** `ranges_plot.tsv` holds the columns `case`, `xmin`, `xmax`, `ymin`, and `ymax`. These are the rows of `ggplot_build(SeedMatchR:::plotRangesDF(ir))$data[[1]]`. The table `eight` is the table of "Facts from the planner probes". The cases:
    - `eight_rows`: `eight`.
    - `subset`: `eight[c(3, 7, 8), ]`.
    - `text_row_names`: `eight[1:3, ]` with the row names `a`, `b`, and `c`.
    - `one_base`: `data.frame(start = c(5, 20), end = c(10, 20))`.
  - Notes of the ranges section:
    - `options`: the title and the coordinate limits of `plotRangesDF(eight, xlim = c(0, 1000), main = "My title", col = "red")`, and the number of layer rows.
    - `fill`, `bar_height 0.8` (ymin and ymax of row 1), and `text_row_names warning`.
    - `plotIRanges arrows`, `plotIRanges titles`, and `plotIRanges one_base warning` of `IRanges(start = c(5, 20, 40), end = c(10, 20, 60))`. The script traces `graphics::arrows` and `graphics::title`, as the planner probe does.
    - `plotIRanges empty` and `plotIRanges data.frame`: the errors of `IRanges::IRanges()` and of `eight`.
  - **PR and ROC.** `benchmark_plot.tsv` holds the columns `case`, `curve`, `segment`, `x0`, `y0`, `x1`, `y1`, and `color`. `curve` is `PR` or `ROC`.
    - The rows are the arguments of `graphics::segments` in `plot(res$PRC)` and `plot(res$ROC)`. The script traces `graphics::segments`, `graphics::title`, `graphics::image.default`, and `graphics::axis`.
    - Do not trace `graphics::layout`. A tracer that reads its default `widths` stops `layout(1)` with `invalid type (NULL) for 'times'`.
    - The case `benchmark_1` is `benchmark_preds(res, padj = 0.15, l2fc = 0)` inside `with_seedmatchr_patch("tools/golden/patches/benchmark_preds.R", ...)`, as in `analyses.R`.
    - The case `one_score_2` is `pr.curve(scores.class0 = rep(2, 20), weights.class0 = rep(c(1, 0), 10), curve = TRUE)`, with only the curve `PR`.
  - Notes of the PR and ROC section:
    - `title PR`, `title ROC`, `labels PR`, and `labels ROC`: the arguments of `title`. A line end becomes ` | `.
    - `lwd`, `usr` (the value of `par("usr")` in the tracer of `segments`), `axis sides`, and `scale` (the colors 1, 50, 53, and 100).
    - `color bar PR` and `color bar ROC`: the first and the last y value and the number of colors of `image`.
    - `title ROC unpatched`: the ROC title of `benchmark_preds(res)` without the patch.
    - `one_score_2 color bar`, and `one score 0`, `one score -3`, and `an Inf score`. Each of the last three gives the error of `plot` and the colors of `getColor`. The inputs are `rep(0, 20)`, `rep(-3, 20)`, and `c(Inf, rep(1, 9), rep(0, 10))`, each with `weights.class0 = rep(c(1, 0), 10)`.
  - **Report table.** `report_table.tsv` holds the columns `row`, `column`, `text`, `background`, and `color`. Make `tab <- SeedMatchReport(res, seqs, guide)$table`. Then read each `<td headers=...>` cell of `gt::as_raw_html(tab, inline_css = FALSE)`.
    - `row` is the stub text of the row. `column` is the column name, without the prefix `a` that gt adds to an ID that starts with a digit.
    - `background` and `color` come from `gt:::dt_styles_get(gt:::build_data(tab, context = "html"))`, or `NA` without a style.
  - Notes of the report section: `title`, `subtitle`, `group`, `columns` (the order of the cells), and `colspan`. Also `spanner Full`, `spanner 18mer`, and `spanner 15mer`, each with its label and its columns. Also `id`, `css selector` (the first selector of the CSS), and `striped rows`.
  - More notes of the report section: `spanner Full of a 21-base guide` (`SeedMatchReport(res, seqs, substr(guide, 1, 21))`), `fmt_number 12345`, and `fmt_percent 0.125`.
- [ ] **Step 2: Pin the versions.** In `check_env.R`, add `ggplot2 = "4.0."`, `gt = "1.4."`, and `scales = "1.4."` to `pinned`, and change the comment above it to name `plots.R`. Run `tools/golden/run.sh tools/golden/R/check_env.R`. Expected: the last line names `ggplot2 4.0.2`, `gt 1.4.0`, and `scales 1.4.0`.
- [ ] **Step 3: Run the script two times.** Expected: the second run of `tools/golden/run.sh tools/golden/R/plots.R` changes no byte, and `git status --short` shows no `Rplots.pdf`. Also expected:
  - `ecdf_lines.tsv`, case `test_data`: 163 rows of `Background` and 20 rows of `mer7m8`. The first row of each group has x `-Inf` and y 0. The second row of `mer7m8` has x -8.16461452145323 and y 0.055555555555555552.
  - `ecdf_lines.tsv`, case `no_group_r`: 3 rows of the group `NA`, with x `-Inf`, 5, and `Inf`.
  - `ecdf_labels.tsv`: the row `0.01` gives `ERROR: object 'pvalue_label' not found`, and the row `NULL` gives `Undefined`.
  - `wasserstein_difference.tsv`, case `test_data`: 100 rows. Row 1 has x -8.1646145214532293 and the difference 0.055555555555555552.
  - The notes `title test_data n_interp 100` and `title test_data n_interp 10000` give the areas 2.197 and 2.196. The case `wass_40x300_ties` gives 0.471, 0.061, and 0.532 for each `n_interp`.
  - The bandwidth notes of `test_data` give 0.50955898190682714 and 0.20292326433035393.
  - `kmer_hist.tsv`: the sum of `kmers` is 16384. `gene_count` 0, 1, and 2 have 915, 1227, and 1252. The note `target CTCTATA` gives 18, and the note `pages` gives 2.
  - `alignment_plot.tsv`: the three lines of "Facts from the planner probes".
  - `ranges_plot.tsv`: the case `eight_rows` has 6 rows, with xmin 10 to 110 and ymin 0.875 to 5.875. The case `subset` has ymin 2.875, 6.875, and 7.875.
  - `ranges_plot.tsv`: the case `text_row_names` has NA in ymin and ymax. Row 2 of the case `one_base` has xmin 20 and xmax 20.
  - The note `options` gives no title and no coordinate limits. The note `plotIRanges empty` gives `need finite 'xlim' values`.
  - `benchmark_plot.tsv`, case `benchmark_1`: 105 rows of `PR` and 3 rows of `ROC`. Row 3 of `ROC` goes from (4/116, 14/63) to (1, 1), with the color `#FF0000`. The test data has 116 negatives.
  - The note `title PR` gives `PR curve | AUC = 0.4723177`, and the note `title ROC unpatched` gives `ROC curve | AUC = 0.3762679`.
  - The note `one score 0` gives `increasing 'x' and 'y' values expected`, and the note `an Inf score` gives `'to' must be a finite number`.
  - `report_table.tsv`: 80 rows. The row `% off-target` and the column `18mer-D4` give `45.83%` and `#FFA68B`. The column `6mer` gives `75.00%`, `#FF6846`, and `#FFFFFF`. The row `In silico predictions` and the column `Total` give `125`.
  - The note `spanner Full of a 21-base guide` gives `Full Guide Strand (g2:g23)`, because the label is a fixed text (`R/SeedMatchReport.R:158`).
- [ ] **Step 4: Document the goldens.** Add a section "Golden files of the plots and the report table" to `tools/golden/README.md`. Give the command, the device rule, the cases, the columns, and the patch of `plotAlignment`. Checkpoint.

### Task 2: The extra `plot`, the plot package, and the ECDF figure

**Files:**
- Modify: `pyproject.toml`, `uv.lock`, `tests/test_package.py`, `docs/r-map.md`, `docs/divergences.md` (`D-063`, `D-064`), `docs/design.md`
- Make: `src/seedmatch/plot/__init__.py`, `src/seedmatch/plot/_style.py`, `src/seedmatch/plot/ecdf.py`, `tests/plot/test_plot_package.py`, `tests/plot/test_plot_ecdf.py`

**Interfaces:**
- Consumes: `EcdfTestResult` and `assign_groups` of `seedmatch.stats.ecdf`, `r_round`, `as_items`, `as_text`, `check_type`, `check_columns`, `check_numeric_columns`, and the goldens `ecdf_lines.tsv` and `ecdf_labels.tsv`.
- Produces: the names of `_style.py` in "Shared interfaces", and:

```python
# Each keyword-only argument is an independent option of the figure of `.plot_ecdf`.
def plot_ecdf(  # noqa: PLR0913
    de: pd.DataFrame,
    groups: Mapping[str, Iterable[str]],
    result: EcdfTestResult,
    *,
    title: str = "ECDF",
    palette: Sequence[str] = PALETTE,
    x_limits: tuple[float, float] = (-1.0, 1.0),
    size: tuple[float, float] = (5.0, 5.0),
) -> Figure: ...
```

Rules of the extra and the package:
- `pyproject.toml`: add `[project.optional-dependencies]` with `plot = ["matplotlib>=3.8"]`. Add `"matplotlib>=3.8"` to the `dev` group. 3.8 is the first version with `Figure.get_suptitle`.
- Add these lines to `constraint-dependencies`: `"matplotlib==3.11.1"`, `"contourpy==1.4.0"`, `"cycler==0.12.1"`, `"fonttools==4.65.0"`, `"kiwisolver==1.5.1"`, `"pillow==12.3.0"`, `"pyparsing==3.3.2"`, `"python-dateutil==2.9.0.post0"`, and `"six==1.17.0"`. These are the versions of the store lock.
- `seedmatch/plot/__init__.py` imports matplotlib first. If the import fails, it raises `ImportError("seedmatch.plot requires matplotlib; install the extra 'plot', for example: uv pip install 'seedmatch[plot]'")` from the first error.
- The package exports `PALETTE` and each plot function. `seedmatch/__init__.py` does not import `seedmatch.plot`.

Rules of `plot_ecdf` (`D-063`, `D-064`):
- `de` goes through `check_columns(de, ("gene_id", "log2FoldChange"))` and `check_numeric_columns(de, ["log2FoldChange"])`.
- `groups` without a key raises `groups has no group; expected at least one group`. Then `assign_groups(de, groups)` gives the group of each row, with its messages and its warning.
- `result` goes through `check_type("result", result, EcdfTestResult)`, and `title` through `as_text("title", title)`.
- `palette` goes through `as_items("palette", palette, as_color, "str")`. Thus a color that `is_color_like` rejects raises `palette[1] 'blu' is not a matplotlib color; expected a color name or a hex code, for example '#3498db'`.
- A palette with fewer colors than groups raises `palette has 1 color(s) for 2 groups; expected a color for each group`. R stops with `Insufficient values in manual scale` for this input.
- `x_limits` goes through `as_limits("x_limits", x_limits)`, and `size` through `as_size`.
- The values of a group are the finite `log2FoldChange` values of its rows. A group without such a value raises `group 'mer7m8' has no row with a finite log2FoldChange in the DE table; expected at least 1 row`.
- A row in no group has no line (`D-064`). The lines come in the order of the keys of `groups`. Line i has the color `palette[i]`.
- The view is `expanded(*x_limits)` for x, and (-0.05, 1.05) for y. Let lo be the lowest of the view and of each value of each group, and hi the highest.
- Let u be the distinct values of a group, sorted, and F the count of the values up to each u divided by n. The line has the points `[lo, *u, hi]` and `[0, *F, 1]`, with `drawstyle="steps-post"` and the width `linewidth_points(1.5)`.
- The label of a line is `f"{name}: {n:,}"`, with n the number of values of the line (`D-064`). The legend has the title `Group`, a bold title of the size 12, `loc="upper left"`, a white face, and an edge of `#2c3e50`.
- `fig.suptitle(title)` has the size 14 and `TEXT_COLOR`. The left title of the Axes is `f"Statistic: {pasted(r_round(result.statistic, 3))}\np-value: {label}"`, with the size 12.
- The p-value label (`D-063`): `Undefined` for `None`, and `NaN` for NaN. A p-value of 0.01 or more gives `pasted(r_round(p, 2))`. From 0.001 to below 0.01, it gives `pasted(r_round(p, 3))`. Below 0.001, it gives `p<0.001`.
- The axis labels are `Log2(Fold Change)` and `Cumulative Fraction`. `style_theme_axes` styles the Axes. The function does not change `de`.
- The time is O(n log n) for n rows, because `np.unique` sorts the values of each group. The figure holds two floats for each distinct value.

The expected results for `plot_ecdf(sirna_de(), {"Background": background, "mer7m8": seed_genes}, ecdf_test(...))` with the KS defaults:
- The figure has one Axes with 2 lines. The x limits are (-1.1, 1.1), and the y limits are (-0.05, 1.05).
- lo is -8.16461452145323 and hi is 3.42272530137105, the lowest and the highest fold change of the data.
- The `mer7m8` line has 20 points. The first two are (lo, 0.0) and (-8.16461452145323, 0.05555555555555555), and the last is (hi, 1.0).
- The `Background` line has 163 points. The second point is (-3.44283360609722, 0.006211180124223602).
- The colors are `#000000` and `#3498db`, and the width is 3.2009 points.
- The legend texts are `Background: 161` and `mer7m8: 18`, with the title `Group`.
- The suptitle is `ECDF`, and the left title is `Statistic: 0.51\np-value: p<0.001`.

- [ ] **Step 1: Change the project.** Change `pyproject.toml` as the rules give. Then run `uv lock` without `--offline`. This command reads the package index and changes only `uv.lock`.
- Then run `uv sync --locked --offline`. It must pass. The planner found a wheel of each store version in the uv cache. Expected:
  - `git diff uv.lock` adds matplotlib 3.11.1 and its dependencies at the versions of the store lock.
  - The diff changes python-dateutil from 4.0.0 to 2.9.0.post0, and no other package changes.
  - pandas uses python-dateutil. The sandbox gets each runtime dependency from the store lock. Thus with 2.9.0.post0, the dev environment agrees with the sandbox.
  - Give the diff summary in the report of the task. If `uv sync --locked --offline` fails, stop and give its error in the report.
- [ ] **Step 2: Write the tests that fail.** In `tests/test_package.py`, add `'seedmatch.plot'` to the prefixes of `test_import_loads_no_fetch_and_no_plot_modules`. In `tests/plot/test_plot_package.py`:
  - `test_plot_package_names_the_extra_without_matplotlib`: a subprocess sets `sys.modules["matplotlib"] = None` and imports `seedmatch.plot`. The error output has `ImportError` and `install the extra 'plot'`.
  - `test_plot_functions_load_no_pyplot`: a subprocess calls each plot function on small synthetic data, and writes each figure to a `BytesIO` as PNG. Then it prints `"matplotlib.pyplot" in sys.modules`, and the output is `False`. Each later task adds its function to the script.
  - `test_plot_functions_keep_the_rc_params`: take the snapshot `{k: v for k, v in dict.items(matplotlib.rcParams) if k != "backend"}` before and after each plot function and its `savefig`. The snapshots are equal.
  - `test_pasted_gives_the_text_of_r`: 0.51, 748.0, 1e-05, 1000.01, and -0.0 give `0.51`, `748`, `1e-05`, `1000.01`, and `0`.
  - `test_linewidth_points_converts_the_ggplot2_width`: 1.5, 1.25, and 0.5 give 3.2009, 2.6674, and 1.0670, within 1e-4.
  - `test_as_size_rejects_a_bad_size(size)`: `(0, 5)`, `(5, -1)`, `(5, math.inf)`, `(True, 5)`, `(5,)`, and `"5x5"`.
- In `tests/plot/test_plot_ecdf.py`, with the KS result above as the default:
  - `test_plot_ecdf_draws_the_r_steps`: the expected results above. Each line without its first and its last point equals the finite golden rows of its group. The golden is the case `test_data` of `ecdf_lines.tsv`, and the relative tolerance is 1e-12. The drawstyle is `steps-post`.
  - `test_plot_ecdf_gives_the_r_titles_and_labels`: the titles, the axis labels, and the limits. The notes `labels`, `subtitle KS`, `coord limits`, and `panel range` of `ecdf_lines.tsv` give the R values.
  - `test_plot_ecdf_gives_the_r_legend`: the note `legend test_data` gives the same texts and colors.
  - `test_plot_ecdf_labels_each_p_value_D063(row)`: parametrized over `ecdf_labels.tsv` without the row `NA`. The result is `dataclasses.replace(result, p_value=p)`, with `None` for `NULL`. A row with a label gives the same label. The row `0.01` gives `0.01`, and the row `NaN` gives `NaN`.
  - `test_plot_ecdf_labels_the_statistic_of_each_test(test, statistic, p_value)`: Wilcoxon gives `748` and `1`, Wasserstein gives `2.196` and `Undefined`, and DTS with `rng=0` gives `74.912` and `p<0.001`. The notes `subtitle Wilcoxen`, `subtitle Wass`, and `subtitle DTS` give the same texts.
  - `test_plot_ecdf_draws_only_the_groups_D064`: a `de` with the row of `ENSNOGROUP` gives 2 lines, and the legend `Background: 161` and `mer7m8: 18`. The case `no_group_r` has 3 groups, and its note gives the key `NA`.
  - `test_plot_ecdf_counts_the_values_of_each_line_D064`: two absent genes and a repeated gene in `Background` give `Background: 161`. The note `legend absent_and_repeated` gives `Background: 164`. A NaN fold change of the first seed gene gives `mer7m8: 17` and 19 points.
  - `test_plot_ecdf_follows_the_order_of_the_keys_D064`: `{"mer7m8": ..., "Background": ...}` gives the legend `mer7m8: 18` and `Background: 161`. The black line has 20 points. The note `legend factor_order_reversed` gives `mer7m8: 161`.
  - `test_plot_ecdf_puts_a_comma_in_a_large_count`: a synthetic `de` of 12,000 and 345 genes gives `Background: 12,000`, as the note `legend big_counts` gives.
  - `test_plot_ecdf_extends_the_steps_to_the_view`: `x_limits=(-10.0, 10.0)` gives the x limits (-11.0, 11.0). The first point of each line has x -11.0, and the last point has x 11.0.
  - `test_plot_ecdf_rejects_a_bad_argument(kwargs, message)`: `palette=("black",)`, `palette=("black", "blu")`, `palette="black"`, `x_limits=(1.0, -1.0)`, `x_limits=(0.0, math.inf)`, `size=(0, 5)`, `title=5`, `result="x"`, and `groups={}`.
  - `test_plot_ecdf_rejects_a_group_without_a_value`: a group of genes that `de` does not have.
  - `test_plot_ecdf_rejects_a_de_without_log2_fold_change`: the message of `check_columns`.
  - `test_plot_ecdf_does_not_change_de`: `pd.testing.assert_frame_equal` with a copy from before the call.
- [ ] **Step 3: Run the tests.** Run `uv run pytest tests/plot tests/test_package.py -q`. Expected: FAIL with `ModuleNotFoundError: No module named 'seedmatch.plot'`.
- [ ] **Step 4: Write the code.** Write `__init__.py`, `_style.py`, and `ecdf.py`.
- [ ] **Step 5: Run the tests.** Expected: PASS. Also run `uv run ruff check`, `uv run ruff format --check`, and `uv run mypy`. Measure `plot_ecdf` and `savefig` to PNG at 100 dpi on 22,000 rows in two groups, and write the time in the docstring.
- [ ] **Step 6: Write the documents.**
  - `docs/r-map.md`, row `deseq_fc_ecdf`: remove "(the statistics and the groups, without the plot)" from the first cell. Add `seedmatch.plot.plot_ecdf` to the second cell. Replace the last sentence of the notes with these facts:
    - `$plot` becomes the `Figure` of `plot_ecdf(de, groups, result)`. `title`, `palette`, and `x_lims` become `title`, `palette`, and `x_limits`. `height` and `width` become `size`.
    - `output_filename` and `dpi` have no counterpart: call `Figure.savefig(path, dpi=...)`. R writes a PDF for each file name, also for `.png` (`R/deseq_fc_ecdf.R:357`).
    - `theme` has no counterpart, because R ignores it (`:349`). Name `D-063` and `D-064`.
  - `docs/r-map.md`, a new row `.plot_ecdf` → `seedmatch.plot.plot_ecdf`:
    - `res` with the column `Group` becomes `de` and `groups`. `test_results` becomes `result`, an `EcdfTestResult`.
    - `legend_labels` has no counterpart, because the legend counts the values of each line (`D-064`).
    - R gives `1e+05` for a statistic of 100000, and the port gives `100000`. Name `D-063`.
  - `docs/r-map.md`, a new row `SeedMatchR.theme`, `SeedMatchR.palette` → `seedmatch.plot.PALETTE`. `PALETTE` has the colors of R. The figure of `plot_ecdf` has the colors of the theme, and sizes near its sizes. There is no theme object.
  - `D-063`: "plot_ecdf labels each p-value."
    - R location: `.plot_ecdf` (`R/deseq_fc_ecdf.R:312-324`), called by `deseq_fc_ecdf` (`:133-145`).
    - R behavior: the branches compare `p < .01 && p >= .001`, `p > .01`, and `p < .001`. No branch takes p = 0.01, thus R stops with `object 'pvalue_label' not found`.
    - R behavior: a NaN or NA p-value stops R with `missing value where TRUE/FALSE needed`. A DTS p-value of 20 of 2000 permutations is exactly 0.01. Evidence: `tests/golden/plots/ecdf_labels.tsv`.
    - New behavior: the label of p = 0.01 is `0.01`, and NaN gives `NaN`. The other labels are the labels of R.
    - Reason: each result of `ecdf_test` must give a figure.
    - Test: `test_plot_ecdf_labels_each_p_value_D063` in `tests/plot/test_plot_ecdf.py`.
  - `D-064`: "plot_ecdf draws only the groups, and its legend counts the values of each line."
    - R location: `deseq_fc_ecdf` (`R/deseq_fc_ecdf.R:100-121`), and `.plot_ecdf` (`:326-336`).
    - R behavior: a row of `res` in no gene list has the group `NA`. `stat_ecdf` draws a gray line for it, with the legend key `NA`. Evidence: the case `no_group_r` of `ecdf_lines.tsv`.
    - R behavior: the count of a legend label is the length of the gene list (`:112`). Thus a gene that `res` does not have, a repeated gene, and a row without a fold change count. Evidence: the notes `legend absent_and_repeated` (`Background: 164`) and `legend na_fold_change`.
    - R behavior: the labels come in the order of `factor_order`, but the counts come in the order of `gene_lists` (`:112`, `:121`). Thus another order gives each group the count of another group. Evidence: the note `legend factor_order_reversed` (`mer7m8: 161`).
    - New behavior: the figure has one line for each key of `groups`, in the order of the keys. The count of a label is the number of finite `log2FoldChange` values of its line.
    - Reason: a legend label must describe its line.
    - Test: `test_plot_ecdf_draws_only_the_groups_D064`, `test_plot_ecdf_counts_the_values_of_each_line_D064`, and `test_plot_ecdf_follows_the_order_of_the_keys_D064`.
  - `docs/design.md`: the rows of `D-063` and `D-064`, and the bound of `plot_ecdf`. Checkpoint.

### Task 3: The Wasserstein figure

**Files:**
- Modify: `src/seedmatch/stats/wasserstein.py`, `src/seedmatch/stats/__init__.py`, `src/seedmatch/__init__.py`, `tests/stats/test_wasserstein.py`, `src/seedmatch/plot/__init__.py`, `tests/plot/test_plot_package.py`, `tests/helpers.py` (`run_core_pipeline`), `docs/r-map.md`, `docs/divergences.md` (`D-065`), `docs/design.md`
- Make: `src/seedmatch/stats/density.py`, `tests/stats/test_density.py`, `src/seedmatch/plot/wasserstein.py`, `tests/plot/test_plot_wasserstein.py`

**Interfaces:**
- Consumes: `check_wasserstein_sample`, `check_wasserstein_range`, `_grid`, `_ecdf`, `_Ecdf`, `_interpolate`, `wasserstein_areas`, `r_quantile`, `as_limit`, `check_at_least`, the names of `_style.py`, and the goldens `wasserstein_ecdf.tsv`, `wasserstein_density.tsv`, and `wasserstein_difference.tsv`.
- Produces: `as_wasserstein_sample`, `Ecdf`, `ecdf`, and the names of `seedmatch.stats.density` in "Shared interfaces", and:

```python
def wasserstein_difference(
    seed: ArrayLike, non_seed: ArrayLike, *, n_interp: int = 10000
) -> pd.DataFrame: ...

def plot_wasserstein(
    seed: ArrayLike,
    non_seed: ArrayLike,
    *,
    n_interp: int = 10000,
    size: tuple[float, float] = (15.0, 5.0),
) -> Figure: ...
```

Rules of `wasserstein_difference`:
- Rename `_as_wasserstein_sample` to `as_wasserstein_sample`, without a change of a message, because `seedmatch.plot.wasserstein` uses it. For the same reason, rename `_ecdf` to `ecdf` and `_Ecdf` to `Ecdf`.
- The checks are the checks of `wasserstein_areas`, with the same messages: `as_wasserstein_sample` of each sample, `check_wasserstein_range`, `as_limit("n_interp", ...)`, and `check_at_least("n_interp", ..., 2)`.
- The grid is `_grid(low, high, n_interp, 0, n_interp)`, from the lowest to the highest pooled value. `seed` and `non_seed` are `_interpolate` of the `ecdf` of each sample on the grid. `difference` is `seed - non_seed`.
- The result has the columns `x`, `seed`, `non_seed`, and `difference` (`float64`), with `n_interp` rows and a `RangeIndex`.
- The time is O(n log n + `n_interp` log n). The result holds 32 bytes for each grid point, thus `n_interp=10000` gives 320 kB.

Rules of `plot_wasserstein` (`D-065`):
- The samples go through `as_wasserstein_sample`, and `wasserstein_difference(seed, non_seed, n_interp=n_interp)` makes the other checks. `size` goes through `as_size`.
- The figure has 3 Axes in one row, in the order of `dens.plot`, `cdf.plot`, and `diff.plot`. Each Axes has the x limits `expanded(low, high)` of the pooled values, and `style_gray_axes`.
- Density Axes: the grid is `np.linspace(low, high, DENSITY_POINTS)`. For the non-seed sample and then the seed sample, `gaussian_density` with `nrd0_bandwidth` of the sample gives the curve.
  - `fill_between` fills each curve with `#F8766D` (non-seed) or `#00BFC4` (seed), alpha 0.25, and the labels `non-seed` and `seed`.
  - A black line of the width `linewidth_points(0.5)` draws each curve. The labels are `log2(FC)` and `Density`, and the legend title is `group`.
- ECDF Axes: let (vlo, vhi) be the x limits, u the `knots` of `ecdf` of a sample, and F its `heights`.
  - The line has the points `[vlo, u[0], *u, vhi]` and `[0, 0, *F, 1]` (`D-065`).
  - The non-seed line is black, and the seed line is red, each with the width `linewidth_points(1.25)`. The labels are `log2(FC)` and `eCDF`, and the legend title is `group`.
- Difference Axes: one black line through `x` and `difference` of `wasserstein_difference`, with the width `linewidth_points(0.5)`. The labels are `log2(FC)` and `Seed - NonSeed eCDF`.
- The title of the figure is the title of `diff.plot`: `f"Wass(- censored) = {neg}   Wass(+ censored) = {pos}   Wass (total) = {total}"`. Each value is `pasted` of a field of `wasserstein_areas(seed, non_seed, n_interp=n_interp)`. The title is too long for one Axes of the figure, thus `fig.suptitle` holds it.
- The time is O(512 · n) for the density of n values, and O(`n_interp` + n log n) for the other Axes. The work memory of the density is 4 MiB. The figure holds 512 points for each curve and `n_interp` points for the difference.

The expected results for the seed and the non-seed fold changes of the test data, with the default `n_interp`:
- The x limits of each Axes are (-8.743981512594443, 4.002092292512264), within 1e-12.
- The seed ECDF line has the first points (-8.743981512594443, 0), (-8.16461452145323, 0), and (-8.16461452145323, 0.05555555555555555). Its last point is (4.002092292512264, 1.0).
- The fill colors are `#f8766d40` and `#00bfc440` with `to_hex(..., keep_alpha=True)`.
- The title is `Wass(- censored) = 2.196   Wass(+ censored) = 0   Wass (total) = 2.196`. With `n_interp=100`, the areas are 2.197.

- [ ] **Step 1: Write the tests that fail.** In `tests/stats/test_wasserstein.py`:
  - `test_wasserstein_difference_matches_r(case)`: parametrized over `test_data` and `wass_40x300_ties` with `n_interp=100`. `x` and `difference` agree with `wasserstein_difference.tsv` with `np.testing.assert_allclose(..., rtol=1e-9, atol=1e-12)`. The x values of R differ in the last bit, thus the test does not compare with `==`.
  - `test_wasserstein_difference_gives_the_areas`: for the test data, the trapezoid sum of `abs(difference)` over `x` differs from `wasserstein_areas(...).absolute` by at most 5e-4.
  - `test_wasserstein_difference_has_the_grid_of_the_areas`: the columns and the dtypes above, `n_interp` rows, `x[0]` the lowest value, `x[-1]` the highest value, and `difference == seed - non_seed`.
  - `test_wasserstein_difference_rejects_what_the_areas_reject(kwargs, message)`: 9 values, one distinct value, NaN, and `n_interp=1`, each with the message of `wasserstein_areas`.
- In `tests/stats/test_density.py`:
  - `test_nrd0_bandwidth_matches_r(case, group)`: the bandwidth notes of `wasserstein_density.tsv`, within a relative tolerance of 1e-15.
  - `test_gaussian_density_is_the_exact_sum`: the values `[0.0, 1.0, 3.0]`, the grid `[-1.0, 0.0, 2.0]`, and the bandwidth 0.5 give the sum of the rule, within 1e-15. With `_BLOCK` monkeypatched to 2, the result is the same.
- In `tests/plot/test_plot_wasserstein.py`:
  - `test_plot_wasserstein_density_agrees_with_r(case)`: for each group, the black line of the density Axes has the x values of the golden, within `rtol=1e-12` and `atol=1e-12`. Its y values differ from the golden by at most 0.005 times the maximum of the golden.
  - `test_plot_wasserstein_draws_the_ecdf_of_the_areas_D065(case)`: the points from the third to the last but one equal the finite rows of `wasserstein_ecdf.tsv`, within a relative tolerance of 1e-12. The first two points have y 0, and the last point has y 1. The golden rows of each group start with x `-Inf`.
  - `test_plot_wasserstein_draws_the_difference`: the line of the difference Axes equals `x` and `difference` of `wasserstein_difference(...)`.
  - `test_plot_wasserstein_gives_the_r_titles_and_labels`: the expected results above, the labels, the legend texts `non-seed` and `seed`, and the colors. The title notes of `wasserstein_density.tsv` give the R titles, and its note `panel x range test_data` gives the x limits within 1e-12.
  - `test_plot_wasserstein_rejects_bad_input(kwargs, message)`: the cases of `wasserstein_areas`, and `size=(0, 5)`.
  - Add `plot_wasserstein` to the script of `test_plot_functions_load_no_pyplot` and to `test_plot_functions_keep_the_rc_params`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/stats/test_wasserstein.py tests/stats/test_density.py tests/plot -q`. Expected: FAIL with `ImportError` for `wasserstein_difference`, for `seedmatch.stats.density`, and for `plot_wasserstein`.
- [ ] **Step 3: Write the code.** Export `wasserstein_difference` from `seedmatch.stats` and from `seedmatch`, and `plot_wasserstein` from `seedmatch.plot`. Add `wasserstein_difference` of the seed and the non-seed fold changes of `sirna_de()` with `n_interp=100` to `run_core_pipeline`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass. Measure `plot_wasserstein` on 20,000 and 2,000 normal values with `n_interp=10000`, without and with `savefig` to PNG at 100 dpi. Give the time and the peak of `tracemalloc` in the docstring.
- [ ] **Step 5: Write the documents.**
  - `docs/r-map.md`, row `wass_dist`: add `seedmatch.wasserstein_difference` and `seedmatch.plot.plot_wasserstein` to the second cell. Replace "The plots `cdf.plot`, `dens.plot`, and `diff.plot` have no counterpart." with these facts:
    - `plot_wasserstein(seed, non_seed, n_interp=...)` draws `dens.plot`, `cdf.plot`, and `diff.plot` as the Axes of one figure, in this order. The title of `diff.plot` is the title of the figure. `wasserstein_difference` gives the data of `diff.plot`.
    - The density is the exact Gaussian sum with the bandwidth of `bw.nrd0`. `density()` of R uses an FFT approximation, and for the test data the difference is at most 0.0014 of the maximum. Name `D-065`.
  - `D-065`: "plot_wasserstein draws the ECDF that the areas use."
    - R location: `wass_dist` (`R/wass_dist.R:89-93`).
    - R behavior: `stat_ecdf(geom = "line")` adds a point at `-Inf` with 0 and a point at `Inf` with 1. `geom_line` joins them to the values. Evidence: `tests/golden/plots/wasserstein_ecdf.tsv`.
    - R behavior: thus the line goes from 0 at the left edge of the panel to the height of the first value, with a slope. The slope changes with the width of the panel. The areas use 0 below the first value (`yleft = 0`, `:70-71`).
    - New behavior: each line is the interpolated ECDF of `wasserstein_areas`. It is 0 up to the first value, and then it goes up to the height of the first value. Straight lines join the values, and the line is 1 after the last value.
    - Reason: the figure must show the functions that the areas integrate.
    - Test: `test_plot_wasserstein_draws_the_ecdf_of_the_areas_D065` in `tests/plot/test_plot_wasserstein.py`.
  - `docs/design.md`: the row of `D-065`, and the bounds of `wasserstein_difference` and `plot_wasserstein`. Checkpoint.

### Task 4: The k-mer histogram

**Files:**
- Modify: `src/seedmatch/kmers.py`, `src/seedmatch/__init__.py`, `tests/test_kmers.py`, `src/seedmatch/plot/__init__.py`, `tests/plot/test_plot_package.py`, `tests/helpers.py` (`run_core_pipeline`), `docs/r-map.md`, `docs/divergences.md` (`D-010`, `D-066`), `docs/design.md`
- Make: `src/seedmatch/plot/kmers.py`, `tests/plot/test_plot_kmers.py`

**Interfaces:**
- Consumes: `kmer_gene_counts`, `as_limit`, `check_at_least`, `check_type`, `as_text`, the names of `_style.py`, and the golden `kmer_hist.tsv`.
- Produces:

```python
def kmer_histogram(gene_counts: pd.Series, *, max_bins: int = 50) -> pd.DataFrame: ...

def plot_kmer_histogram(
    gene_counts: pd.Series,
    target: str,
    *,
    title: str = "K-mer gene frequencies",
    max_bins: int = 50,
    size: tuple[float, float] = (7.0, 5.0),
) -> Figure: ...
```

Rules of `kmer_histogram` (`D-066`):
- `gene_counts` goes through `check_type("gene_counts", gene_counts, pd.Series)`. A dtype that is not an integer dtype raises `gene_counts has dtype float64; expected an integer dtype, as kmer_gene_counts gives`.
- A missing value raises `gene_counts has a missing value for k-mer 'C'; expected a count`. The k-mer is the index label of the first missing value. A nullable `Int64` Series can hold it.
- A Series without a value raises `gene_counts has no value; expected at least one k-mer`. A negative value raises `gene_counts holds -1 for k-mer 'AAA'; expected 0 or more`. The k-mer is the index label of the first negative value.
- `max_bins` goes through `as_limit("max_bins", ...)` and `check_at_least("max_bins", ..., 1)`.
- Let lo and hi be the lowest and the highest count. The width w is the first value of 1, 2, 5, 10, 20, 50, 100, and so on, with `hi // w - lo // w + 1 <= max_bins`. w is a Python `int`.
- Bin b, from `lo // w` to `hi // w`, holds the counts from `b * w` to `b * w + w - 1`. Each bin holds w integer counts, also the first bin.
- The result has the columns `low`, `high`, and `kmers` (`int64`), one row for each bin, and a `RangeIndex`. A bin without a k-mer has the row with 0.
- The time is O(n) for n k-mers, with one `np.bincount`. The memory is O(n + `max_bins`).

Rules of `plot_kmer_histogram` (`D-066`):
- `kmer_histogram(gene_counts, max_bins=max_bins)` gives the bins. `title` and `target` go through `as_text`, and `size` through `as_size`.
- A target that is not in the index of `gene_counts` raises `target 'CTCTATAA' is not a k-mer of gene_counts; expected a k-mer of 7 letters, for example 'AAAAAAA'`. The example is the first index label, and the length is its length.
- The figure has one Axes. `Axes.bar` draws one bar for each bin, with the left edge `low - 0.5`, the width `high - low + 1`, and the height `kmers`. The bars have the face `#E5E5E5` (`grey90`), a black edge, and the edge width 0.75.
- `axvline` draws the target line at `gene_counts[target]`, with the color `firebrick` and the width 2.25.
- The title is `title`, bold, with the size 14.4. The axis labels are `Number of genes with k-mer` and `Frequency`. The Axes has no top and no right edge line.
- The time is the time of `kmer_histogram`, and the figure holds one bar for each bin.

The expected results for `kmer_gene_counts(kmer_counts(sirna_seqs()))`:
- `kmer_histogram` gives 44 rows. The first rows are (0, 1, 2142) and (2, 3, 2242), and the last row is (86, 87, 1). The sum of `kmers` is 16384.
- `max_bins=10` gives the width 10 and 9 rows, with the counts 8484, 6021, 1570, 225, 59, ..., 4, 1, and 1. `max_bins=3` gives (0, 49, 16359) and (50, 99, 25).
- For the counts `[7, 12, 12, 30]` with the index `["A", "C", "G", "T"]`, `max_bins=50` gives 24 rows from (7, 7, 1) to (30, 30, 1). `max_bins=5` gives (0, 9, 1), (10, 19, 2), (20, 29, 0), and (30, 39, 1).
- `plot_kmer_histogram(..., "CTCTATA")` gives 44 bars. The first bar has the x -0.5, the width 2, and the height 2142. The target line is at 18.

- [ ] **Step 1: Write the tests that fail.** In `tests/test_kmers.py`:
  - `test_kmer_histogram_bins_the_integer_counts_D066`: the 44 rows above. The note `hist counts` of `kmer_hist.tsv` starts with 3394, the first bar of R.
  - `test_kmer_histogram_sums_the_counts_of_r`: each bin holds the sum of `kmers` of `kmer_hist.tsv` for the values from `low` to `high`.
  - `test_kmer_histogram_widens_the_bins_to_max_bins`: the results for `max_bins=10` and `max_bins=3`.
  - `test_kmer_histogram_starts_at_the_bin_of_the_lowest_count`: the results for `[7, 12, 12, 30]`.
  - `test_kmer_histogram_rejects_bad_input(kwargs, message)`: a list, a `float64` Series, an empty `int64` Series, a value -1, `max_bins=0`, and `max_bins=True`.
- In `tests/plot/test_plot_kmers.py`:
  - `test_plot_kmer_histogram_draws_the_bins_and_the_target_D066`: the results above. The figure has one Axes and the title `K-mer gene frequencies`. The notes `pages` and `title calls` give 2 pages and the second title of R.
  - `test_plot_kmer_histogram_rejects_a_target_that_is_not_a_kmer_D066(target)`: `CUCUAUA`, `ctctata`, `CTCTATAA`, and `NOTAKMER`. The notes give `NA` for each target in R.
  - `test_plot_kmer_histogram_rejects_a_bad_argument(kwargs, message)`: `title=5`, `target=5`, `size=(0, 5)`, and `max_bins=0`.
  - Add `plot_kmer_histogram` to the two tests of `tests/plot/test_plot_package.py`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/test_kmers.py tests/plot -q`. Expected: FAIL with `ImportError` for `kmer_histogram` and for `plot_kmer_histogram`.
- [ ] **Step 3: Write the code.** Export `kmer_histogram` from `seedmatch`, and `plot_kmer_histogram` from `seedmatch.plot`. Add `kmer_histogram(kmer_gene_counts(kmer_counts(seqs, k=3)))` to `run_core_pipeline`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass. Measure `kmer_histogram` on the 16,384 counts of k = 7, and on 4^12 counts of k = 12. Give the times in the docstring.
- [ ] **Step 5: Write the documents.**
  - `D-010`: replace "The plot module of a later plan draws the histogram." with "`plot_kmer_histogram` draws the histogram (`D-066`)." In the table of `docs/design.md`, replace "The plot module makes the histogram." with "`plot_kmer_histogram` draws the histogram."
  - `docs/r-map.md`, row `kmer_hist`: add `seedmatch.kmer_histogram` and `seedmatch.plot.plot_kmer_histogram` to the second cell. Replace "`title` and the plot have no counterpart in this version, and a later version draws the histogram." with these facts:
    - `plot_kmer_histogram(gene_counts, target.sequence, title=...)` draws the histogram, and `kmer_histogram` gives its bins (`D-066`).
    - R gives `lw = 3` to `abline`, and R reads it as `lwd`. The line of the port has the width 2.25 points, the width of the R line.
  - `D-066`: "plot_kmer_histogram draws one figure with bins of integer counts and the target line."
    - R location: `kmer_hist` (`R/kmer_hist.R:30-34`).
    - R behavior: `plot(graphics::hist(...))` draws the histogram two times. Page 1 has the title and the label of `kmer_hist`. Page 2 has `Histogram of gene.counts` and `gene.counts`, and `abline` draws the target line on page 2 only. Evidence: the notes `pages` and `title calls` of `tests/golden/plots/kmer_hist.tsv`.
    - R behavior: `breaks = 50` gives the breaks of `pretty`, and each bin is closed on the right. The first bin also holds its lower break. Thus for the breaks 0, 2, ..., 88, the first bin holds the counts 0, 1, and 2, and each other bin holds two counts.
    - R behavior: on `Test_seqs.fa`, the first bar of R is 3394, but only 2142 k-mers have 0 or 1 genes. For other counts, `pretty` gives more than 50 bins, or bins of a fractional width. For example, the range 1 to 7 gives 60 bins of width 0.1.
    - R behavior: a target that is not a k-mer of the counts gives `NA`, for example a lowercase target or a target of another length. R then draws no line and gives no message. Evidence: the target notes of `kmer_hist.tsv`.
    - New behavior: `kmer_histogram` gives bins of w integer counts each. w is the first of 1, 2, 5, 10, 20, 50, and so on, that gives at most `max_bins` bins.
    - New behavior: `plot_kmer_histogram` draws one Axes with the bars, the title, and the target line. A target that is not a k-mer raises `InputError`.
    - Reason: each bar must hold the same number of integer counts, and the target line must be on the histogram.
    - Test: `test_kmer_histogram_bins_the_integer_counts_D066` in `tests/test_kmers.py`. `test_plot_kmer_histogram_draws_the_bins_and_the_target_D066` and `test_plot_kmer_histogram_rejects_a_target_that_is_not_a_kmer_D066` in `tests/plot/test_plot_kmers.py`.
  - `docs/design.md`: the row of `D-066`, and the bounds of `kmer_histogram` and `plot_kmer_histogram`. Checkpoint.

### Task 5: The alignment figure

**Files:**
- Modify: `src/seedmatch/plot/__init__.py`, `tests/plot/test_plot_package.py`, `docs/r-map.md`, `docs/divergences.md` (`D-067`), `docs/design.md`
- Make: `src/seedmatch/plot/alignment.py`, `tests/plot/test_plot_alignment.py`

**Interfaces:**
- Consumes: `alignment_text`, `as_text`, and the golden `alignment_plot.tsv`.
- Produces:

```python
def plot_alignment(text: str) -> Figure: ...
```

Rules of `plot_alignment` (`D-067`):
- `text` goes through `as_text("text", text)`. Text that does not have 3 lines raises `text has 2 line(s); expected the 3 lines of alignment_text`.
- The label has the three lines in order, each after its name and one space. The names are `siRNA `, six spaces, and `target`. A line end comes between two lines. `paste` of `R/plotAlignment.R:63` joins the lines in the same way, but with `miRNA ` (`D-067`).
- Let W be the length of the longest line of the label. The figure has the size (0.1 · W + 0.4, 0.9) inches, a white face, and no Axes. This is the `theme_void` of R.
- `Figure.text` puts the label at (0.5, 0.5) with `ha="center"`, `va="center"`, and `multialignment="left"`. The text has the family `monospace`, the weight `bold`, and the size 11, the default size of `geom_text`.
- The time and the memory are O(W).

- [ ] **Step 1: Write the tests that fail** in `tests/plot/test_plot_alignment.py`:
  - `test_plot_alignment_shows_the_r_lines_D067`: the text of `alignment_text` for the first `mer7m8` match of `sirna_seqs()`. The figure has no Axes and one text. Line 1 of the text is line 1 of `alignment_plot.tsv` with `siRNA ` in place of `miRNA `. Lines 2 and 3 are equal to the golden. The note `data.frame row names` gives `[siRNA ]`.
  - `test_plot_alignment_uses_a_bold_monospace_font`: the family, the weight, and the size of the text.
  - `test_plot_alignment_fits_the_text`: with `left_pad=20` and `right_pad=20`, the window extent of the text from `FigureCanvasAgg(fig).get_renderer()` is inside the box of the figure.
  - `test_plot_alignment_rejects_bad_text(text, message)`: 5, and `"a\nb"`.
  - Add `plot_alignment` to the two tests of `tests/plot/test_plot_package.py`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/plot -q`. Expected: FAIL with `ImportError` for `plot_alignment`.
- [ ] **Step 3: Write the code.** Export `plot_alignment` from `seedmatch.plot`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass. If the text of `test_plot_alignment_fits_the_text` is not inside the figure, make the factor 0.1 of the width larger. Keep the font size 11, and give the extent and the factor in the report of the task.
- [ ] **Step 5: Write the documents.**
  - `docs/r-map.md`, a new row `plotAlignment` with `output.type = "ggplot"` → `seedmatch.plot.plot_alignment`. The function takes the text of `alignment_text`, thus the arguments of R go to `alignment_text`. `aln.name` has no counterpart. Name `D-045`, `D-046`, and `D-067`.
  - `docs/r-map.md`, row `plotAlignment` with `output.type = "data.frame"`: replace "`output.type = "ggplot"` has no counterpart in this version." with "`seedmatch.plot.plot_alignment` draws the figure of `output.type = "ggplot"`."
  - `D-067`: "plot_alignment names the guide line siRNA."
    - R location: `plotAlignment` (`R/plotAlignment.R:63`), and `.create_alignment_df` (`:135-137`).
    - R behavior: the figure names the guide line `miRNA`, but the data.frame names it `siRNA `. Evidence: `tests/golden/plots/alignment_plot.tsv` and its note `data.frame row names`.
    - New behavior: the names are `siRNA `, six spaces, and `target`, as in the data.frame.
    - Reason: the guide of SeedMatchR is the guide strand of an siRNA. One thing must have one name.
    - Test: `test_plot_alignment_shows_the_r_lines_D067` in `tests/plot/test_plot_alignment.py`.
  - `docs/design.md`: the row of `D-067`, and the bound of `plot_alignment`. Checkpoint.

### Task 6: The range figure

**Files:**
- Modify: `src/seedmatch/plot/__init__.py`, `tests/plot/test_plot_package.py`, `docs/r-map.md`, `docs/divergences.md` (`D-068`), `docs/design.md`
- Make: `src/seedmatch/plot/ranges.py`, `tests/plot/test_plot_ranges.py`

**Interfaces:**
- Consumes: `check_columns`, `check_integer_columns`, `as_text`, `as_number`, the names of `_style.py`, `find_matches`, and the golden `ranges_plot.tsv`.
- Produces:

```python
# Each keyword-only argument is an option of plotIRanges or of plotRangesDF.
def plot_ranges(  # noqa: PLR0913
    ranges: pd.DataFrame,
    *,
    x_limits: tuple[float, float] | None = None,
    title: str = "Ranges",
    color: str = "blue",
    bar_height: float = 0.25,
    size: tuple[float, float] = (7.0, 5.0),
) -> Figure: ...
```

Rules of `plot_ranges` (`D-068`):
- `ranges` goes through `check_columns(ranges, ("start", "end"), name="ranges", label="ranges table")`. Then it goes through `check_integer_columns(ranges, ["start", "end"], label="ranges table")`.
- Thus a result of `find_matches` is a valid table. Another column has no effect.
- A table without a row raises `ranges table has no row; expected at least one range`. For an empty `IRanges`, R stops with `need finite 'xlim' values`.
- A row with `end < start` raises `ranges table has start 9 and end 7 in row 2; expected end >= start`. The row is the position from 0 of the first such row.
- `x_limits` is `None`, or it goes through `as_limits("x_limits", x_limits)`. `title` goes through `as_text`, `color` through `as_color("color", color)`, and `size` through `as_size`.
- `bar_height` goes through `as_number`. A value that is not above 0 and at most 1 raises `bar_height 0.0 is out of range; expected a number above 0 and at most 1`. NaN gives the same message.
- The limit 1 of `bar_height` is the distance between two rows. A higher bar covers the bar of the next row.
- Let n be the number of rows. Row i, from 1 in the order of the table, gets a bar from `start - 0.5` to `end + 0.5`. The bar goes from `i - bar_height / 2` to `i + bar_height / 2` in y.
- Thus a range of one base has the width 1. The row labels of the table have no effect.
- One `PolyCollection` holds the bars, with the face color and the edge color `color`, and the line width 0.
- The x limits are `x_limits`, or else `expanded(min(start) - 0.5, max(end) + 0.5)`. The y limits are (0, n + 1), as in `plotIRanges`.
- The axis labels are `Position` and `Ranges`, and the title of the Axes is `title`. `style_gray_axes` styles the Axes, as the default theme of `plotRangesDF`.
- The function does not change `ranges`. It has no row limit, because a limit drops rows without a message, as `head()` of R does.
- The time and the memory are O(n). One `PolyCollection` holds all bars. It holds 4 vertices for each row, and it adds no artist for a row.
- The planner measured 42 µs and 451 bytes for each row with `savefig` to PNG, on 10^5 rows. The docstring and "Performance and bounds" give the time and the memory for each row.

The expected results for the table `eight` of "Facts from the planner probes":
- The figure has one Axes with one `PolyCollection` of 8 paths. The first path has the vertices (9.5, 0.875), (16.5, 0.875), (16.5, 1.125), and (9.5, 1.125).
- The x limits are (2.15, 163.85) within 1e-12, and the y limits are (0, 9). The face color is `#0000ff`.
- The title is `Ranges`, and the axis labels are `Position` and `Ranges`.

- [ ] **Step 1: Write the tests that fail** in `tests/plot/test_plot_ranges.py`:
  - `test_plot_ranges_draws_each_row_D068`: the expected results above. The case `eight_rows` of `ranges_plot.tsv` has 6 rows. For each golden row, xmin is the left edge of the bar of the port plus 0.5. ymin is the low edge of the bar.
  - `test_plot_ranges_uses_the_options_D068`: `x_limits=(0.0, 1000.0)`, `title="My title"`, `color="red"`, and `bar_height=0.8` give the limits (0, 1000) and the title `My title`.
  - In the same test, the face color is `#ff0000`, and the bar of row 1 goes from 0.6 to 1.4 in y. The note `options` gives no title and no coordinate limits in R.
  - `test_plot_ranges_puts_each_row_at_its_position_D068`: `eight.iloc[[2, 6, 7]]` gives the y centers 1, 2, and 3. The same rows with the index `["a", "b", "c"]` give the same centers.
  - In the same test, the case `subset` gives ymin 2.875, 6.875, and 7.875 in R, and the case `text_row_names` gives NA.
  - `test_plot_ranges_shows_a_range_of_one_base_D068`: start 20 and end 20 give a bar from 19.5 to 20.5. The case `one_base` gives xmin 20 and xmax 20, and the note `plotIRanges one_base warning` gives the warning of R.
  - `test_plot_ranges_draws_a_find_matches_result`: `find_matches(Guide(GUIDE).target(Seed.MER7M8), sirna_seqs())` gives 20 rows and 20 bars. The first bar goes from 72.5 to 79.5, for the range 73 to 79 of `ENSRNOG00000016275`.
  - `test_plot_ranges_rejects_bad_input(kwargs, message)`: a list, a table without `end`, and a `float64` column `start`. Also a table without a row (`eight.iloc[:0]`), and a row with start 9 and end 7.
  - More cases of the same test: `x_limits=(1.0, -1.0)`, `title=5`, `color="blu"`, `bar_height=0.0`, `bar_height=1.5`, `bar_height=math.nan`, `bar_height=True`, and `size=(0, 5)`.
  - `test_plot_ranges_does_not_change_ranges`: `pd.testing.assert_frame_equal` with a copy from before the call.
  - Add `plot_ranges` to the two tests of `tests/plot/test_plot_package.py`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/plot -q`. Expected: FAIL with `ImportError` for `plot_ranges`.
- [ ] **Step 3: Write the code.** Export `plot_ranges` from `seedmatch.plot`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass. Measure `plot_ranges` on 10^5 ranges, without and with `savefig` to PNG at 100 dpi. Divide the time and the peak of `tracemalloc` by the number of rows. Give the totals and the values for each row in the docstring.
- [ ] **Step 5: Write the documents.**
  - `docs/r-map.md`, a new row `plotIRanges`, `plotRangesDF` → `seedmatch.plot.plot_ranges`:
    - `ir` becomes `ranges`, a table with the integer columns `start` and `end`, for example a result of `find_matches`. A data.frame stops `plotIRanges` with `invalid time series parameters specified (1)`, because it takes only an `IRanges`.
    - `xlim`, `main`, and `col` become `x_limits`, `title`, and `color`. The default title is `Ranges`, not `IRanges Plot` or `Ranges DF Plot`.
    - `bar_height` keeps its name. R accepts each value of `bar_height`, but the port accepts only a value above 0 and at most 1. A higher bar covers the bar of the next row.
    - `...` has no counterpart. `NAMESPACE` does not export the R functions. Name `D-068`.
  - `D-068`: "plot_ranges draws each range at the position of its row, with the options."
    - R location: `plotIRanges` (`R/plotIRanges.R:3-14`) and `plotRangesDF` (`:16-21`).
    - R behavior: `plotRangesDF` draws only the first 6 rows, because of `head()` (`:20`). It ignores `xlim`, `main`, and `col`. Evidence: the case `eight_rows` and the note `options` of `tests/golden/plots/ranges_plot.tsv`.
    - R behavior: the y value of a row is `as.integer(rownames(ir))` (`:17-18`). Thus a subset keeps the old positions of its rows, and text row names give NA. Evidence: the cases `subset` and `text_row_names`.
    - R behavior: the bar of a range of one base has the width 0, because `geom_rect` goes from `start` to `end`. `plotIRanges` skips such a range with a warning. Evidence: the case `one_base` and the note `plotIRanges one_base warning`.
    - New behavior: `plot_ranges` draws one bar for each row, at the position of the row in the table. A bar goes from `start - 0.5` to `end + 0.5`, and each option has an effect.
    - Reason: the figure must show each range, and a range of one base covers one position.
    - Test: `test_plot_ranges_draws_each_row_D068`, `test_plot_ranges_uses_the_options_D068`, `test_plot_ranges_puts_each_row_at_its_position_D068`, and `test_plot_ranges_shows_a_range_of_one_base_D068` in `tests/plot/test_plot_ranges.py`.
  - `docs/design.md`: the row of `D-068`. In "Performance and bounds", the bound of `plot_ranges` gives O(n), the measured time and memory for each row, and the absence of a row limit. Checkpoint.

### Task 7: The PR and ROC figure

**Files:**
- Modify: `src/seedmatch/plot/__init__.py`, `tests/plot/test_plot_package.py`, `docs/r-map.md`, `docs/divergences.md` (`D-069`), `docs/design.md`
- Make: `src/seedmatch/plot/benchmark.py`, `tests/plot/test_plot_benchmark.py`

**Interfaces:**
- Consumes: `Benchmark`, `benchmark_predictions`, `pr_curve`, and `roc_curve` of `seedmatch.stats`, `check_type`, `as_size`, `new_figure`, and the golden `benchmark_plot.tsv`.
- Produces:

```python
def plot_benchmark(
    benchmark: Benchmark, *, size: tuple[float, float] = (12.0, 5.0)
) -> Figure: ...
```

Rules of `plot_benchmark` (`D-069`):
- `benchmark` goes through `check_type("benchmark", benchmark, Benchmark)`, and `size` through `as_size`. Thus a `PrCurve` alone raises the message of `check_type`.
- The figure has the PR Axes and the ROC Axes in one row, in this order, with a color bar to the right of each. `fig.axes` holds the PR Axes, the ROC Axes, and then the color bars.
- The points of the PR Axes are `recall` and `precision` of `benchmark.pr.curve`, in the order of the rows. The points of the ROC Axes are `fpr` and `tpr` of `benchmark.roc.curve`.
- Each Axes has one `LineCollection` with a segment from each row to the next row. The line width is 2.25 points, the width of `lwd = 3`.
- `_SCALE` holds 100 colors: `hsv_to_rgb` of the hue `0.8 * j / 99`, for j from 0 to 99, with the saturation 1 and the value 1. This is `scale.color` of PRROC.
- Let t be the threshold of a row, and lo and hi the lowest and the highest finite threshold of the curve.
- If hi > lo, the position of a row is `np.clip((99 * (t - lo)) / (hi - lo), 0, 99)`, with the operations in this order. The color index from 0 is `np.rint(1 + position) - 1`.
- Thus an infinite threshold gets the first or the last color. A segment has the color of its first row, as in R.
- If the curve does not have two distinct finite thresholds, each segment gets `_SCALE[49]`, `#00ff60`. Let c be the finite threshold, or 0. The color bar then goes from c - 0.5 to c + 0.5.
- The color bar is `fig.colorbar(ScalarMappable(norm=Normalize(lo, hi), cmap=ListedColormap(_SCALE)), ax=axes)`. It shows the 100 colors, as `image` of R does.
- The view of each Axes is (-0.04, 1.04) for x and for y. This is the default `xaxs = "r"` of R for the limits 0 and 1.
- The titles are `f"PR curve\nAUC = {benchmark.pr.auc_integral:.7g}"` and `f"ROC curve\nAUC = {benchmark.roc.auc:.7g}"`. `.7g` gives the 7 significant digits of `format` of R.
- The PR area is `auc_integral`, the default `auc.type` of `plot.PRROC`. For many negatives, it differs from PRROC (`D-050`).
- The axis labels are `Recall` and `Precision`, and `FPR` and `Sensitivity`. The Axes keep the default style of matplotlib, because the plots of R use base graphics.
- The time and the memory are O(R) for R rows of the two curves. The PR curve has at most 2K + 2 + max(100, P) rows for K distinct scores and P positives.
- The planner measured 45 µs and 393 bytes for each curve row with `savefig`, on 262,604 rows. The figure has no row limit, as `plot_ranges` has none.

The expected results for `plot_benchmark(benchmark_predictions(sirna_de()))`. The test data has 63 positives and 116 negatives:
- `fig.axes` has 4 Axes. The PR Axes has a `LineCollection` of 105 segments, and the ROC Axes has one of 3 segments.
- The first PR segment goes from (1.0, 63/179) to the same point, with the color `#ff0000`. The last PR segment has the color `#cc00ff`. The PR colors are `#ff0000`, `#00ff60`, and `#cc00ff`.
- The ROC segments go from (0, 0) to (2/116, 0) in `#cc00ff`, then to (4/116, 14/63) in `#00ff60`, and then to (1, 1) in `#ff0000`.
- The titles are `PR curve\nAUC = 0.4723177` and `ROC curve\nAUC = 0.591954`. The line width is 2.25, and the limits are (-0.04, 1.04).
- Each color bar goes from 0 to 2.

- [ ] **Step 1: Write the tests that fail** in `tests/plot/test_plot_benchmark.py`:
  - `test_plot_benchmark_draws_the_r_segments`: the expected results above. The segments of each curve equal the case `benchmark_1` of `benchmark_plot.tsv` with a relative tolerance of 1e-12. The colors equal the golden colors in lowercase.
  - `test_plot_benchmark_gives_the_r_titles_and_labels`: the titles, the axis labels, the line width, and the limits. The notes `title PR`, `title ROC`, `labels PR`, `labels ROC`, `lwd`, and `usr` give the R values.
  - In the same test, the note `title ROC unpatched` gives `AUC = 0.3762679`. R gives the labels to `scores.class1` (`D-006`), thus the AUC of R is not the AUC of the port.
  - `test_plot_benchmark_has_a_color_bar_for_each_curve`: after `fig.canvas.draw()`, the face colors of the `QuadMesh` of each color bar Axes are 100 colors from `#ff0000` to `#cc00ff`, and the range is 0 to 2. The notes `scale` and `color bar PR` give the same colors and range.
  - `test_plot_benchmark_colors_one_distinct_score_D069(score)`: parametrized over 0.0, -3.0, and 2.0. The input is a `Benchmark` of `pr_curve` and `roc_curve` of `np.full(20, score)` with the labels `[1, 0] * 10`.
  - In the same test, each segment is `#00ff60`, and the color bar goes from score - 0.5 to score + 0.5. The notes `one score 0` and `one score -3` give the R errors.
  - In the same test, the case `one_score_2` gives `#00FF85`, and the note `one_score_2 color bar` gives 1.8 to 2.18.
  - `test_plot_benchmark_colors_an_infinite_score_D069`: the scores `[inf] + [1.0] * 9 + [0.0] * 10` with the labels `[1, 0] * 10`. The thresholds 0, 1, and `inf` give `#ff0000`, `#cc00ff`, and `#cc00ff`.
  - In the same test, the color bar goes from 0 to 1. The note `an Inf score` gives the R error.
  - `test_plot_benchmark_rejects_bad_input(kwargs, message)`: `benchmark="x"`, `benchmark=benchmark.pr`, and `size=(0, 5)`.
  - Add `plot_benchmark` to the two tests of `tests/plot/test_plot_package.py`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/plot -q`. Expected: FAIL with `ImportError` for `plot_benchmark`.
- [ ] **Step 3: Write the code.** Export `plot_benchmark` from `seedmatch.plot`.
- [ ] **Step 4: Run the tests.** Expected: PASS, and ruff and mypy pass.
- Measure `plot_benchmark` on 10^5 distinct scores, and on 10^6 scores of 8 values with 30% positives. Measure without and with `savefig` to PNG at 100 dpi.
- Give the times, the peak of `tracemalloc`, and the values for each curve row in the docstring.
- [ ] **Step 5: Write the documents.**
  - `docs/r-map.md`, row `benchmark_preds`: add `seedmatch.plot.plot_benchmark` to the second cell. Add these facts to the notes:
    - R documents the result as two figures (`R/benchmark.R:11-13`). But it returns the curve objects of PRROC (`:32-40`), and it draws nothing.
    - In R, `plot(res$PRC)` and `plot(res$ROC)` draw the figures. `plot_benchmark(benchmark)` draws the two curves in one figure.
  - `docs/r-map.md`, a new row `plot.PRROC` (PRROC 1.4, for the curves of `benchmark_preds`) → `seedmatch.plot.plot_benchmark`:
    - The function takes a `Benchmark`, not one curve object. The figure has the defaults of `plot.PRROC`, and the other arguments of `plot.PRROC` have no counterpart.
    - For some AUC values, R gives a short scientific form, for example `1e-04` for 0.0001. The port gives `0.0001`. Name `D-050` and `D-069`.
  - `D-069`: "plot_benchmark colors the curve of one distinct score and of an infinite score."
    - R location: `plot.PRROC` and `getColor` of PRROC 1.4, for the curves of `benchmark_preds` (`R/benchmark.R:32-40`).
    - R behavior: for one distinct score s, the color range goes from 0.9 s to s + 0.09 s. For s = 0, each color is NA.
    - R behavior: for s = 0 and for s < 0, `image` stops R with `increasing 'x' and 'y' values expected`. Evidence: the notes `one score 0` and `one score -3` of `tests/golden/plots/benchmark_plot.tsv`.
    - R behavior: an infinite score gives the color NA to its threshold, and R stops with `'to' must be a finite number`. Evidence: the note `an Inf score`.
    - New behavior: the color range comes from the finite thresholds, and an infinite threshold gets the color of the nearest end of the scale.
    - New behavior: without two distinct finite thresholds, each segment gets the middle color `#00ff60`, and the color bar goes from c - 0.5 to c + 0.5. For s = 2, R gives `#00FF85` (the case `one_score_2`).
    - Reason: each result of `benchmark_predictions` must give a figure. A seed count of 0 for each gene is a valid input.
    - Test: `test_plot_benchmark_colors_one_distinct_score_D069` and `test_plot_benchmark_colors_an_infinite_score_D069` in `tests/plot/test_plot_benchmark.py`.
  - `docs/design.md`: the row of `D-069`. In "Performance and bounds", the bound of `plot_benchmark` gives O(R) and the measured time and memory for each curve row. Checkpoint.

### Task 8: The HTML table of the report

**Files:**
- Make: `src/seedmatch/report/html.py`, `tests/report/test_report_html.py`
- Modify: `src/seedmatch/report/__init__.py`, `src/seedmatch/__init__.py`, `tests/helpers.py` (`run_core_pipeline`), `docs/r-map.md`, `docs/divergences.md` (`D-070`), `docs/design.md`

**Interfaces:**
- Consumes: `SeedMatchReport`, `REPORT_CATEGORIES`, `Guide`, `Seed`, `seed_window`, `check_type`, `check_columns`, `check_numeric_columns`, and the golden `report_table.tsv`.
- Produces:

```python
def report_html(report: SeedMatchReport, guide: Guide) -> str: ...
```

Rules of `report_html` (`D-070`):
- `report` goes through `check_type("report", report, SeedMatchReport)`, and `guide` through `check_type("guide", guide, Guide)`.
- `report.report` goes through `check_columns(report.report, ("category", *REPORT_CATEGORIES, "Total", "group"), name="report.report", label="report table")`. Then `check_numeric_columns` of the categories and `Total`.
- The columns of the table come in the order of `gt`. First come `Full-D0` to `Full-D4`, `18mer-D0` to `18mer-D4`, and `15mer-D0` to `15mer-D4`. Then come `8mer`, `7mer-m8`, `7mer-A1`, `6mer`, and `Total`. With the stub, the table has 21 columns.
- The result is one `<style>` element and one `<table class="seedmatch-report">`. The table has no `id`. Each CSS selector starts with `table.seedmatch-report`, and the text color is black.
- `<thead>` has these rows:
  - the title `In-silico siRNA Binding Prediction`, with `colspan="21"` and the class `title`
  - the subtitle `Identifying siRNA hits in the transcriptome`, with `colspan="21"` and the class `subtitle`
  - an empty stub cell with `rowspan="2"`, three spanners with `colspan="5"` and `scope="colgroup"`, and the cells `8mer`, `7mer-m8`, `7mer-A1`, `6mer`, and `Total` with `rowspan="2"`
  - the labels `D0` to `D4` of each spanner, with `scope="col"`
- The spanner labels come from `seed_window(seed, len(guide.sequence))` (`D-070`). `Seed.FULL` gives `Full Guide Strand (g2:gN)`, `Seed.MER18` gives `18-mer (g2:g19)`, and `Seed.MER15` gives `15-mer (g2:g16)`.
- The docstring says that `guide` must be the guide of the call of `seed_match_report`.
- Each group of `report.report["group"]`, in the order of its first row, has a `<tbody>`. The first row of the `<tbody>` is the group name, with `colspan="21"`, `scope="colgroup"`, and the class `group`.
- Each row of the group follows, with the category in `<th scope="row" class="stub">` and one `<td>` for each column. The rows 2 and 4 of a group have the class `striped`, as in `gt`.
- A cell of the row `% off-target` has the text `f"{value:.2f}%"`. A cell of another row has the text `f"{value:,.0f}"`. These are the formats of `fmt_percent` and `fmt_number`.
- A category cell of the row `% off-target` has `style="background-color: #FF{c:02X}{c:02X}"`, but the `Total` cell has no style. The value c is `round(2.55 * (100 - v))`, with v clipped to 0 to 100.
- A NaN cell has the text `NaN`, as `fmt_percent` of gt gives, and no style. Evidence: `tmp/claude/plan-fix/nan_pct.R` of the controller.
- `html.escape` escapes each text. Then `.encode("ascii", "xmlcharrefreplace").decode("ascii")` makes each character outside ASCII a character reference. Thus the result is ASCII.
- The time and the memory are O(R) for R rows of `report.report`.

The expected results for `report_html(seed_match_report(sirna_de(), sirna_seqs(), Guide(GUIDE)), Guide(GUIDE))`:
- The cell texts equal the column `text` of `report_table.tsv` for each row and column. For example, `In silico predictions` gives `125` for `Total`, and `% off-target` gives `45.83%` for `18mer-D4` and `33.60%` for `Total`.
- The spanners are `Full Guide Strand (g2:g23)`, `18-mer (g2:g19)`, and `15-mer (g2:g16)`. The group row is `SeedMatchReport`.
- The backgrounds of the row `% off-target`: 0 gives `#FFFFFF`, and 100 gives `#FF0000`. These are the colors of `gt`.
- More backgrounds: 33.333333 gives `#FFAAAA`, 45.833333 gives `#FF8A8A`, 28.767123 gives `#FFB6B6`, 22.222222 gives `#FFC6C6`, and 75 gives `#FF4040`. `gt` gives `#FFBFAA`, `#FFA68B`, `#FFC8B5`, `#FFD5C6`, and `#FF6846`.

- [ ] **Step 1: Write the tests that fail** in `tests/report/test_report_html.py`. A class of the test file reads the HTML with `html.parser.HTMLParser`, and it keeps the tag, the attributes, and the text of each cell.
  - `test_report_html_gives_the_r_cell_texts`: the cell texts and the column order above. The note `columns` gives the same order.
  - `test_report_html_gives_the_r_headers`: the title, the subtitle, the group row, the labels `D0` to `D4`, and `colspan="21"`. The notes `title`, `subtitle`, `group`, `colspan`, `spanner Full`, and `spanner 18mer` give the same texts.
  - `test_report_html_labels_the_windows_from_the_guide_D070`: the spanner `15-mer (g2:g16)`, and the note `spanner 15mer` gives `15-mer (g2:g19)`. A report of `Guide(GUIDE[:21])` gives `Full Guide Strand (g2:g21)`. The note `spanner Full of a 21-base guide` gives `(g2:g23)`.
  - `test_report_html_has_no_fixed_id_D070`: the HTML has no `id=`, the table has the class `seedmatch-report`, and each CSS selector starts with `table.seedmatch-report`. The note `id` gives `mygt`.
  - `test_report_html_colors_the_percent_row`: the backgrounds above, and no style on the `Total` cell or on a count cell.
  - `test_report_html_stripes_rows_2_and_4`: rows 2 and 4 of each group have the class `striped`, and rows 1 and 3 have no class. The note `striped rows` gives `2,4`.
  - `test_report_html_puts_a_comma_in_a_large_count`: a copy of the report frame with `Total` 12345 in the row `In silico predictions` gives `12,345`. The note `fmt_number 12345` gives the same text.
  - `test_report_html_escapes_text`: the group `<b>&é` gives `&lt;b&gt;&amp;&#233;`, and the result is ASCII.
  - `test_report_html_rejects_bad_input(args, message)`: the report `"x"`, the guide `"UUA"`, and a report frame without `Total`.
- [ ] **Step 2: Run the tests.** Run `uv run pytest tests/report/test_report_html.py -q`. Expected: FAIL with `ImportError` for `report_html`.
- [ ] **Step 3: Write the code.** Export `report_html` from `seedmatch.report` and from `seedmatch`. Add `report_html` of the report of `run_core_pipeline` to `run_core_pipeline`.
- [ ] **Step 4: Run the tests.** Run `uv run pytest tests/report tests/test_offline.py tests/test_package.py -q`. Expected: PASS, and ruff and mypy pass.
- [ ] **Step 5: Write the documents.**
  - `docs/design.md`, "Not identical by design": replace "The report table uses a pandas `Styler` in HTML, not `gt`." with "The report table is the HTML of `report_html`, not `gt`. pandas 2.3.3 gives a `Styler` only with jinja2 3.1.2 or later, and the store lock has jinja2 3.0.3."
  - `docs/r-map.md`, row `SeedMatchReport`: replace "The `gt` table (`table`) has no counterpart." with "`report_html(report, guide)` gives the HTML of the `gt` table (`table`) (`D-070`)."
  - `docs/r-map.md`, a new row `.make_table` → `seedmatch.report_html`:
    - R returns a `gt_tbl` object, and the port returns HTML text. The columns, the labels, the formats, the title, and the subtitle are the same.
    - The color of a percent cell is a linear blend from white to red, with black text. `gt` blends in another color space and gives white text on a dark cell. 0 and 100 give the colors of `gt`. Name `D-070`.
  - `D-070`: "report_html labels the windows from the guide, and gives no fixed element ID."
    - R location: `.make_table` (`R/SeedMatchReport.R:142-196`), called by `SeedMatchReport` (`:135`).
    - R behavior: the spanner labels are fixed texts: `Full Guide Strand (g2:g23)` (`:158`), `18-mer (g2:g19)` (`:162`), and `15-mer (g2:g19)` (`:166`).
    - R behavior: the `15mer` window is g2 to g16 (`R/get_seed.R:121-123`), thus the 15-mer label is wrong. For a guide that does not have 23 bases, the `Full` label is wrong. Evidence: the spanner notes of `tests/golden/plots/report_table.tsv`.
    - R behavior: `gt::gt(..., id = "mygt")` (`:145`) gives each table the element ID `mygt`, and each CSS rule selects `#mygt`. Thus two tables in one HTML page have the same ID and the same rules.
    - New behavior: `seed_window` of `Seed.FULL`, `Seed.MER18`, and `Seed.MER15` for the guide gives the labels. The table has the class `seedmatch-report` and no ID, and each CSS rule selects this class.
    - Reason: a label must give the window of the search, and an ID must be unique in a page.
    - Test: `test_report_html_labels_the_windows_from_the_guide_D070` and `test_report_html_has_no_fixed_id_D070` in `tests/report/test_report_html.py`.
  - `docs/design.md`: the row of `D-070`, and the bound of `report_html`. Checkpoint.

### Task 9: The README and the design text

**Files:**
- Modify: `README.md`, `tests/test_readme.py`, `docs/design.md`, `docs/plans/2026-11-20-network.md`, `docs/plans/2026-11-20-small-analyses.md`

**Interfaces:**
- Consumes: the public names of Tasks 2 to 8, and the golden row `readme` of `tests/golden/stats/ecdf.tsv`.

- [ ] **Step 1: Write the README part.** Add a section "## Plots and the HTML report" after "Other analyses". It has a paragraph and a Python block. The paragraph gives these facts:
  - The functions of `seedmatch.plot` work only with the extra `plot`.
  - A plot function returns a `matplotlib.figure.Figure` and writes no file. `savefig` writes the file.
  - `report_html` is in the core, and its result is ASCII HTML.
  - Name `D-063`, `D-064`, `D-066`, `D-070`, `D-068`, and `D-069`.
- The block uses the names of the examples above:

```python
from pathlib import Path

from seedmatch import report_html
from seedmatch.plot import (
    plot_alignment,
    plot_benchmark,
    plot_ecdf,
    plot_kmer_histogram,
    plot_ranges,
    plot_wasserstein,
)

# The ECDF of the two gene sets of the ecdf_test example, with the result of the test.
ecdf_figure = plot_ecdf(de, {"Background": background, "mer7m8": seed_genes}, ecdf_result)
ecdf_figure.savefig("ecdf.png", dpi=300)

# The density, the ECDF, and the ECDF difference of the Wasserstein areas.
wasserstein_figure = plot_wasserstein(
    de.loc[de["mer7m8"] >= 1, "log2FoldChange"], de.loc[de["mer7m8"] == 0, "log2FoldChange"]
)

# The number of genes with each 7-mer, with a line at the target.
kmer_figure = plot_kmer_histogram(genes_with_kmer, target.sequence)

# The alignment text of the example above.
alignment_figure = plot_alignment(alignment)

# One bar for each match range of the first example.
ranges_figure = plot_ranges(ranges)

# The PR curve and the ROC curve of the benchmark example.
benchmark_figure = plot_benchmark(benchmark)

# The SeedMatchReport table as HTML.
Path("report.html").write_text(report_html(report, guide), encoding="utf-8")
```

- [ ] **Step 2: Write the README test.** Add `test_readme_plot_example_draws_the_figures` to `tests/test_readme.py`:
  - The legend of `ecdf_figure` is `Background: 121` and `mer7m8: 17`. The left title is `Statistic: 0.46\np-value: 0.001`, from the row `readme`, KS, `greater` of `ecdf.tsv`.
  - `ecdf.png` starts with the PNG signature `b"\x89PNG\r\n\x1a\n"`.
  - The title of `wasserstein_figure` starts with `Wass(- censored) = 2.033`, the Wass statistic of the row `readme`.
  - The target line of `kmer_figure` is at `genes_with_kmer[target.sequence]`.
  - The `PolyCollection` of `ranges_figure` has one path for each row of `ranges`.
  - The title of the PR Axes of `benchmark_figure` is `f"PR curve\nAUC = {benchmark.pr.auc_integral:.7g}"`.
  - `report.html` holds `In-silico siRNA Binding Prediction`.
- [ ] **Step 3: Change the other parts of the README.**
  - "Installation": add the command of the extra, `uv pip install "seedmatch[plot] @ git+https://github.com/inflexa-ai/seedmatch.git"`.
  - "Scope": add "the plots of the ECDF, the Wasserstein areas, the k-mer counts, the alignment, the match ranges, and the PR and ROC curves" to the first list. Also add "the HTML table of the report". Remove the sentence "These parts of SeedMatchR are not in this version:" and its list.
  - "Differences from SeedMatchR": add one item. `plot_ecdf` gives a figure for each p-value, and its legend counts the values of each line. R stops for a p-value of 0.01, and its legend counts the gene lists (D-063, D-064).
- [ ] **Step 4: Change the design text.**
  - "Packages": the row of `seedmatch.report` becomes "Full search, best match, report, and the HTML table of the report." The row of `seedmatch.plot` becomes "Figures from results: the ECDF, the Wasserstein areas, the k-mer histogram, the alignment, the ranges, and the PR and ROC curves. Extra `plot`."
  - "Rules": after the rule of the plot functions, add "A plot function uses the object API of matplotlib. It does not import `matplotlib.pyplot`, and it does not change `matplotlib.rcParams`."
  - "Tests": add a part "### Plots" before "### Documentation examples". Text: "The tests of `seedmatch.plot` read the data of the artists: the points of each line, the bars, the texts, and the colors. They do not compare pixels. `tools/golden/R/plots.R` gives the layer data of ggplot2 and the cells of `gt`."
- [ ] **Step 5: Close the plot lines of the earlier plans.**
  - In `docs/plans/2026-11-20-network.md`, "Later plans", replace "- The plot of a network (delivery step 9)." with "- The plot of a network: closed. SeedMatchR has no function that draws a network, thus there is no port (`docs/plans/2026-11-21-plot-and-report.md`)."
  - In `docs/plans/2026-11-20-small-analyses.md`, "Later plans", replace "- The k-mer histogram, the `ggplot` output of `plotAlignment`, and the PR and ROC plots (delivery step 9)." with this line:
  - "- The k-mer histogram, the `ggplot` output of `plotAlignment`, and the PR and ROC plots: closed. `docs/plans/2026-11-21-plot-and-report.md` ports them."
- [ ] **Step 6: Run the tests.** Run `uv run pytest tests/test_readme.py -q`. Expected: PASS. Checkpoint.

### Task 10: Verification of the plan

- [ ] **Step 1: Run the checks one time.** Run `uv run ruff check`, `uv run ruff format --check`, `uv run mypy`, and `PYTHONPATH=tools/coverage NUMBA_DISABLE_JIT=1 uv run pytest -p serial_pool --cov`. Expected: each passes, and the coverage is not below 95%. Give a reason for each line of the new core code without coverage.
- [ ] **Step 2: Compare the documents with the code.**
  - Each public name of `seedmatch.plot`, and each new core function, has a row in `docs/r-map.md`.
  - Each `D-NNN` of this plan has an entry, a table row, and a test.
  - Each README example runs in `tests/test_readme.py`.
  - `tests/test_offline.py` passes with the new calls in `run_core_pipeline`, and `import seedmatch` loads no module of matplotlib.
  - `git status --short` shows no `Rplots.pdf` and no figure file outside `tmp/`.
  - `grep -rn "plot of a network\|PR and ROC plots" docs README.md` finds only the closed lines of the network plan and of the small-analyses plan, and each closed line names this plan.
  - Checkpoint.
