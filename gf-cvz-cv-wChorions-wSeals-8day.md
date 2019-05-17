Zebrafish behavior
================

#### Comparing Germ-free, Conventionalized, and Conventional Zebrafish in Sealed 96-well Plates

### Day 1 - HMAT

![Red lines indicate the window of movement we will analyze
statistically](gf-cvz-cv-wChorions-wSeals-8day_files/figure-gfm/hmat-all-movement-plot-1.png)

![Post-hoc pairwise comparisons conducted with Dunn’s Kruskal-Wallis
Multiple Comparisons Test with False Discovery
Rate](gf-cvz-cv-wChorions-wSeals-8day_files/figure-gfm/hmat-auc-stats-1.png)

![AUCS broken down by both plate and
treatment](gf-cvz-cv-wChorions-wSeals-8day_files/figure-gfm/hmat-plot-with-plates-1.png)

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  AUC by interaction(Treatment, plate.id)
    ## Kruskal-Wallis chi-squared = 40.028, df = 5, p-value = 1.474e-07

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">

<caption>

Significant (P.adj \< 0.05) pairwise comparisons

</caption>

<thead>

<tr>

<th style="text-align:left;">

Comparison

</th>

<th style="text-align:right;">

Z

</th>

<th style="text-align:right;">

P.unadj

</th>

<th style="text-align:right;">

P.adj

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

CV.1 - CV.2

</td>

<td style="text-align:right;">

3.578148

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

0.002

</td>

</tr>

<tr>

<td style="text-align:left;">

CV.1 - CVZ.1

</td>

<td style="text-align:right;">

4.982573

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

0.000

</td>

</tr>

<tr>

<td style="text-align:left;">

CV.1 - CVZ.2

</td>

<td style="text-align:right;">

2.253341

</td>

<td style="text-align:right;">

0.024

</td>

<td style="text-align:right;">

0.045

</td>

</tr>

<tr>

<td style="text-align:left;">

CVZ.1 - CVZ.2

</td>

<td style="text-align:right;">

\-2.730676

</td>

<td style="text-align:right;">

0.006

</td>

<td style="text-align:right;">

0.014

</td>

</tr>

<tr>

<td style="text-align:left;">

CVZ.1 - GF.1

</td>

<td style="text-align:right;">

\-3.423794

</td>

<td style="text-align:right;">

0.001

</td>

<td style="text-align:right;">

0.002

</td>

</tr>

<tr>

<td style="text-align:left;">

CV.1 - GF.2

</td>

<td style="text-align:right;">

5.043540

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

0.000

</td>

</tr>

<tr>

<td style="text-align:left;">

CVZ.2 - GF.2

</td>

<td style="text-align:right;">

2.803811

</td>

<td style="text-align:right;">

0.005

</td>

<td style="text-align:right;">

0.013

</td>

</tr>

<tr>

<td style="text-align:left;">

GF.1 - GF.2

</td>

<td style="text-align:right;">

3.493212

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

0.002

</td>

</tr>

</tbody>

</table>

### Day 5 - LPR

#### Second Epoch Only

![LPR results from epoch 2 by Treatment. Yellow line indicates light
cycle (vs dark
cycle)](gf-cvz-cv-wChorions-wSeals-8day_files/figure-gfm/lpr-ep2-by-treat-plot-1.png)

<table class="table table-striped table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Comparing GF and CVZ to CV (Kolmogorov-Smirnov Tests)

</caption>

<thead>

<tr>

<th style="text-align:left;">

Treatment

</th>

<th style="text-align:right;">

AUC

</th>

<th style="text-align:right;">

Pval

</th>

<th style="text-align:right;">

RelativeRatio

</th>

<th style="text-align:left;">

significance p\<0.01

</th>

<th style="text-align:left;">

activity

</th>

<th style="text-align:left;">

Sig

</th>

</tr>

</thead>

<tbody>

<tr grouplength="3">

<td colspan="7" style="border-bottom: 1px solid;">

<strong>Light Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CV (Ctrl)

</td>

<td style="text-align:right;">

63.296

</td>

<td style="text-align:right;">

1.000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ

</td>

<td style="text-align:right;">

33.649

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

\-0.4683875

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    ">HYPO</span>

</td>

<td style="text-align:left;">

<span style="     color: green !important;">YES</span>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:right;">

64.448

</td>

<td style="text-align:right;">

0.357

</td>

<td style="text-align:right;">

0.0182102

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

<tr grouplength="3">

<td colspan="7" style="border-bottom: 1px solid;">

<strong>Dark Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CV (Ctrl)

</td>

<td style="text-align:right;">

96.368

</td>

<td style="text-align:right;">

1.000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ

</td>

<td style="text-align:right;">

90.380

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

\-0.0621372

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:right;">

93.664

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

\-0.0280591

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

</tbody>

</table>

<table class="table table-striped table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Comparing GF and CVZ to CV (Kolmogorov-Smirnov Tests)

</caption>

<thead>

<tr>

<th style="text-align:left;">

Treatment

</th>

<th style="text-align:right;">

AUC

</th>

<th style="text-align:right;">

Pval

</th>

<th style="text-align:right;">

RelativeRatio

</th>

<th style="text-align:left;">

significance p\<0.01

</th>

<th style="text-align:left;">

activity

</th>

<th style="text-align:left;">

Sig

</th>

</tr>

</thead>

<tbody>

<tr grouplength="2">

<td colspan="7" style="border-bottom: 1px solid;">

<strong>Light Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ (Ctrl)

</td>

<td style="text-align:right;">

33.649

</td>

<td style="text-align:right;">

1.000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:right;">

64.448

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

0.9153241

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    ">HYPER</span>

</td>

<td style="text-align:left;">

<span style="     color: green !important;">YES</span>

</td>

</tr>

<tr grouplength="2">

<td colspan="7" style="border-bottom: 1px solid;">

<strong>Dark Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ (Ctrl)

</td>

<td style="text-align:right;">

90.380

</td>

<td style="text-align:right;">

1.000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:right;">

93.664

</td>

<td style="text-align:right;">

0.001

</td>

<td style="text-align:right;">

0.0363360

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

NA

</td>

<td style="text-align:left;">

<span style="     color: red !important;">NO</span>

</td>

</tr>

</tbody>

</table>

![LPR results from epoch 2 by Treatment and Plate. Yellow line indicates
light cycle (vs dark
cycle)](gf-cvz-cv-wChorions-wSeals-8day_files/figure-gfm/lpr-ep2-by-plate-plot-1.png)

#### All\* Epochs

\* 2nd thru 4th. The first is dropped in the script given to my by Lisa.
![Lines indicate treatment means, shading indicates the bootstrapped 95%
CI for each mean. Yellow segments at bottom indicate light intervals (vs
dark
intervals).](gf-cvz-cv-wChorions-wSeals-8day_files/figure-gfm/lpr-all-eps-plot-1.png)

<table class="table table-striped table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Comparing GF and CVZ to CV (Kolmogorov-Smirnov Tests)

</caption>

<thead>

<tr>

<th style="text-align:left;">

Treatment

</th>

<th style="text-align:left;">

cycle

</th>

<th style="text-align:right;">

AUC

</th>

<th style="text-align:right;">

PercDiff

</th>

<th style="text-align:right;">

Pval

</th>

<th style="text-align:left;">

significance p\<0.01

</th>

<th style="text-align:left;">

significance p\<0.01 and 50% cutoff

</th>

<th style="text-align:left;">

activity

</th>

</tr>

</thead>

<tbody>

<tr grouplength="3">

<td colspan="8" style="border-bottom: 1px solid;">

<strong>Light Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CV (Ctrl)

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

26690.21

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

1.000

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

10036.44

</td>

<td style="text-align:right;">

\-62.397

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    ">HYPO</span>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

24917.96

</td>

<td style="text-align:right;">

\-6.640

</td>

<td style="text-align:right;">

0.008

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr grouplength="3">

<td colspan="8" style="border-bottom: 1px solid;">

<strong>Dark Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CV (Ctrl)

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

174850.61

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

1.000

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

111408.11

</td>

<td style="text-align:right;">

\-36.284

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

172678.94

</td>

<td style="text-align:right;">

\-1.242

</td>

<td style="text-align:right;">

0.227

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

</tr>

</tbody>

</table>

<table class="table table-striped table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Comparing GF to CVZ (Kolmogorov-Smirnov Tests)

</caption>

<thead>

<tr>

<th style="text-align:left;">

Treatment

</th>

<th style="text-align:left;">

cycle

</th>

<th style="text-align:right;">

AUC

</th>

<th style="text-align:right;">

PercDiff

</th>

<th style="text-align:right;">

Pval

</th>

<th style="text-align:left;">

significance p\<0.01

</th>

<th style="text-align:left;">

significance p\<0.01 and 50% cutoff

</th>

<th style="text-align:left;">

activity

</th>

</tr>

</thead>

<tbody>

<tr grouplength="2">

<td colspan="8" style="border-bottom: 1px solid;">

<strong>Light Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ (Ctrl)

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

10036.44

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

24917.96

</td>

<td style="text-align:right;">

148.275

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    ">HYPER</span>

</td>

</tr>

<tr grouplength="2">

<td colspan="8" style="border-bottom: 1px solid;">

<strong>Dark Interval</strong>

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

CVZ (Ctrl)

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

111408.11

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NO

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left; padding-left: 2em;" indentlevel="1">

GF

</td>

<td style="text-align:left;">

Cycle

</td>

<td style="text-align:right;">

172678.94

</td>

<td style="text-align:right;">

54.997

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

YES

</td>

<td style="text-align:left;">

<span style=" font-weight: bold;    ">HYPER</span>

</td>

</tr>

</tbody>

</table>

![Lines indicate plate means, shading indicates the bootstrapped 95% CI
for each mean. Yellow segments at bottom indicate light intervals (vs
dark
intervals).](gf-cvz-cv-wChorions-wSeals-8day_files/figure-gfm/lpr-all-eps-plot-by-plate-1.png)

### Day 5 - LSR
