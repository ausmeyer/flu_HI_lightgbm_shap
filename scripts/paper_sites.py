#!/usr/bin/env python3
"""Reference HA1 site lists used for paper-comparison outputs."""

from __future__ import annotations


NEHER2016_H3_SITE_ROWS = [
    {"site": 62, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
    {"site": 121, "paper_category": "named_substitution", "paper_note": "FU02 cluster change attributed partly to N121T"},
    {"site": 135, "paper_category": "named_substitution", "paper_note": "Repeated effects include K135E; SI87 to BE89 includes G135N"},
    {"site": 140, "paper_category": "named_substitution", "paper_note": "Repeated effects include K140E"},
    {"site": 144, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
    {"site": 145, "paper_category": "koel7_and_cluster_transition", "paper_note": "Koel 7 site; SI87 to BE89 involves N145K"},
    {"site": 155, "paper_category": "koel7", "paper_note": "Koel 7 site listed in paper"},
    {"site": 156, "paper_category": "koel7_and_cluster_transition", "paper_note": "Koel 7 site; WU95 to SY97 set includes K156Q; FU02 change includes Q156H"},
    {"site": 158, "paper_category": "koel7_and_named_substitution", "paper_note": "Koel 7 site; repeated effects include K158R; WU95 to SY97 set includes E158K"},
    {"site": 159, "paper_category": "koel7_and_named_substitution", "paper_note": "Koel 7 site; repeated effects include Y159F; later text discusses position 159"},
    {"site": 186, "paper_category": "cluster_transition", "paper_note": "SI87 to BE89 includes I186S"},
    {"site": 189, "paper_category": "koel7_and_named_substitution", "paper_note": "Koel 7 site; repeated effects include K189N; text notes S189N can be small"},
    {"site": 193, "paper_category": "koel7_and_cluster_transition", "paper_note": "Koel 7 site; SI87 to BE89 includes N193S"},
    {"site": 196, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
    {"site": 276, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
]


# Primary Harvey comparison: all HA1 positions with posterior inclusion probability
# (PIP) >= 0.95 in the structurally-aware model (Harvey et al. 2023, Results;
# PLoS Comput Biol 19:e1010885, p. 15). The paper reports n = 14:
# site A 131, 135, 144, 145; site B 157, 158, 159, 189, 193; site D 173;
# receptor-binding site 194, 225; near the receptor-binding site 138, 223.
HARVEY2023_PIP95_SITE_ROWS = [
    {"site": 131, "paper_category": "structurally_aware_pip_ge_0.95_site_A", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site A."},
    {"site": 135, "paper_category": "structurally_aware_pip_ge_0.95_site_A", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site A."},
    {"site": 138, "paper_category": "structurally_aware_pip_ge_0.95_near_rbs", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; located close to the receptor-binding site."},
    {"site": 144, "paper_category": "structurally_aware_pip_ge_0.95_site_A", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site A."},
    {"site": 145, "paper_category": "structurally_aware_pip_ge_0.95_site_A", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site A."},
    {"site": 157, "paper_category": "structurally_aware_pip_ge_0.95_site_B", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site B."},
    {"site": 158, "paper_category": "structurally_aware_pip_ge_0.95_site_B", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site B."},
    {"site": 159, "paper_category": "structurally_aware_pip_ge_0.95_site_B", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site B."},
    {"site": 173, "paper_category": "structurally_aware_pip_ge_0.95_site_D", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site D."},
    {"site": 189, "paper_category": "structurally_aware_pip_ge_0.95_site_B", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site B."},
    {"site": 193, "paper_category": "structurally_aware_pip_ge_0.95_site_B", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; antigenic site B."},
    {"site": 194, "paper_category": "structurally_aware_pip_ge_0.95_rbs", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; defined as belonging to the receptor-binding site."},
    {"site": 223, "paper_category": "structurally_aware_pip_ge_0.95_near_rbs", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; located close to the receptor-binding site."},
    {"site": 225, "paper_category": "structurally_aware_pip_ge_0.95_rbs", "paper_note": "Harvey et al. 2023 structurally-aware model; PIP at least 0.95; defined as belonging to the receptor-binding site."},
]

# Restricted 15-site subset used in an earlier draft for figure/table comparability.
# Mixes naive-model PIP>=0.95 sites with top-25 posterior-inclusion sites in canonical
# antigenic regions. It is NOT Harvey's complete identification and omits 225.
# Keep only for labeled sensitivity/footnote comparisons.
HARVEY2023_RESTRICTED15_SITE_ROWS = [
    {"site": 53, "paper_category": "restricted15_top25_site_C", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site C."},
    {"site": 121, "paper_category": "restricted15_top25_site_D", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site D."},
    {"site": 126, "paper_category": "restricted15_top25_site_A", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site A."},
    {"site": 131, "paper_category": "restricted15_pip_ge_0.95_site_A", "paper_note": "Restricted 15-site subset: PIP at least 0.95; antigenic site A."},
    {"site": 135, "paper_category": "restricted15_pip_ge_0.95_site_A", "paper_note": "Restricted 15-site subset: PIP at least 0.95; antigenic site A."},
    {"site": 137, "paper_category": "restricted15_top25_site_A", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site A."},
    {"site": 144, "paper_category": "restricted15_top25_site_A", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site A."},
    {"site": 145, "paper_category": "restricted15_pip_ge_0.95_site_A", "paper_note": "Restricted 15-site subset: PIP at least 0.95; antigenic site A."},
    {"site": 157, "paper_category": "restricted15_pip_ge_0.95_site_B", "paper_note": "Restricted 15-site subset: PIP at least 0.95; antigenic site B."},
    {"site": 158, "paper_category": "restricted15_top25_site_B", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site B."},
    {"site": 159, "paper_category": "restricted15_top25_site_B", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site B."},
    {"site": 160, "paper_category": "restricted15_top25_site_B", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site B."},
    {"site": 173, "paper_category": "restricted15_top25_site_D", "paper_note": "Restricted 15-site subset (not the primary Harvey rule): top-25 posterior inclusion; antigenic site D."},
    {"site": 189, "paper_category": "restricted15_pip_ge_0.95_site_B", "paper_note": "Restricted 15-site subset: PIP at least 0.95; antigenic site B."},
    {"site": 193, "paper_category": "restricted15_pip_ge_0.95_site_B", "paper_note": "Restricted 15-site subset: PIP at least 0.95; antigenic site B."},
]

# Backward-compatible alias: primary Harvey set used by comparison scripts.
WIC2023_H3_SITE_ROWS = HARVEY2023_PIP95_SITE_ROWS


SHAH2024_H3_SITE_ROWS = [
    {"site": 45, "paper_category": "aggregated_top20_epitope_C", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope C."},
    {"site": 62, "paper_category": "aggregated_top20_epitope_E", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope E."},
    {"site": 122, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 131, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 135, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 138, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 140, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 142, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 144, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 145, "paper_category": "aggregated_top20_epitope_A", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope A."},
    {"site": 158, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 159, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 171, "paper_category": "aggregated_top20_epitope_D", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope D."},
    {"site": 173, "paper_category": "aggregated_top20_epitope_D", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope D."},
    {"site": 183, "paper_category": "aggregated_top20_unknown", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; outside established epitope set."},
    {"site": 186, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 189, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 193, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 194, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 196, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 197, "paper_category": "aggregated_top20_epitope_B", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope B."},
    {"site": 208, "paper_category": "aggregated_top20_epitope_D", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope D."},
    {"site": 213, "paper_category": "aggregated_top20_epitope_D", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope D."},
    {"site": 219, "paper_category": "aggregated_top20_epitope_D", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope D."},
    {"site": 223, "paper_category": "aggregated_top20_unknown", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; outside established epitope set."},
    {"site": 225, "paper_category": "aggregated_top20_unknown", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; outside established epitope set."},
    {"site": 241, "paper_category": "aggregated_top20_unknown", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; outside established epitope set."},
    {"site": 261, "paper_category": "aggregated_top20_epitope_E", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope E."},
    {"site": 269, "paper_category": "aggregated_top20_unknown", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; outside established epitope set."},
    {"site": 311, "paper_category": "aggregated_top20_epitope_C", "paper_note": "Aggregated across seasonal top-20 sites in Fig. 4a; epitope C."},
]
