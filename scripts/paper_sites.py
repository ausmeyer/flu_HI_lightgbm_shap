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


WIC2023_H3_SITE_ROWS = [
    {"site": 53, "paper_category": "top25_site_C", "paper_note": "Top-25 posterior inclusion position; antigenic site C in Fig 2 discussion."},
    {"site": 121, "paper_category": "top25_site_D", "paper_note": "Top-25 posterior inclusion position; antigenic site D in Fig 2 discussion."},
    {"site": 126, "paper_category": "top25_site_A", "paper_note": "Top-25 posterior inclusion position; antigenic site A in Fig 2 discussion."},
    {"site": 131, "paper_category": "high_confidence_site_A", "paper_note": "Posterior inclusion probability at least 0.95; antigenic site A."},
    {"site": 135, "paper_category": "high_confidence_site_A", "paper_note": "Posterior inclusion probability at least 0.95; antigenic site A."},
    {"site": 137, "paper_category": "top25_site_A", "paper_note": "Top-25 posterior inclusion position; antigenic site A in Fig 2 discussion."},
    {"site": 144, "paper_category": "top25_site_A", "paper_note": "Top-25 posterior inclusion position; antigenic site A in Fig 2 discussion."},
    {"site": 145, "paper_category": "high_confidence_site_A", "paper_note": "Posterior inclusion probability at least 0.95; antigenic site A."},
    {"site": 157, "paper_category": "high_confidence_site_B", "paper_note": "Posterior inclusion probability at least 0.95; antigenic site B."},
    {"site": 158, "paper_category": "top25_site_B", "paper_note": "Top-25 posterior inclusion position; antigenic site B in Fig 2 discussion."},
    {"site": 159, "paper_category": "top25_site_B", "paper_note": "Top-25 posterior inclusion position; antigenic site B in Fig 2 discussion."},
    {"site": 160, "paper_category": "top25_site_B", "paper_note": "Top-25 posterior inclusion position; antigenic site B in Fig 2 discussion."},
    {"site": 173, "paper_category": "top25_site_D", "paper_note": "Top-25 posterior inclusion position; antigenic site D in Fig 2 discussion."},
    {"site": 189, "paper_category": "high_confidence_site_B", "paper_note": "Posterior inclusion probability at least 0.95; antigenic site B."},
    {"site": 193, "paper_category": "high_confidence_site_B", "paper_note": "Posterior inclusion probability at least 0.95; antigenic site B."},
]


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
