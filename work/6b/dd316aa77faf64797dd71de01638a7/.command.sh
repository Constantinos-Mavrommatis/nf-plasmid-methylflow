#!/bin/bash -ue
extract_collapse.R     --in-parquet EGF_met_1.parquet     --out-parquet EGF_met_1.extract_collapsed.parquet     --sample-id EGF_met_1     --plasmid-id PlasmidToy     --treatment A     --enzyme-conc 1x     --replicate 1     --run-id toy_run
