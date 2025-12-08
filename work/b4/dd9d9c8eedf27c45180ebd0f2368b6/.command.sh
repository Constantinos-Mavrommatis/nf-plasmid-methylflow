#!/bin/bash -ue
pileup_one.R     --in-parquet EGF_met_1.parquet     --out-parquet EGF_met_1.pileup.parquet     --mod-prob 0.8     --mod-code m     --sample-id EGF_met_1     --plasmid-id PlasmidToy     --treatment A     --enzyme-conc 1x     --replicate 1     --run-id toy_run
