#!/usr/bin/env bash

python run_analysis.py -o results equivalent_width -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7;
python run_analysis.py -o results equivalent_width -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -f;
python run_analysis.py -o results lc_color         -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7;
