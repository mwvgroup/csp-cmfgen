#!/usr/bin/env bash

python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results equivalent_width -b 0;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results equivalent_width -b 1;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_15;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq -s -10 -e 20;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq -s -5 -e 20;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq -s 5 -e 50;
