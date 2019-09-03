#!/usr/bin/env bash

python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results equivalent_width -b 0;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results equivalent_width -b 1;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_15;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq -s -10 -e 20;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq -s -5 -e 20;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results color_chisq -s 5 -e 50;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results spec_chisq -e .03 -t .1 -b csp_dr3_u csp_dr3_g csp_dr3_r csp_dr3_i csp_dr3_B csp_dr3_V csp_dr3_Y csp_dr3_J csp_dr3_H;
python run_analysis.py -m CMFGEN,1.04 CMFGEN,1.02 CMFGEN,1.4 CMFGEN,1.7 -o results spec_chisq -e .03 -f pW1 pW2 pW3 pW4 pW5 pW6 pW7 pW8;
