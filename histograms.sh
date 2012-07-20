#!/bin/bash

python histograms.py "initial" "hx"
python histograms.py "initial" "hy"
python histograms.py "initial" "hz"
python histograms.py "initial" "dm"
python histograms.py "initial" "hc"


python histograms.py "final" "hx"
python histograms.py "final" "hy"
python histograms.py "final" "hz"
python histograms.py "final" "dm"
python histograms.py "final" "hc"