#!/bin/bash

python unitsphere.py "initial"
python unitsphere.py "final"

python field.py "initial"
python field.py "final"

./histograms.sh
./animate.sh