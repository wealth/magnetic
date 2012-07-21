#!/bin/bash

python unitsphere.py "initial"
python unitsphere.py "final"

python field.py "initial"
python field.py "final"

python black-field.py "initial"
python black-field.py "final"

./histograms.sh
./animate.sh