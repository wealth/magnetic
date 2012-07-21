#!/bin/bash

python diameters.py "initial"
python diameters.py "final"

ffmpeg -sameq -i images/animate/initial-dms-animated-%03d.png video/initial-dms.mp4

ffmpeg -sameq -i images/animate/final-dms-animated-%03d.png video/final-dms.mp4
