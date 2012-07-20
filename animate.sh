#!/bin/bash

python diameters.py "initial"
python diameters.py "final"

ffmpeg -sameq -i images/animate/initial-dms-animated-%03d.png output/initial-dms.mp4
ffmpeg -i output/initial-dms.mp4 -pix_fmt rgb24 -s 640x360 -loop 0 output/initial-dms.gif

ffmpeg -sameq -i images/animate/final-dms-animated-%03d.png output/final-dms.mp4
ffmpeg -i final-dms.mp4 -pix_fmt rgb24 -s 640x360 -loop 0 final-dms.gif