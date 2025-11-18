#!/bin/bash
### Converts mp3 format videos into webm
### Must include formats (mp3 and webm) in arguments

# Arguments
input=$1
output=$2
# FFMPEG command
ffmpeg -i $input -c:v libvpx -crf 10 -b:v 2M -c:a libvorbis $output
exit 0
