# Ensure this script fails if any command fails
set -e

eos-analysis \
    draw-figure \
    -f btopilnu.analysis \
    -b /tmp/btopilnu \
    CKM-Vub

eos-analysis \
    draw-figure \
    -f btopilnu.analysis \
    -b /tmp/btopilnu \
    CKM-Vub-v-FF
