export EOS_BASE_DIRECTORY=/tmp/btopilnu

## Sample from posterior 'th+exp' defined in 'btopilnu.analysis'

##  - Output samples to the directory /tmp/btopilnu
##  - Sample using Markov chains, set index '0' to set up the seed
##    for the random number generation.
##  - Carry out 5 prerun steps with 500 samples per step (-p 10 -n 500).
##  - Produce further 5000 samples in the main run (-N 5000)
eos-analysis                            \
    sample-pmc                          \
    -f btopilnu.analysis                \
    th+exp                              \
    $EOS_BASE_DIRECTORY/th+exp/clusters \
    -s 5                                \
    -n 1000                             \
    -N 5000
