export EOS_BASE_DIRECTORY=/tmp/btopilnu

## Find clusters within samples of posterior 'th+exp' defined in 'btopilnu.analysis'

eos-analysis                            \
    find-clusters                       \
    -f btopilnu.analysis                \
    $EOS_BASE_DIRECTORY/th+exp/clusters \
    $EOS_BASE_DIRECTORY/th+exp/mcmc-*   \
    -t 1.5                              \
    -c 2
