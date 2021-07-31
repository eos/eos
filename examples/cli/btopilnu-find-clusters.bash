export EOS_BASE_DIRECTORY=/tmp/btopilnu

## Find clusters within samples of posterior 'th+exp' defined in 'btopilnu.analysis'

eos-analysis                            \
    find-clusters                       \
    th+exp                              \
    -f btopilnu.analysis                \
    -t 1.5                              \
    -c 2
