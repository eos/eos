export EOS_BASE_DIRECTORY=/tmp/btopilnu

## Find clusters within samples of posterior 'th+exp' defined in 'btopilnu.analysis'

eos-analysis                             \
    find-mode                            \
    -f btopilnu.analysis                 \
    th+exp                               \
    -c 0000
