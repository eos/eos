export EOS_BASE_DIRECTORY=/tmp/btopilnu

## Plot samples of posterior 'th+exp' defined in 'btopilnu.analysis'

##  - Use samples in the directories /tmp/btopilnu/{mcmc-*,pmc}
eos-analysis                            \
    plot-samples                        \
    -f btopilnu.analysis                \
    th+exp
