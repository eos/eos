export EOS_BASE_DIRECTORY=/tmp/btopilnu

## Predict observables for each posterior sample.

eos-analysis                            \
    predict-observables                 \
    -f btopilnu.analysis                \
    differential                        \
    $EOS_BASE_DIRECTORY/th+exp/pmc
