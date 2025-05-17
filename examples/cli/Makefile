.PHONY: examples

# Note that some examples have to be executed in order, since they may use
# a previous examples output files stored in /tmp.
LINUX_EXAMPLES = \
    btopilnu-evaluate.bash \
    btopilnu-sample-mcmc.bash \
    btopilnu-find-mode.bash \
    btopilnu-find-clusters.bash \
    btopilnu-sample-pmc.bash \
    btopilnu-plot-samples.bash \
    btopilnu-predict-observables.bash \
    btopilnu-plot.bash \
    btopilnu-draw-figure.bash

MACOSX_EXAMPLES = \
    btopilnu-evaluate.bash \
    btopilnu-sample-mcmc.bash \
    btopilnu-find-mode.bash \
    btopilnu-find-clusters.bash \
    btopilnu-sample-pmc.bash \
    btopilnu-predict-observables.bash

ifeq ($(shell uname),Darwin)
    EXAMPLES=$(MACOSX_EXAMPLES)
else
    EXAMPLES=$(LINUX_EXAMPLES)
endif

examples: $(EXAMPLES)
	for x in $(EXAMPLES) ; do \
	    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" ; \
	    echo "# Running CLI example '$$x'" ; bash $$x || exit 1 ; \
	    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" ; \
	done
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" ; \
	echo "# All CLI examples ran successfully " ; \
	echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" ; \
