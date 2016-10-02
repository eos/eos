.PHONY: examples

# Note that some examples have to be executed in order, since they may use
# a previous examples output files stored in /tmp.
EXAMPLES = \
    evaluate-btopilnu-integrated.bash \
    sample-mcmc-btopi-ff.bash \
    sample-pmc-btopi-ff.bash \
    find-mode-btopi-ff.bash \
    propagate-uncertainty-btotaunu.bash \
    propagate-uncertainty-btopilnu.bash

examples: $(EXAMPLES)
	for x in $(EXAMPLES) ; do \
	    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" ; \
	    echo "# Running example '$$x'" ; bash $$x || exit 1 ; \
	    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" ; \
	done
