.PHONY: examples

# Note that some examples have to be executed in order, since they may use
# a previous examples output files stored in /tmp.
LINUX_EXAMPLES = \
	basics.ipynb \
	predictions.ipynb \
	inference.ipynb \
	analysis-organisation.ipynb \
	simulation.ipynb

MACOSX_EXAMPLES =

ifeq ($(shell uname),Darwin)
    EXAMPLES=$(MACOSX_EXAMPLES)
else
    EXAMPLES=$(LINUX_EXAMPLES)
endif

%.nbconvert.ipynb: %.ipynb
	jupyter nbconvert --to notebook --execute $<

examples: $(patsubst %.ipynb,%.nbconvert.ipynb,$(EXAMPLES))
