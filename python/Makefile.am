CLEANFILES = *~

EXTRA_DIST = \
	setup.py.in \
	eos/__init__.py \
	eos/analysis.py \
	eos/analysis_TEST.d \
	eos/analysis_TEST.py \
	eos/analysis_file.py \
	eos/analysis_file_TEST.d \
	eos/analysis_file_TEST.py \
	eos/analysis_file_description.py \
	eos/config.py \
	eos/constraint.py \
	eos/datasets.py \
	eos/datasets_file_description.py \
	eos/deserializable.py \
	eos/ipython.py \
	eos/observable.py \
	eos/parameter.py \
	eos/reference.py \
	eos/signal_pdf.py \
	eos/tasks.py \
	eos/tasks_TEST.py \
	eos/pyhf_likelihood.py \
	eos/data/__init__.py \
	eos/data/native.py \
    eos/data/native_TEST.d \
	eos/figure/__init__.py \
	eos/figure/common.py \
	eos/figure/config.py \
	eos/figure/data.py \
	eos/figure/figure.py \
	eos/figure/item.py \
	eos/figure/plot.py \
	eos/plot/__init__.py \
	eos/plot/config.py \
	eos/plot/plotter.py

AM_TESTS_ENVIRONMENT = \
			 export EOS_TESTS_CONSTRAINTS="$(top_srcdir)/eos/constraints"; \
			 export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters"; \
			 export PYTHONPATH="$(top_builddir)/python/.libs/:$(top_srcdir)/python"; \
			 export SOURCE_DIR="$(abs_srcdir)";

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS =

if EOS_ENABLE_PYTHON

eosdir = $(pkgpythondir)
eos_PYTHON =
eos_SCRIPTS = \
	eos/__init__.py \
	eos/analysis.py \
	eos/analysis_file.py \
	eos/analysis_file_description.py \
	eos/config.py \
	eos/constraint.py \
	eos/datasets.py \
	eos/datasets_file_description.py \
	eos/deserializable.py \
	eos/ipython.py \
	eos/observable.py \
	eos/parameter.py \
	eos/reference.py \
	eos/signal_pdf.py \
	eos/tasks.py \
	eos/pyhf_likelihood.py

eosdatadir = $(pkgpythondir)/data
eosdata_SCRIPTS = \
	eos/data/__init__.py \
	eos/data/native.py

eosplotdir = $(pkgpythondir)/plot
eosplot_SCRIPTS = \
	eos/plot/__init__.py \
	eos/plot/config.py \
	eos/plot/plotter.py

eosfiguredir = $(pkgpythondir)/figure
eosfigure_SCRIPTS = \
	eos/figure/__init__.py \
	eos/figure/common.py \
	eos/figure/config.py \
	eos/figure/data.py \
	eos/figure/figure.py \
	eos/figure/item.py \
	eos/figure/plot.py

_eosdir = $(pythondir)
_eos_LTLIBRARIES = _eos.la
_eos_la_SOURCES = \
	_eos.cc \
	_eos/external-log-likelihood-block.cc _eos/external-log-likelihood-block.hh \
	_eos/external-observable.cc _eos/external-observable.hh \
	_eos/log.cc _eos/log.hh \
	_eos/version.cc _eos/version.hh \
	_eos/wrappers.cc _eos/wrappers.hh
_eos_la_CXXFLAGS = @BOOST_PYTHON_CXXFLAGS@ @PYTHON_CXXFLAGS@ @AM_CXXFLAGS@ @GSL_CXXFLAGS@
_eos_la_LDFLAGS = -module -avoid-version -export-symbols-regex PyInit__eos @BOOST_PYTHON_LDFLAGS@
_eos_la_LIBADD = $(top_builddir)/eos/libeos.la $(top_builddir)/eos/utils/libeosutils.la -lboost_python$(BOOST_PYTHON_SUFFIX)

TESTS = \
	eos_TEST.py \
	eos/analysis_TEST.py \
	eos/analysis_file_TEST.py \
	eos/tasks_TEST.py \
	eos/data/native_TEST.py \
	eos/observable_TEST.py \
	eos/parameter_TEST.py \
	eos/figure/figure_TEST.py \
	eos/figure/item_TEST.py \
	eos/figure/plot_TEST.py \
	eos/plot/plotter_TEST.py

EXTRA_DIST += $(TESTS)

PYPI_SHARED_OBJECTS = 'libboost|libeos|libgfortran|libgsl|libquadmath|libsatlas|libyaml|libz'
eoshep-before: setup.py all
	echo preparation to build wheel in "$(DESTDIR)" for version $(PYPI_VERSION) using "$(PYTHON)" \(`$(PYTHON) -V`\) ; \
    cp $(abs_top_builddir)/python/setup.py $(DESTDIR)/ ; \
    cp $(abs_top_builddir)/python/.libs/_eos.so $(DESTDIR)/ ; \
    echo "$(PYPI_SHARED_OBJECTS)" ; \
    for f in `ldd $(DESTDIR)/_eos.so | grep '=>' | grep -E $(PYPI_SHARED_OBJECTS) | awk '{print $$3}'` ; \
    do \
        echo copying $$f ; \
        cp $$f $(DESTDIR)/ ; \
    done ; \
    chrpath -r '$$ORIGIN/_eos_libs' $(DESTDIR)/_eos.so ; \
    for f in $(DESTDIR)/libeos*.so* ; \
    do \
        chrpath -l $$f || continue ; \
        chrpath -r '$$ORIGIN' $$f ; \
    done ; \
    cp -r $(abs_top_srcdir)/python/eos $(DESTDIR)/ ; \
    mkdir $(DESTDIR)/_eos_data ; \
    for d in references.yaml constraints parameters ; \
    do \
        cp -r $(abs_top_srcdir)/eos/$$d $(DESTDIR)/_eos_data/ ; \
    done ; \
    echo 'is_wheel = True' >> $(DESTDIR)/eos/config.py ; \
    echo 'int foo() { return -1; }' > $(DESTDIR)/dummy.cc ; \
    mkdir $(DESTDIR)/bin ; \
    for f in $(abs_top_builddir)/src/scripts/* ; \
    do \
        cp -r $$f $(DESTDIR)/bin/ ; \
    done ; \
    cp $(abs_top_srcdir)/README.md $(DESTDIR)/ ; \
    pushd $(DESTDIR) ; \
    $(PYTHON) $(DESTDIR)/setup.py build ; \
    popd ; \
    for f in `find $(DESTDIR)/build/ -name "_eos.*.so"` ; \
    do \
        cp $(DESTDIR)/_eos.so $$f ; \
    done

endif
