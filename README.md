mirna_stats
===========

Compute sensitivity and specificity of miRNA target prediction tools

Note on testing: Tests are provided using py.test. You need to set PYTHONPATH so that *import compute_stats* works for the tests, e.g.

    $ export PYTHONPATH=.
    $ py.test

Currently only MiRanda is supported, will support RNAhybrid and MicroTar in future
