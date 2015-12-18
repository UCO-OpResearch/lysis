#! /bin/csh -f
set noglob
set verbose
libtbx.python "/Users/bpaynter/ownCloud/_School/_Research/Computation/BloodClotting/code/fable_sources/libtbx/run_tests.py" $*
libtbx.python "/Users/bpaynter/ownCloud/_School/_Research/Computation/BloodClotting/code/fable_sources/fable/run_tests.py" $*
