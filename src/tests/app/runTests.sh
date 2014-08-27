#!/bin/bash

set -e
rm -rf out
wclean
wmake
tests --gtest_throw_on_failure
lcov -c --directory . --output-file coverage.info
genhtml coverage.info --output-directory out --demangle-cpp --no-function-coverage
rm coverage.info
kde-open out/index.html
