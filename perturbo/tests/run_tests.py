#!/usr/bin/env python3

from perturbopy import tests_use
import sys

result = tests_use.do_tests(sys.argv[1:])
sys.exit(result)
