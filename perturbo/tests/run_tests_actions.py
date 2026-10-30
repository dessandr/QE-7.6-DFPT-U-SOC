from perturbopy import tests_use
import sys
import os

test_number = os.environ.get("TEST_NUMBER")
config_machine = os.environ.get("CONFIG_MACHINE")
regime = os.environ.get("REGIME")
if regime == 'perturbo':
    result = tests_use.do_tests(["-s", '--config_machine', config_machine, '--epr', f'epr{test_number}'])
elif regime == 'qe2pert':
    result = tests_use.do_tests(["-s", '--config_machine', config_machine, '--run_qe2pert', '--epr', f'epr{test_number}'])
sys.exit(result)