## Tests for Perturbo

This folder contains the files necessary to run the tests of the Perturbo package. You can read in detail how to run the tests in the <a href="https://perturbopy.readthedocs.io/en/latest/testsuite.html" target="_blank"> Perturbopy documentation</a>.

### Required steps:

1. Compile Perturbo
2. Install Perturbopy (details on installation can be found <a href="https://perturbopy.readthedocs.io/en/latest/usage.html" taget="_blank">here</a>)

3. Navigate to the `perturbo/tests` folder
    1. To run only the tests for perturbo.x, run the following command:
    ```bash
    ./run_tests.py
    ```
    2. To run tests for both perturbo.x and qe2pert.x, run:
    ```bash
    ./run_tests.py --run_qe2pert
    ```
