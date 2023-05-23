# Readme

This project is for the simulation of SUSY quantum mechanics systems using quantum computing.  The most interesting things for users are going to be in the tutorials directory, which features jupyter notebooks with combinations of latex and code.

## Testing

At the moment the tests feature some example based testing.  To run all the tests and see the code coverage you will need the python packages `pytest` and `coverage`.  Then you can run the following two commands

    coverage run -m --source=src pytest
    coverage report -m
