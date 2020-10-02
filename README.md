# progenitor_grid
=================

Tool for processing a calculated grid of MESA sdB progenitors.

## Installation
Install by cloning the repository, `cd` into it and then execute

    pip install .
    
to  install the package on your system.

## Uninstallation
Uninstall by executing

    pip uninstall progenitor_grid

## Sample usage

A directory with LOG directories containing MESA history.data and .mod files for progenitors with names formatted as f"model_{model_number:05d}.mod":

    `logs = "logs_test"`

Output .txt file for the final grid:

    `output_file = "test_grid.txt"`

Output directory for the selected .mod files:

    `output_dir = "test_dir"`

Initialize a grid:

    `g = progenitor_grid.ProgenitorGrid(logs, output_file, output_dir)`

Evaluate the grid. This method evaluated the grid, saves the output to a file and copies the selected .mod files.

    `g.evaluate_initial_grid()`