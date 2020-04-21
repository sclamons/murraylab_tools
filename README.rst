.. image:: https://img.shields.io/travis/sclamons/murraylab_tools/master.svg
   :alt: Testing with TravisCI
   :target: https://travis-ci.org/sclamons/murraylab_tools/branches
.. image:: https://img.shields.io/codecov/c/github/sclamons/murraylab_tools/master.svg
   :alt: TravisCI test coverage
   :target: https://codecov.io/github/sclamons/murraylab_tools/

Tools for Common Murray Lab Protocols
=====================================

Code for common tasks in the Murray lab. Supports Python >=3.6 only.

Install with the following (probably requires sudo): ``pip install git+git://github.com/biocircuits/murraylab_tools.git@master``

To update, run (also with superuser privilege): ``pip install --upgrade --no-deps git+git://github.com/biocircuits/murraylab_tools.git@master``

Currently includes the following subpackages:

utilities: Simple Quality-Of-Life Improvement Tools
===================================================

* A code block that will automatically text or email you when your code is done running.

echo: Echo Setup for TX-TL Experiments (and others)
=============================================

Code to produce Echo source plate setup instructions and picklists for one of the following cases:

* A TX-TL experiment with a TX-TL setup spreadsheet.
* Programmatically-constructed TX-TL experiments, including easy 2D dilution series setup plus arbitrary additions.
* Fully general well contents, as defined by a spreadsheet of source materials and a spreadsheet listing what goes in each.

For details, see the "Echo Setup Usage Examples" ipython notebook under "examples".

biotek: Tidying and Analysis of Biotek Plate Reader Data
========================================================

Code for converting raw excel/CSV data from a Biotek plate reader into tidy format, plus a class to make plotting Biotek data traces quick and easy (BiotekCellPlotter).

Also contains convenience functions for analysis of Biotek timecourse data:

* Fluorescence calibration against lab data (TX-TL data only)
* Background subtraction
* Endpoint averaging
* Data smoothing (spline fit or moving average)
* Smoothed derivative calculation
* OD normalization
* Growth curve summarazation with logistic curve fits

For details, see the "Biotek Analysis Usage Examples" ipython notebook under "examples".
