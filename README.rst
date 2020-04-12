.. image:: https://img.shields.io/travis/sclamons/murraylab_tools/master.svg
   :alt: Testing with TravisCI
   :target: https://travis-ci.org/sclamons/murraylab_tools/branches
.. image:: https://img.shields.io/codecov/c/github/sclamons/murraylab_tools/master.svg
   :alt: TravisCI test coverage
   :target: https://codecov.io/github/sclamons/murraylab_tools/

Tools for Common Murray Lab Protocols
=====================================

Code for common tasks in the Murray lab. Supports Python >=3.6 only.

Install with the following (probably requires sudo): ``pip install git+git://github.com/sclamons/murraylab_tools.git@master``

To update, run (also with superuser privilege): ``pip install --upgrade --no-deps git+git://github.com/sclamons/murraylab_tools.git@master``

Currently includes the following subpackages:

echo: Echo Setup for TX-TL Experiments (and others)
=============================================

Code to produce Echo source plate setup instructions and picklists for one of the following cases:

* A TX-TL experiment with a TX-TL setup spreadsheet.
* A 2-dimensional dilution series, in TX-TL.
* Association of one or more substances on one or more source plates (not in TX-TL)

For details, see the "Echo Setup Usage Examples" ipython notebook under "examples".

biotek: Tidying and Analysis of Biotek Plate Reader Data
========================================================

Code for converting raw excel/CSV data from a Biotek plate reader into tidy format.

Also contains convenience functions for analysis of Biotek timecourse data:

* Background subtraction
* Endpoint averaging

For details, see the "Biotek Analysis Usage Examples" ipython notebook under "examples".
