Tools for Common Murray Lab Protocols
=====================================

Code for common tasks in the Murray lab. Currently supported in Python 2.7.

Install with the following (probably requires sudo): pip install git+git://github.com/sclamons/murray_echo.git@master

Currently includes the following subpackages:

echo: Echo Setup for TX-TL Experiments (and others)
=============================================

Code to produce Echo source plate setup instructions and picklists for one of the following cases:

* A TX-TL experiment with a TX-TL setup spreadsheet.
* A 2-dimensional dilution series, in TX-TL.
* Association of one or more substances on one or more source plates (not in TX-TL)

See the example usage ipython notebook for usage details.

biotek: Tidying and Analysis of Biotek Plate Reader Data
========================================================

Code for converting raw excel/CSV data from a Biotek plate reader into tidy format.