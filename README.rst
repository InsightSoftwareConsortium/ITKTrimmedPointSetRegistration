ITKTrimmedPointSetRegistration
=================================

.. image:: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_apis/build/status/itktrimmedpointsetregistration?branchName=master
    :target: https://dev.azure.com/InsightSoftwareConsortium/ITKModules/_build/latest?definitionId=8&branchName=master
    :alt:    Build Status

.. image:: https://img.shields.io/pypi/v/itk-trimmedpointsetregistration.svg
    :target: https://pypi.python.org/pypi/itk-trimmedpointsetregistration
    :alt: PyPI Version

.. image:: https://img.shields.io/badge/License-Apache%202.0-blue.svg
    :target: https://github.com/InsightSoftwareConsortium/ITKTrimmedPointSetRegistration/blob/master/LICENSE)
    :alt: License

Overview
--------

This is a module tfor trimmed point set registration

This module implmenets a trimmed euclidean point set metric to only use a subset of the data points for registration.


A decorator to the PointSetMetric that overrides the accumulation of the value and derivative computation to use a trimmed 
number of points is in the module as well but is as of yet not functional.


A Simple Example
----------------

Comaring Jensenm Euclidean and trimmed Euclinean metric on 2d toy example

.. image:: Documentation/Figures/jensen.png
.. image:: Documentation/Figures/euclidean.png
.. image:: Documentation/Figures/trimmed-euclidean.png
