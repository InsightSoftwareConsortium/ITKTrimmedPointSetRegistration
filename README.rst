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

This module implements a decorator to the PointSetMetric that overrides the accumulation of the value and derivative computation to use a trimmed number of points.

![jensen.pdf](Documentation/Figures/jensen.pdf)
![euclidean.pdf](Documentaiton/Figures/euclidean.pdf)
![trimmed-euclidean.pdf](Documentation/Figures/trimmed-euclidean.pdf)
