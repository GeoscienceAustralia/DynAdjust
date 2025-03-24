[![cmake workflow](https://github.com/geoscienceaustralia/DynAdjust/actions/workflows/cmake_release.yml/badge.svg)](https://github.com/geoscienceaustralia/DynAdjust/actions/workflows/cmake_release.yml)
[![cmake tests](https://github.com/geoscienceaustralia/DynAdjust/actions/workflows/test_coverage.yml/badge.svg)](https://github.com/geoscienceaustralia/DynAdjust/actions/workflows/test_coverage.yml)
[![docker build](https://github.com/GeoscienceAustralia/DynAdjust/actions/workflows/docker.yml/badge.svg)](https://github.com/GeoscienceAustralia/DynAdjust/actions/workflows/docker.yml)
[![codacy badge](https://img.shields.io/codacy/grade/a3944cda0c72445f8a13b1f82b64f714)](https://app.codacy.com/gh/GeoscienceAustralia/DynAdjust/dashboard)
[![cii best practices](https://img.shields.io/badge/cii%20best%20practices-passing-success)](https://bestpractices.coreinfrastructure.org/projects/4894)
[![Coverage Status](https://coveralls.io/repos/github/GeoscienceAustralia/DynAdjust/badge.svg)](https://coveralls.io/github/GeoscienceAustralia/DynAdjust)
[![codecov status](https://img.shields.io/codecov/c/github/GeoscienceAustralia/dynadjust)](https://codecov.io/gh/GeoscienceAustralia/DynAdjust)
[![documentation (user's guide)](https://img.shields.io/badge/docs-usersguide-blue.svg)](https://github.com/geoscienceaustralia/DynAdjust/raw/master/resources/DynAdjust%20Users%20Guide.pdf)
[![gitHub license](https://img.shields.io/badge/license-Apache-blue.svg)](https://raw.githubusercontent.com/geoscienceaustralia/DynAdjust/master/LICENSE)
[![github releases](https://img.shields.io/github/v/release/geoscienceaustralia/DynAdjust)](https://github.com/geoscienceaustralia/DynAdjust/releases)
[![github issues](https://img.shields.io/github/issues/geoscienceaustralia/DynAdjust.svg)](https://github.com/geoscienceaustralia/DynAdjust/issues)
[![github downloads](https://img.shields.io/github/downloads/geoscienceaustralia/DynAdjust/total)](https://tooomm.github.io/github-release-stats/?username=geoscienceaustralia&repository=DynAdjust&target=_blank)
[![Docker Pulls](https://img.shields.io/docker/pulls/agrsga/dynadjust)](https://hub.docker.com/r/agrsga/dynadjust)
[![lines of code](https://tokei.rs/b1/github/geoscienceaustralia/dynadjust)](https://github.com/geoscienceaustralia/DynAdjust/tree/master/dynadjust)




[![DynAdjust](https://github.com/GeoscienceAustralia/DynAdjust/raw/master/resources/img/dynadjust-banner-sml.png)](https://github.com/GeoscienceAustralia/dynadjust/releases)

## Contents

- [Contents](#contents)
- [Overview](#overview)
- [Installing DynAdjust](#installing-dynadjust)
- [Building DynAdjust](#building-dynadjust)
- [Contributing to DynAdjust](#contributing-to-dynadjust)
- [Feedback](#feedback)
- [User's guide](#users-guide)
- [Creating images](#creating-images)
- [Licence details](#licence-details)

## Overview

DynAdjust is a rigorous, high performance least squares adjustment application. It has been designed
to estimate 3D station coordinates and uncertainties for both small and extremely large geodetic networks,
and can be used for the adjustment of a large array of Global Navigation Satellite System
(GNSS) and conventional terrestrial survey measurement types. On account of the phased adjustment
approach used by DynAdjust, the maximum network size which can be adjusted is effectively
unlimited, other than by the limitations imposed by a computer’s processor, physical memory and
operating system memory model. Example projects where DynAdjust can and has been used include
the adjustment of small survey control networks, engineering surveys, deformation monitoring
surveys, national and state geodetic networks and digital cadastral database upgrade initiatives.

DynAdjust provides the following capabilities:

- Import of data in geographic, cartesian and/or projection (UTM) coordinates contained in DNA, DynaML and SINEX data formats;
- Input of a diverse range of measurement types;
- Transformation of station coordinates and measurements between several static and dynamic reference frames;
- Rigorous application of geoid–ellipsoid separations and deflections of the vertical;
- Simultaneous (traditional) and phased adjustment modes;
- Automatic segmentation and adjustment of extremely large networks in an efficient manner;
- Rigorous estimation of positional uncertainty for all points in a network;
- Detailed statistical analysis of adjusted measurements and station corrections;
- Production of high quality network plots;
- Automated processing and analysis with minimal user interaction.

## Installing DynAdjust

To install DynAdjust on Linux, Mac or Windows, please visit the [installation page](resources/INSTALLING.md).

## Building DynAdjust

To build DynAdjust from source code on Linux, Mac or Windows, please visit the [building page](resources/BUILDING.md).

## Contributing to DynAdjust

Pleae refer to the [Contributing guidelines](./.github/CONTRIBUTING.md) for information on how to contribute to DynAdjust.

## Feedback

To suggest an enhancement to the functionality of DynAdjust, or to report a defect or unexpected behaviour, please submit your query via the [issue tracker](https://github.com/GeoscienceAustralia/dynadjust/issues).

## User's guide

A comprehensive User's Guide can be found in the [resources](https://github.com/GeoscienceAustralia/DynAdjust/tree/master/resources) folder.  The User's Guide provides information about the history of DynAdjust, its architecture and algorithms, its usage, supported file format specifications and basic command-line examples.

## Creating images

DynAdjust provides a capability to generate publication-quality images of raw station and measurement data and adjustment results, including shift vectors and estimated uncertainty, in a variety of projection types. This capability is made available by the program **plot**, the command line reference for which is documented in Appendix A.7 of the [User's guide](#users-guide).

[![DynAdjust plot examples](https://raw.githubusercontent.com/GeoscienceAustralia/DynAdjust/master/resources/img/dynadjust-plot-images.png)](https://github.com/GeoscienceAustralia/dynadjust/releases)
  
DynAdjust uses the **Generic Mapping Tools (GMT)**, available from the [GMT website](https://www.generic-mapping-tools.org/download/), and **gnuplot**, available from the [gnuplot homepage](http://www.gnuplot.info/). The current version of DynAdjust supports GMT version 6 and gnuplot version 5.4. To install GMT, please refer to the [install instructions](https://github.com/GenericMappingTools/gmt/blob/master/INSTALL.md), or the [build instructions](https://github.com/GenericMappingTools/gmt/blob/master/BUILDING.md) to build GMT from source. To install gnuplot, please refer to the [gnuplot download](http://www.gnuplot.info/download.html) page.

## Licence details

DynAdjust has an Apache 2.0 Licence - http ://www.apache.org/licenses/LICENSE-2.0.

NOTE: DynAdjust makes use of Boost C++, Apache's Xerces-C++ XML Parser (Apache 2.0 Licence) and CodeSynthesis XSD code. Hence, the following licence agreements will also need to be taken into account with the Apache 2.0 Licence.

- <https://www.boost.org/users/license.html>
- <https://www.codesynthesis.com/products/xsd/license.xhtml>

The free licence of CodeSynthesis XSD is GPL2, which requires any software that uses it to also be open source.  However, the CodeSynthesis site above states that you can use XSD generated code in proprietary applications provided that the lines of code do not exceed 10,000 lines.  The files generated for DynaML from XSD (dnaparser_pimpl.hxx/cxx and dnaparser_pskel.hxx/cxx), all of which have been heavily modified after they were originally generated, contain less than 10,000 lines.

DynAdjust also makes use of Intel's performance libraries - Math Kernel Library (MKL) and Threaded Building Blocks (TBB). Given the size of the associated libraries, it is not possible to upload the binaries to this site. Please visit the following website to download the latest versions:

- <https://software.seek.intel.com/performance-libraries>

As stated on this website, Intel's performance libraries are _free to use for personal and commercial applications_.
