# Contributing guidelines

## Contents

- [Contributing guidelines](#contributing-guidelines)
  - [Contents](#contents)
  - [Preamble](#preamble)
  - [How to contribute](#how-to-contribute)
    - [Getting started](#getting-started)
    - [Raising and addressing an issue](#raising-and-addressing-an-issue)
    - [General discussion](#general-discussion)
    - [Private reports](#private-reports)
  - [Core functionality, algorithms and formulae](#core-functionality-algorithms-and-formulae)
  - [Coding conventions and guidelines](#coding-conventions-and-guidelines)
    - [Programming language](#programming-language)
    - [Coding style](#coding-style)
    - [Cross platform compatibility](#cross-platform-compatibility)
    - [Folder and file and structure](#folder-and-file-and-structure)
    - [Best practice](#best-practice)
    - [Release schedule](#release-schedule)
  - [Automated test suite](#automated-test-suite)
    - [Continuous integration](#continuous-integration)
    - [Test code coverage](#test-code-coverage)
    - [Static code analysis](#static-code-analysis)

## Preamble

Thank you for considering making a contribution to DynAdjust.

DynAdjust forms an integral component of public geospatial infrastructure and commercial software. It is used for the establishment and routine maintenance of Australia's Geospatial Reference System, the maintenance of International Terrestrial Reference Frame local ties, the computation of digital cadastre (or land parcel) fabrics, and a diverse array of engineering and geospatial projects. For these reasons, great emphasis is placed upon ensuring the DynAdjust code base is maintained using high quality coding standards.

The purpose of this page is to provide some contributing guidelines to help maintain a high quality code base, and to foster a positive and prosperous collaborative development culture.

## How to contribute

### Getting started

To make any changes to the DynAdjust code base, suggest new features, or improve the documentation using GitHub, you will need a [GitHub account](https://github.com/signup/free).

Please familiarise yourself with [GitHub's pull request process](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).  Another useful resource is Aaron Meurer's [tutorial](https://www.asmeurer.com/git-workflow/) on the git workflow.

If you plan to make changes to the code base, please ensure you have read the [installation instructions](./INSTALLING.md) and have obtained all the essential prerequisites.

### Raising and addressing an issue

[![GitHub Issues](https://img.shields.io/github/issues/geoscienceaustralia/DynAdjust.svg)](https://github.com/geoscienceaustralia/DynAdjust/issues)

The general process for addressing issues in DynAdjust is as follows.

1. Search the list of open issues in the [issue tracker](https://github.com/geoscienceaustralia/DynAdjust/issues) for any occurrence(s) of the issue you would like to address.  
   1. If it is closely related to an existing open issue, please add a comment to the open issue rather than to create a new issue.
   2. By default, the issue tracker displays open issues. Perhaps your issue relates to a defect previously marked as fixed, but remains unresolved. In this case, remove the `is:open` switch to search through [all issues](https://github.com/geoscienceaustralia/DynAdjust/issues?q=is%3Aissue).  If you feel the issue persists, reopen the issue and provide comments (as described in the next step).
   3. If the issue has not been raised before, continue to the next step.
2. Create a new issue, or reopen a closed issue.
   1. **Defects**. If it is a defect, please describe how to reproduce it and, if possible, provide sample command line arguments and snippets of the output which demonstrate the defect. Please select an appropriate [label](https://github.com/geoscienceaustralia/DynAdjust/issues/labels) that characterises the issue.
   2. **Enhancements**. If it is an enhancement, please provide the desired or expected behaviour not presently delivered by DynAdjust. Please cite any relevant technical documents (e.g. journal articles, reference texts or other publications) that may help the developers address the issue. Wherever possible, provide some test cases demonstrating expected behaviour.
   3. **Compilation**. If your issue relates to compiler errors or warnings, please provide full details of the operating system and compiler version, and provide a snippet of the compiler error or warning produced.
   4. **Feedback**. Perhaps you'd like to share your thoughts on how DynAdjust can be improved. In this context, you can:
      - Contribute code you already have. If it is not ready for production, we'd be glad to hear from you and to help you get it ready for the next release.
      - Propose a new function or suggest an alternative formula/algorithm.
      - Contribute suggestions, corrections or updates to the [documentation](https://github.com/geoscienceaustralia/DynAdjust/blob/master/resources/DynAdjust%20Users%20Guide.pdf).
      - Share some sample data.
3. Create your own separate [fork](https://github.com/geoscienceaustralia/DynAdjust/network/members) of DynAdjust, or rebase your forked copy of DynAdjust to the latest version.
4. Prepare for a specific pull request and begin making your changes.
   1. Ensure your changes relate to the discrete issue you're attempting to address.  That is, don't try to solve or fix everything in the one change request.
   2. For new DynAdjust features, create a new test script that tests your work and provides satisfactory [code coverage](https://coveralls.io/github/geoscienceaustralia/DynAdjust).
   3. Make sure all tests pass using the supplied [`make_dynadjust_gcc.sh`](https://github.com/geoscienceaustralia/DynAdjust/blob/master/resources/make_dynadjust_gcc.sh) script. For instance, to build and test your changes, run the script as: `$ make_dynadjust_gcc.sh -c -a -n -t` This will build the code base without cloning a fresh copy, no user interaction, no installation and will execute all cmake tests. For more information, please refer to the [installation instructions](./INSTALLING.md).
   4. Issue commits in logical packets of work, providing a clear message that describes what you're doing.
   5. When all commits for your change request are complete, issue a pull request. Please provide concise, clear and simple descriptions of your changes (and don't forget to check your spelling!). If required, please assign a reviewer you feel would be best placed to review your changes.
5. People with sufficient permissions will review your changes and approve the request when they are satisfied.  They may respond to you with questions or requests for further improvements.

### General discussion

Most discussions happen in the [issue tracker](https://github.com/geoscienceaustralia/DynAdjust/issues) or within pull requests. For general questions, please feel free to post to the repository's [discussion page](https://github.com/geoscienceaustralia/DynAdjust/discussions).

We encourage the use of language that is objective, professional, respectful and considerate of the wide user/developer audience. DynAdjust's audience has a diverse background, including theoretical geodesists, GIS users, software developers and senior level experts.  Not everyone will understand the level of detail being discussed, however, it should be simple enough for those with the relevant domain knowledge to understand and respond to. In all cases, we encourgage you to model positive communication behaviours.

### Private reports

Primarily, all issues, defects, enhancements and queries should be tracked publicly using the repository's [issue tracker](https://github.com/geoscienceaustralia/DynAdjust/issues). However, if you would like to submit a request or report privately (i.e. for cases when information should not be released publically), please direct your request ot report to [geodesy@ga.gov.au](mailto:geodesy@ga.gov.au).

## Core functionality, algorithms and formulae

[![Documentation (User's Guide)](https://img.shields.io/badge/docs-usersguide-red.svg)](https://github.com/geoscienceaustralia/DynAdjust/raw/master/resources/DynAdjust%20Users%20Guide.pdf)

DynAdjust implements a wide range of specialist geodetic and surveying algorithms and formulae. These have been sourced from journal articles, reference texts, published standards and other peer-reviewed publications. All functionality, algorithms and formulae have been documented in the [DynAdjust User's Guide](https://github.com/geoscienceaustralia/DynAdjust/blob/master/resources/DynAdjust%20Users%20Guide.pdf) and are appropriately cited with the full biblographic reference.  When proposing a new feature that implements an undocumented algorithm or formula, please cite its full reference with your comments.

Periodically, a DynAdjust Steering Committee meets to review the functionality of DynAdjust and to discuss potential enhancements.  If you would like the Steering Committee to consider your thoughts on the future of DynAdjust, please submit your comments to [geodesy@ga.gov.au](mailto:geodesy@ga.gov.au).

## Coding conventions and guidelines

### Programming language

DynAdjust development has evolved over 14 years from C++98 to C++14, with heavy reliance upon the Standard Template Library (STL). Wherever possible, [Boost](https://www.boost.org/) has been used to incorporate portable libraries unavailable within the C++ Standard.  Over time, some of Boost's libraries features have become a part of the C++ Standard.

Current development is focussed on C++14.  Various (older) parts of the code base is in C. Over time, these portions are being converted to C++.  While there are many desirable features of C++20, such features aren't being pursued at the present time so as to provide greater compatibility with older compilers.

### Coding style

DynAdjust has been in development since late 2008.  Consequently, the style of the code has evolved with the capability and experience of the developers. For this reason, different projects will exhibit different styles; so too will different functions within the same project. Where possible, please endeavour to conform with the style of the surrounding code.

### Cross platform compatibility

Please consider writing code that can be implemented on the widest number of platforms and systems (e.g. UNIX, Linux, Windows, Apple), rather than contributing code or using a library of functions that is specific to one particular system. [Boost](https://www.boost.org/) has a range of functions that can be implemented on various systems. In certain cases, it will be inevitable to provide code to work around the system-specific commands and behaviours using `#if defined` (or other means).

### Folder and file and structure

The folder and file structure of the code base and supporting files is as follows:

- Code that is specific to a particular DynAdjust program is located in the [`dynadjust/dynadjust/`](https://github.com/geoscienceaustralia/DynAdjust/tree/master/dynadjust/dynadjust) folder. The projects named `dna...wrapper` are programs with user interfaces that can be executed from the command line. The projects without the `...wrapper` suffix are libraries that provide the specific functionality relating to that project.
- Code that is general and consumed by all projects is located in the [`dynadjust/include/`](https://github.com/geoscienceaustralia/DynAdjust/tree/master/dynadjust/include) folder.
- User documentation and build scripts are located in the [`resources/`](https://github.com/geoscienceaustralia/DynAdjust/tree/master/resources) folder.
- Sample data that you can use to test DynAdjust functionality is located in the [`sampleData/`](https://github.com/geoscienceaustralia/DynAdjust/tree/master/sampleData) folder.

### Best practice

The use of consistent, high quality coding standards is key to delivering high quality software and minimising defects. We encourage you to adopt and adhere to coding standards that are well respected or acknowledged by industry experts.

For the most part, coding standards and principles for best practice adopted for DynAdjust have been inspired from many sources, some of which include (in alphabetical order):

- Josuttis, N. M. (1999). The C++ Standard Library - A Tutorial and Reference. 1st edition. Addison-Wesley.
- Lakos, J. (1996). Large-Scale C++ Software Design. Pearson Education Limited.
- Meyers, S. (2001). Effective STL: 50 specific ways to improve your use of the standard template library. Addison-Wesley Longman Ltd.
- Meyers, S. (2005). Effective C++: 55 Specific Ways to Improve Your Programs and Designs. 3rd edition. Addison-Wesley Professional.
- Meyers, S. (2014). Effective Modern C++: 42 Specific Ways to Improve Your Use of C++11 and C++14. O'Reilly Media, Inc.
- Schäling, B. (2014). The Boost C ++ libraries. XML Press. Available online: [https://boost.org](https://boost.org).
- Stroustrup, B. (2013). C++ Programming Language. 4th edition. Addison Wesley.
- Williams, A. (2012). C++ Concurrency in Action: Practical Multithreading. 1st edition. Manning Publications.

### Release schedule

[![GitHub Releases](https://img.shields.io/github/v/release/geoscienceaustralia/DynAdjust.svg)](https://github.com/geoscienceaustralia/DynAdjust/releases)

Although code changes can occur frequently, we endeavour to deliver new [releases](https://github.com/geoscienceaustralia/DynAdjust/releases) on a regular basis to help developers of other software packages and products that depend upon DynAdjust to plan accordingly.  The currently adopted release schedule cadence is six months.

## Automated test suite

The file [`make_dynadjust_gcc.sh`](https://github.com/geoscienceaustralia/DynAdjust/blob/master/resources/make_dynadjust_gcc.sh) provides a way to build the code and run tests on local copies and on GitHub via continuous integration (CI) script.

With every commit and pull request, the project is tested with a unit test suite managed by cmake (via [`CMakeLists.txt`](https://github.com/geoscienceaustralia/DynAdjust/blob/master/dynadjust/CMakeLists.txt)).

### Continuous integration

[![cmake tests](https://img.shields.io/github/workflow/status/geoscienceaustralia/dynadjust/Build,%20test%20and%20code%20coverage?label=cmake%20tests)](https://github.com/geoscienceaustralia/DynAdjust/actions/workflows/test.yml)

The DynAdjust repository uses an automated test suite managed by GitHub Actions and [Travis](https://travis-ci.org/github/geoscienceaustralia/DynAdjust) to ensure new or changed code is built and tested in an automated way. In this context, continuous integration (CI) has been set up to invoke builds for Linux using the script [`test_coverage.yml`](https://github.com/geoscienceaustralia/DynAdjust/actions/workflows/test_coverage.yml).

### Test code coverage

[![coveralls status](https://img.shields.io/coveralls/github/GeoscienceAustralia/DynAdjust)](https://coveralls.io/github/GeoscienceAustralia/DynAdjust)

To ensure the test suite adequately covers the code base (as much as practically possible), code coverage analysis is run on every commit and pull request.  The code coverage analysis tool used for DynAdjust is [coveralls.io](https://coveralls.io/github/geoscienceaustralia/DynAdjust). The developers are currently working on increasing the coverage provided by the cmake test suite (via [`CMakeLists.txt`](https://github.com/geoscienceaustralia/DynAdjust/blob/master/dynadjust/CMakeLists.txt)).

### Static code analysis

[![codecov status](https://img.shields.io/codecov/c/github/icsm-au/dynadjust)](https://codecov.io/gh/icsm-au/DynAdjust)
[![codacy badge](https://img.shields.io/codacy/grade/a3944cda0c72445f8a13b1f82b64f714)](https://app.codacy.com/gh/icsm-au/DynAdjust/dashboard)

With every commit and pull request, [Codecov](https://app.codecov.io/gh/icsm-au/DynAdjust) and [Codacy](https://app.codacy.com/gh/icsm-au/DynAdjust/dashboard) are run to perform static code analysis. The code base is undergoing regular review to address the remaining issues (some of which are insignificant).
