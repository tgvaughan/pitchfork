Pitchfork
=========

[![Build Status](https://github.com/tgvaughan/pitchfork/workflows/Unit%2Fintegration%20tests/badge.svg)](https://github.com/tgvaughan/pitchfork/actions?query=workflow%3A%22Unit%2Fintegration+tests%22)

Pitchfork is a [BEAST 2](https://www.beast2.org) package for performing
Bayesian phylogenetic inference under models supporting explicit polytomies.
It is a work in progress, so results should not be trusted!

Installation
------------

To install Pitchfork:

1. Download and install [BEAST 2](https://www.beast2.org).
2. Launch the BEAUti application distributed with BEAST.
3. From the File menu select "Manage Packages".
4. Click the "Package repositories" button at the bottom of the Package Manager dialog box.
5. Select "Add URL" and enter the following URL:
   `https://tgvaughan.github.io/pitchfork/package.xml`.
6. Click the "Done" button, then select "pitchfork" from the list of packages.
7. Click the "Install/Upgrade" button.  Once installation is complete, restart BEAUti.

You should now be able to set up a tree analysis under the beta coalescent.
To do this, proceed as for any other BEAST analysis, but select the Beta coalescent
tree prior from the priors panel.  The resulting analysis XML will perform
Bayesian phylogenetic inference under this model, jointly inferring the effective
population size Nₑ and the α parameter of the beta distribution.

License
-------

Pitchfork is distributed under the terms of version 3 of the GNU
General Public License.  A copy of the license can be found in this
directory in the file COPYING.
