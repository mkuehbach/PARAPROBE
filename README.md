# paraprobe-toolbox

This repository summarizes the initial work on what has now become a software package called the
paraprobe-toolbox. This is a collection of software for strong-scaling distributed- and shared-memory
CPU-parallelized tools for mining atom probe tomography data as well as GPU-parallelized tools
for atom probe crystallography. The software was developed at the Max-Planck-Institut für
Eisenforschung GmbH (MPIE) in Düsseldorf.

This repository contains the relevant code of the initial work on the project.
I consider the version here v0.1. This version contained a number of post-processing routines
that were integrated into a single C++ application. With more and more requests and questions to
address, maintaining this implementation became no longer efficient. Consequently, I separated the individual
functionalities into own tools, added a Python scripting interface, all with the aim to make the project
easier for me to maintain, to develop the tools functionally further, and thus improve their value
for the community.

Fast forward to 2022, the project has grown with now a number of case studies seen which
lead to a functional improvement. Some of them were reported in publications that have lead
to version v0.2 of the software. Version v0.2 is documented with the following
repository http://gitlab.mpcdf.mpg.de/mpie-aptfim-toolbox/paraprobe/

It was nice to grow with the project as a programmer; and while doing so, the time in fact great
and very instructive time at the MPIE came to an end. I am thankful for having had the chance of
working together with many very supportive colleagues, especially Franz Roters, Baptiste Gault,
Andrew Breen, and Dierk Raabe to name a few. With the German National Research Data Infrastructure
in the building kicking off at around 2020 I felt it was time to move on with my journey in the
field of data-centric research at the interface of (applied) condensed-matter physics and
materials engineering.

To assure that the benefits of the paraprobe-toolbox get not lost and - most importantly - 
to remain able in the future to offer the software via an MPCDF/MPIE-independent repository,
I have decided to communicate the new developments of the software via a different repository:

https://gitlab.com/paraprobe/paraprobe-toolbox

Therefore, the specific code in this repository will no longer be maintained but the majority
of the functionalities have been migrated and developed further. A preprint with a summary
of the newer developments to functionalities of the paraprobe-toolbox is available here:

https://arxiv.org/abs/2205.13510

Thanks to the continued support and many atom probers, I am thankful that the project can
continue. With the integration into the research data management system NOMAD OASIS
of the FAIRmat project, I am convinced that we are heading for a bright and interesting future
with the paraprobe-toolbox to continue supporting scientists who wish to apply complementary
computational methods on their atom probe data which are off-the-beaten path and which show
a glimpse of what is possible with open-source community-driven software.

Please consider the code in the new repository the relevant one in the future.
As usual feel free to contact me with specific questions about the paraprobe-toolbox
and questions about research data management related to e.g. atom probe, electron microscopy,
and multi-scale materials modeling techniques:
* https://www.fair-di.eu/fairmat/about-fairmat/team-fairmat