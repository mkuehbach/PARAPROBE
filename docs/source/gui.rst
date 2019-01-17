**Graphical Front End**
=======================
To assist APT practitioners as much as possible in setting up PARAPROBE runs, we have developed a Python/Bokeh-based graphical user interface (GUI). 
It is accessible as a Bokeh wheel through the **scripts** folder and enables to generate an XML control file through a browser GUI.
Its main purpose is to assist the user in generating functional, i.e. format conformant control settings file.

The user interacts with the GUI to enter the analyses of interest. The GUI parses the user input and writes a conformant XML file.
Subsequently, this file can be used to run a PARAPROBE job.

We are currently working on an integration of this GUI into a virtual machine workflow that dispatches PARAPROBE jobs to a workstation queing system. Eventually, this will constitute the world's first fully automatized open source solution for instructing strongly scaling APT data mining operations without having the users to worry about the HPC details buzzing in the background.

.. figure:: ../images/PARAPROBEGUI_01.png
..   :scale: 20%
..   :align: left
..   :target: https://github.com/mkuehbach/PARAPROBE
