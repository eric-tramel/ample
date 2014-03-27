ample
=====

A General implementation of Approximate Message Passing (AMP) for Matlab. 

***Note:*** This repository is still under active development. There may be some
hiccups.

Intro
-----
The goal of this problem is to provide a quick-to-use and easily generalizable
implementation of the approximate message passing (AMP) algorithm....
* ***Quick-to-use:*** At its core, one only needs to provide the set of obstervations,
the system which obtained those observations, and the type of signal prior. However,
`ample` contains a number of options to allow the user to tailor AMP to suit their
particular needs.
* ***Generalizable:*** `ample` provides a framework by which the user can easily 
implement their own priors to fit AMP to their specific application.

Acknowledgements
----------------
Many thanks to Jean Barbier and Andre Manoel for their moment and update calculations
for the supplied priors.
