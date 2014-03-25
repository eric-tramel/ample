ample
=====

General Implementation of Approximate Message Passing for Matlab

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
