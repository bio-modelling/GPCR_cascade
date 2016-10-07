Drosophila Phototransduction Simulator
======================================

by Konstantin Nikolic and Joaquim Loizu 

Code developed in MATLAB using the Parallel Computing Toolbox, which allows computations 
on multicore computers and computer clusters. An appropriate graphic user interface was 
developed which gives very convenient and instructive presentation of the parameters 
used in the modelling.

Code covered by the BSD License (Matlab Central, File exchange), please see 
BSD_License.txt file or visit: 
http://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=22942

Please cite papers [1] and [2] if you use this code.

User Guide
----------

The main file is "Drosophila_Phototransduction.m" which when you run it opens the GUI.
The GUI contains the preloaded values for parameters and should be mostly 
self-explanatory in combination with the paper [1]. It is possible to run stochastic 
and deterministic models for single and multiple runs. The output consists of a table 
for some average values and three graphs (for the average quantum bump, latency times 
distribution and peak current distribution). It is possible to produce some other plots 
currently commented out in gui_averageQBall.m and gui_singleQB_multi.m files.

In the stochastic mode QBs are filtered with (100Hz) low pass filter.  
Averaging: QBs aligned by their half bump rise and fall times S. R. Henderson, H. Reuss
and R. C. Hardie, J. Physiol. 524, 2000, pp. 179-194, which was used as a source of 
experimental values.

time: time period in ms you want to simulate (e.g. 120)  
tstep: time step in ms for numerical solution (<=0.1)  
flagM/G/P/D/T: to set deterministic mode (0) or stochastic mode (1) for each step 
of the cascade  

Any additional questions (but please note that the support for the software can
NOT be guaranteed): k.nikolic@imperial.ac.uk


References
----------

1. Nikolic K, Loizu J, Degenaar P, Toumazou C: "A stochastic model of the 
   single photon response in Drosophila photoreceptors", Integrative Biology, 
   2010, Vol:2, Pages:354-370. [doi: 10.1039/C0IB00031K](https://doi.org/10.1039/C0IB00031K)
2. Nikolic K, Loizu J, Degenaar P: "Computational Modelling of the Drosophila 
   Phototransduction Cascade", Biophysical Journal, Vol:98, Issue 3, Supplement 1, 
   January 2010, Page 495a. [doi: 10.1016/j.bpj.2009.12.2695](http://dx.doi.org/10.1016/j.bpj.2009.12.2695)
