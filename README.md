### Modification Notes on NWS=7
NWS=7 is designed to read in wind stress, instead of wind speed, on a structured grid into the ADCIRC model in the same way as 
NWS=6. However the original code is still implemented as if NWS=2 (on a triangular mesh).


Correction has been applied to source codes in **prep**, **src** folders to make NWS=7 to perform as it was described on the ADCIRC manual. 
However, somehow the NWS=2 can not work properly in this version I uploaded and I couldn't locate the error source. I checked that I didn't make modification to any code related to NWS=2.


#### test cases;
two test cases have been provided under **testcase** for confirmation that the NWS=7 (cold start) works properly after correction.
