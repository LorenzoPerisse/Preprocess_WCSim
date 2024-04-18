This file can be used to transform a WCSim output into a regular ROOT file.
It extracts all hit and interaction vertex data and put them into a TTree.
Additionaly, it can patch signal events over darknoise events, with specific index values for noise hits (index=0) and signal hits (index=1).
The patched events are then split into windows event that can be used to train machine learning software.
