This is the visual object tracker KCFDP presented in:
[1] "Enable Scale and Aspect Ratio Adaptability in Visual Tracking with Detection Proposals", BMVC, 2015,
    Dafei Huang, Lei Luo, Mei Wen, Zhaoyun Chen and Chunyuan Zhang.

The implementation is built upon:
Color-name feature integration and model updating scheme:
[2] "Adaptive Color Attributes for Real-Time Visual Tracking", CVPR, 2014,
    Martin Danelljan, Fahad Shahbaz Khan, Michael Felsberg and Joost van de Weijer.
[3] "Learning color names for real-world applications", TIP, 18(7):1512-1524, 2009,
    J. van de Weijer, C. Schmid, J. J. Verbeek and D. Larlus.

Original CSK and KCF tracking framework:
[4] "Exploiting the circulant structure of tracking-by-detection with kernels", ECCV, 2012,
    J. F. Henriques, R. Caseiro, P. Martins and J. Batista.
[5] "High-Speed Tracking with Kernelized Correlation Filters", TPAMI, 2014,
    J. F. Henriques, R. Caseiro, P. Martins and J. Batista.
    http://www.isr.uc.pt/~henriques/circulant/
	
Structured Forests edge detector and Edge Boxes detection proposal generator:
[6] "Structured Forests for Fast Edge Detection", ICCV, 2013,
    P. Dollar and C. Zitnick.
[7] "Edge Boxes: Locating Object Proposals from Edges", ECCV, 2014,
    C. Zitnick and P. Dollar.
	
The IoU calculation code and example sequence along with annotations:
[8] "Online Object Tracking: A Benchmark", CVPR, 2013,
    Y. Wu, J. Lim and M.-H. Yang.
    http://visual-tracking.net/
	
Additional tools needed when running the code:
[9] "Piotr's Image and Video Matlab Toolbox (PMT)",
    P. Dollar.
    http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html

Codes above are integrated and modified by Dafei Huang.

Quick Start Guide:
Running the code directly on Girl[8] sequence:
1. Download and compile PMT[9];
2. Modify the path in Line 59 of run_tracker.m to your PMT path;
3. Go to the root directory of this code and run "run_tracker" in Matlab.

Integrating the code into OTB[8] tracking benchmark suite:
1. Download and prepare your environment according to [8];
2. Download and compile PMT[9];
3. Modify the path in Line 58 of run_KCFDP.m to your PMT path;
4. Copy the whole directory of this code into OTB_ROOT_PATH/tackers/;
5. Add a new line "struct('name','KCFDP','namePaper','KCFDP'),..." into
   the "trackers1" array in OTB_ROOT_PATH/util/configTrackers.m;
6. Run the OTB benchmark suite according to [8].

NOTE: 
1. For your convenience we have generated the binary files of [6] and [7]
   for 64-bit MAC OS, Windows, and Linux. 
   Please recompile the codes in ./private/ if needed.
2. The following files are part of Structured Forests[6] and Edge Boxes[7], and are provided for
   convenience only:
   edgeBoxes.m, edgesChns.m, edgesDetect.m, modelBsds.mat,
   edgeBoxesMex.cpp, edgesDetectMex.cpp, edgesNmsMex.cpp, spDetectMex.cpp and their relevant 
   binary files.
   These files are under the license specified in license_Structured_Forests_and_Edge_Boxes.txt.
3. The following files are part of ACT[2], and are provided for convenience only:
   im2c.m, get_feature_map.m, w2crs.mat.
   Please refer to readme_ACT.txt for the authorship information.
4. The tracking framework utilized here is from KCF[5] under the license specified in 
   license_KCF.txt.
5. calcRectInt.m and the example sequence along with annotations are from OTB[8] under the 
   GNU-GPL license.
6. The rest parts of KCFDP are distributed under the BSD license:

	  
Contact:
Dafei Huang
huangdafei1012@163.com






