#HMD Distortion Correction Toolbox for MATLAB

This toolbox provides a non-parametric distortion correction method for **Optical See-Through** Head-Mounted Display (OST-HMDs).

## How to Use it
On MATLAB, You can simply run our evaluation code as :
```Matlab
>> main_undistorty_optics
```
(We tested the code on MATLAB 2015b with Statistics Toolbox.)
If the program runs successfully, you would see an estimation result of an optical distortion learned from a sample dataset:

<img src="https://cloud.githubusercontent.com/assets/7195124/11901820/fa866912-a5ad-11e5-8c24-f0db24561fb8.jpg" width="350"/>


## Problem Background

When we wear an OST-HMD, optical aberration from the display's optics distorts our field of view: 

<img src="https://cloud.githubusercontent.com/assets/7195124/11900674/acf0eb24-a5a7-11e5-9452-b062d1366a3c.jpg" width="800"/>
 
This results in a misaligned Augmented Reality rendering, and our method aims to correct it:

<img src="https://cloud.githubusercontent.com/assets/7195124/11900867/96c26cb4-a5a8-11e5-9269-322140308ad0.jpg" width="500"/>

Our method is based on geometrical optics. We measure light rays in the scene trace how an optics of the display distorts the ray when we place an OST-HMD in front of the view point. By correcting many of these correspondences, we learn a mapping between the original light rays to the distorted. In out code, we achieve this by using a Kernel regression.

<img src="https://cloud.githubusercontent.com/assets/7195124/11900960/2f3cd344-a5a9-11e5-9fba-493ee54ec5b0.jpg" width="600"/>

For more detail, please refer to our papers listed in the bottom of this document .



## References:
The original paper is: <br>
[**Pre-print**](http://campar.in.tum.de/pub/itoh2015vr/itoh2015vr.pdf), 
[**Talk-slides**](http://campar.in.tum.de/pub/itoh2015vr/itoh2015vr.slides.pdf) <br>
```latex
@article{DBLP:journals/tvcg/ItohK15,
   author = {Yuta Itoh and Gudrun Klinker},
   title = {Light-Field Correction for Spatial Calibration of Optical See-Through Head-Mounted Displays},
   year = {2015},
   month = {April},
   journal = {IEEE Transactions on Visualization and Computer Graphics (Proceedings Virtual Reality 2015)},
   volume = {21},
   number = {4},
   pages = {471--480},
   doi = {http://doi.ieeecomputersociety.org/10.1109/TVCG.2015.2391859},
}
```

Follow-up paper which extends this calibration method to the image distortions (not included in the current code): <br>
[**Pre-print**](http://campar.in.tum.de/pub/itoh2015ismar2/itoh2015ismar2.pdf), 
[**Talk-slides**](http://campar.in.tum.de/pub/itoh2015ismar2/itoh2015ismar2.slides.pdf) <br>
```latex
@InProceedings{DBLP:conf/ismar/ItohK15,
  Title                    = {{Simultaneous Direct and Augmented View Distortion Calibration of Optical See-Through Head-Mounted Displays}},
  Author                   = {Yuta Itoh and Maksym Dzitsiuk and Toshiyuki Amano and Gudrun Klinker},
  Booktitle                = {14th {IEEE} International Symposium on Mixed and Augmented Reality, {ISMAR} 2015, FUkuoka, Japan, Sep. 29 - Oct. 3, 2015},
  Year                     = {2015},
  Pages                    = {43--48},
  Bibsource                = {dblp computer science bibliography, http://dblp.org},
  Crossref                 = {DBLP:conf/ismar/2015}
}
```
