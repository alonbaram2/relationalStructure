Surfing is a Matlab toolbox for surface-based information mapping of the cerebral cortex, intended for fMRI data. 

Nikolaas N. oosterhof [NNO]              <n.oosterhof@bangor.ac.uk>
           Tobias Wiestler [TW]          <tobias.wiestler.09@acl.ac.uk>
           Joern Diedrichsen [JD]        <j.diedrichsen@acl.ac.uk>

Version 0.3, Nov 2012.

Quick start:
- add surfing and surfing/misc folders to Matlab path
- download Peyre's "Fast Marching" toolbox from 
  http://www.ceremade.dauphine.fr/~peyre/matlab/fast-marching/content.html, and add the folder 
  and subfolders to the matlab path
- view the surfing_voxelselection.m function, or download the example data and run misc/example1.m

Surfing Toolbox 
Introduction
Surfing is a Matlab toolbox that implements voxel selection for surface-based information 
mapping of the cerebral cortex, intended for multivoxel pattern analysis (MVPA) of functional 
magnetic resonance imaging (fMRI) data. Unlike traditional information mapping (Kriegeskorte 
et al, 2006), where voxels are selected using a spherical mask, our toolbox uses a circular mask 
that follows the folded cortex and allows for grey matter masking. 
Briefly, a given set of nodes (N), the toolbox main function surfing_voxelselection selects voxels 
from functional data that are within a certain (geodesic) distance radius r (for each node 
separately), and returns the indices of the selected voxels. Alternatively, this function can also 
select a fixed number of k voxels, where the radius is variable across nodes. Importantly, 
functional data is not interpolated and/or projected on the surface (except for possible 
preprocessing such as slice time and motion correction).
To calculate the geodesic distance along the cortical surface Surfing uses the Fast Marching 
toolbox by Gabriel Peyre (2008), http://www.ceremade.dauphine.fr/~peyre/download/ 

Quick overview
What is surfing:
¥	voxel selection based on cortical surface reconstruction.
¥	support for Euclidian and geodesic distance metrics.
¥	selection of voxels across the brain based on either a fixed radius or a fixed number of 
voxels.

What surfing is not:
¥	no GUI: commands are entered in the Matlab command window, or in a matlab script.
¥	no full processing pipeline: preprocessing, surface reconstruction, classification, 
statistical analyses, and visualization are not part of the toolbox. Many other programs 
can be used for this, for example FSL, AFNI, SPM for preprocessing and regression; 
Freesurfer and Caret for surface reconstruction; several online toolboxes for 
classification and statistical analyses; and AFNI, Freesurfer and Caret for visualization. 
Hopefully one day our code may be compatible and/or integrated with with SPM, 
pyMVPA and/or the Princeton MVPA toolbox.
¥	no volume-based voxel selection (but you could easily write your own, and we may 
include this functionality in the future)

Requirements:
¥	Matlab, necessary to run the toolbox. 
¥	Freesurfer or Caret, for surface reconstruction. Freesurfer is preferred because it 
reconstructs inner (white) and outer (pial) surfaces, while Caret only provides a single 
surface (around layer 4). 
¥	Some familiarity with programming and fMRI analysis. 
Compatibility
We tested the toolbox under Matlab R2009b (7.9.0.529) for Windows and Matlab R2008a 
(7.6.0.324) on Mac OS 10.5 Leopard, using Freesurfer for surface reconstruction, AFNI and SPM 
for preprocessing and regression, and Caret and Freesurfer for visualization. 
Instructions
To find out more about the structure and functions in the Surfing toolbox, see below for a quick 
tutorial.
Current version
The present version is 0.1alpha, which is the very first release. In itÕs present form, all features are 
experimental and feedback and bug reports are welcomed. 
License
The toolbox is licensed under the permissive Òsimplified BSD licenseÓ. 

Copyright (c) 2009-2010 Nikolaas N. Oosterhof, Tobias Wiestler, Joern Diedrichsen.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted 
provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of 
conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of 
conditions and the following disclaimer in the documentation and/or other materials provided with 
the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY 
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Citation
If you use this toolbox for a scientific publication, please cite:
Oosterhof, N.N., Wiestler, T, Downing, P.E., & Diedrichsen, J. (in press). A comparison of volume-
based and surface-based information mapping. Neuroimage. Available online 
http://dx.doi.org/10.1016/j.neuroimage.2010.04.270
Contact
For any questions, comments, bug-fixes, feature-requests, compliments etcetera, please 
contact us at n.oosterhof [at]  bangor.ac.uk and tobias.wiestler.09 [at]  acl.ac.uk 


Tutorial 
Preparation
Download the toolbox from http://surfing.sourceforget.net, put the files in a suitable directory 
and add the main directory and (optionally) the misc directory to your path. Also, make sure 
you download the fast marching toolbox (google for: fast marching peyre), add those files to 
your path, and compile the mex functions if necessary (on Linux).
Input files
Before you can use the toolbox, you need reconstructed surfaces (inner and outer, i.e. white-
matter/grey matter and grey matter/pial) that are aligned to the functional volumes. This can 
be accomplished in many ways, and the authors used the following: 
¥	Usual preprocessing of functional data
¥	Surface reconstruction using freesurfer and one or more anatomical (T1) volumes
¥	Resample freesurfer surfaces to standard topology (using mapicosahedron)
¥	Find the (affine) transformation from the freesurfer output anatomical image and the 
functional volumes (possibly indirectly through , and apply this transformation to the 
resampled surfaces
¥	For group analyses, use surface in common Talairach space and average the coordinates 
across participants. Use of the standard topology (see above) ensures that data can be 
analyses node-wise without interpolation.
AFNI users can use the experimental script examples/afni_align_surfaces.sh.
2. Main function 
The core function of the toolbox is surfing_voxelselection.m, which uses the complete 
functionality of Surfing to select voxels along the cortical surface. A full example is given in 
examples/example1.m; in this section we briefly describe the most important steps.
The main function is called by
>>n2v=surfing_voxelselection(c1,c2,f,circledef,voldef,nodeidxs,linedef,distancemetric,progressstep)

This function takes two surfaces c1 and c2 with shared topology f, typically the pial and a 
surface of the white-matter boundary. Voldef defines the volume dimensions, nodeidxs 
contains the nodes indices for which a center is taken and voxels around it are selected, and 
circledef defines either a fixed radius r (form [r 0]), or a fixed number of voxels k (form [0 k]). 
The other arguments (linedef, distancemetric, progressstep) are optional and define how voxels 
are selected; linedef is described in more detail below. The output n2v contains the result, 
described at the end. 
The task at hand is now to find nodes that have a certain distance r from a center node and to 
collect the corresponding voxels for these nodes. To make sure that we collect all voxels 
between the to surfaces we define new coordinates between the layers. This is done by the 
function 
>>surfing_nodeidxs2coords(c1,c2,1:nverts,linedef);
which takes the coordinates of the two surfaces an index parameters that refer to all nodes and 
a matrix of three parameters that defines the starting and end point of each line that is ÒdrawnÓ 
between the corresponding nodes of the surfaces. A value of 0 corresponds to the first surface, 
a value of 1 corresponds to the second surface. The default is MN=0,MX=1.
Examples: 
[1,0.5,0.5] 	returns an intermediate surface in between the two input surfaces 
[10,0,1]    	returns ten coordinates per vertex along the line connecting the two surfaces    
[2,0,1]     	returns the same coordinates as the input surfaces
As an output the function will give back a 3xPxN array of coordinates (P steps, N nodes). Because 
ultimately we are interested in the corresponding voxels, we pass the coordinates to 
>>surfing_coords2linvoxelidxs(coords,voldef);
which takes the coordinates of interest and the volume information, and returns the linear 
voxel indices. Voxels outside the volume are set to NaN.
The output of surfing_voxelselection is an n2v cell that contains a mapping from node indices to 
linear voxel indices. The linear indices in each cell element refer to the grey matter surrounding 
the corresponding node. Use surfing_inds2subs (or the builtin ind2sub) to convert the linear to 
sub indices.
