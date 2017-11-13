CSCI 420 Computer Graphics
Assignment 2: Roller Coaster
Author: Aakash Shanbhag
USC ID:3205699915
------------------------------------------------------------------------------------------------
The code in assign1.cpp was built on MAC-OS Sierra with 
OpenGL Version: 2.1 INTEL-10.25.17
OpenGL Renderer: Intel(R) HD Graphics 6000
Shading Language Version: 1.20
------------------------------------------------------------------------------------------------
The Make file can be utilised to build this project with the commands for MAC.

> make
> ./assign2 track.txt

The first argument is the track file and for the this case we select rollercoaster.sp as the points we want to reconstruct the spline.

The Jpeg files attached are also with all the levels implemented at 15fps. Total 1000 images have been attached which display the results as expected.

------------------------------------------------------------------------------------------------

The Code follows the requirements specified with:
Complete all levels.
•Properly render Catmull-Rom splines to represent your track.
•Render a texture-mapped ground and sky.
•Render a rail cross-section.
•Move the camera at a reasonable speed in a continuous path and orientation along the coaster.
•Render the coaster in an interesting manner (good visibility, realism).

Keys:	  a -> Begins animation
	  x -> Saves a single screenshot
	  q -> Saves 1000 screenshots
	  l -> Display the texture mapped rail with support beams and cross sections
	  r -> Zoom in
	  d -> Move center to the right
	  g -> Move center to the left
	  f -> Zoom out
	  y -> Scale
	  27 -> Exit 
------------------------------------------------------------------------------------------------
For Extra credits
•Render double rail (like in real railroad tracks).
•Add OpenGL lighting to make your coaster look more realistic (SPOT light with light0 is used).

------------------------------------------------------------------------------------------------
•By default any spline can be drawn. Using rollercoaster.sp for optimal results since its not planar. 
•Important to note that for the animation of the camera, in order to avoid entering black regions line 659 ensures that we start from the 10th point once we reach the end of the spline(-50 points).
•The rollercoaster.sp gives the best result for level 5 with quads being used to render the output.(Click l key for creating the realistic spline.)
• Rotation and scaling can be done in the same manner as mentioned in the previous assignment.





	 

