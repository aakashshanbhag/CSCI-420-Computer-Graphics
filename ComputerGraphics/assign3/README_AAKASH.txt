CSCI 420 Computer Graphics
Assignment 3: Ray Tracer
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
Usage1: ./assign3 <scenefile> <softShadowFlag> <lightMotionFlag> <motionBlurFlag> <jpegname> <jpegFlag>

Flag=1, enable
Flag=0, disable

eg for soft shadows: ./assign3 screenfile.txt 1 0 0 img.jpg 1

eg for light motion: ./assign3 screenfile.txt 0 1 0 0 0 //No saving of frames.
eg for light motion: ./assign3 screenfile.txt 0 1 0 0 1 //Saving of frames.

eg for animation: ./assign3 screenfile.txt 1 1 1 0 0//No saving of frames.
eg for animation: ./assign3 screenfile.txt 1 1 1 0 1//Saving of frames.


eg for motion blur: ./assign3 screenfile.txt 0 0 1 0 0//No saving of frames.
eg for motion blur: ./assign3 screenfile.txt 0 0 1 0 1//Saving of frames.


Usage2: ./assign3 <scenefile> [jpegname]

------------------------------------------------------------------------------------------------

Three kinds of displays:
1. When motionBlur=0 and lightMotion=0 : A still image is displayed after ray tracing
2. When motionBlur=0 and lightMotion=1 : An animation is produced as the lights in the scene move
3. When motionBlur=1 and lightMotion=0 : An animation is produced showcasing the motion blur effect

Keys:	 
	  27 -> Exit 


------------------------------------------------------------------------------------------------

The Code follows the requirements specified with:

------------------------------------------------------------------------------------------------
Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                    YES

2) Ray tracing sphere                       YES

3) Triangle Phong Shading                   YES

4) Sphere Phong Shading                     YES

5) Shadows rays                             YES

6) Still images                             YES
   
7) Extra Credit  			YES
 
 ------------------------------------------------------------------------------------------------
 For motion blur and soft shadows and animation:

 >Utilised the accumulation buffer in all the above cases.
 >Moved the lights in the scene by an offset for the animation with refresh in the buffers and redisplay per frame.
 >For motion blur, blur updates saved in the accumulation buffer per frame and then translated the blur by a factor of 2 every frame.
 > Refreshed color and accumulation buffers for the animations and other features.

 >Tutorials from http://nehe.gamedev.net/ & Caleb Piercy's youtube channel helped in the animation and usage of the accumulation buffer. 

> Outputs for the soft shadows for all the model scenes are shown for fov 60 and 90.
  












	 

