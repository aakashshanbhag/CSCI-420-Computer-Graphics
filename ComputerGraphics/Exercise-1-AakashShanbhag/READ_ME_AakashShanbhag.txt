CSCI 420 Computer Graphics
Assignment 1: Height Fields
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
> ./assign1 SantaMonicaMountains-256.jpg spiral.jpg

The first argument is the image (SantaMonicaMountains-256.jpg) whose height field needs to be created and the the second argument produces the color map(spiral.jpg) for the height map image.

If the second argument size does not match with the height map image then heightmap and the color map are produced from the same image, which is the default case. The heightmap and color map extra credit has been attempted.

The Jpeg files attached are also for the above case with SantaMonicaMountains producing the height maps and the Spiral.jpg producing the color map for it. Total 300 images have been attached which display the results as expected.

------------------------------------------------------------------------------------------------

The Code follows the requirements specified with:

•Handle at least a 256x256 image for your height field at interactive frame rates (window size of 640x480). Height field manipulations should run smoothly.

•Be able to render the height field as points, lines("wireframe"), or solid triangles (with keys for the user to switch between the three).

Keys: a -> Renders points
	  s -> Renders wireframes
	  d -> Renders solid triangles
	  x -> Saves a single screenshot
	  q -> Saves 300 screenshots 

•Render as a perspective view, utilizing GL's depth buffer for hidden surface removal using the gluperspective.

•Use input from the mouse to spin the heightfield around using glRotate.

•Use input from the mouse to move the heightfield around using glTranslate.
Keys: SHIFT -> along with the mouse input.

•Use input from the mouse to change the dimensions of the heightfield using glScale.
Keys: y -> along with the mouse input since MAC doesnt register ALT and CTRL keys.

•Color the vertices using some smooth gradient with the GLsmooth function.




