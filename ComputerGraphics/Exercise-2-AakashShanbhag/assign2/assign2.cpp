/*
  CSCI 420 Computer Graphics
  Assignment 2: Roller Coaster
  Author: Aakash Shanbhag
  USC ID:3205699915
*/

#include <stdio.h>
#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "pic.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include "math.h"
#include <time.h>
using namespace std;

Pic * g_pHeightData=NULL;

/* represents one control points along the spline */
struct points 
{
   double x;
   double y;
   double z;
};

/* spline struct which contains how many control points, and an array of control points */
struct spline 
{
   int numControlpoints;
   struct points *points;
};

// Additional features for key press
bool ride=0;
bool displayrail=0;
// Setting up the screenshot counter and the total number of frames.
const int totalframes = 1001;
int screenshotCounter = 0;
bool Screenshots_flag = false;

// Global variable description
// Initilization of the mouse and the keyboard parameters.
int g_vMousePos[2] = {0, 0};
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;
// Control states for mouse control
typedef enum { ROTATE, TRANSLATE, SCALE } CONTROLSTATE;
CONTROLSTATE g_ControlState = ROTATE;

/* state of the world */
float g_vLandRotate[3] = {0.0, 0.0, 0.0};
float g_vLandTranslate[3] = {0.0, 0.0, 0.0};
float g_vLandScale[3] = {1.0, 1.0, 1.0};

/* Gluint objects needed for textures*/
GLuint texture_ground;
GLuint texture_sky;
GLuint texture_top;
GLuint railTexture;
GLuint texture_wood;

/* Images needed for textures*/
Pic *groundImage_1024;
Pic *skyImage_1024;
Pic *topImage_1024;
Pic *imageRail;
Pic *woodTexture;

/* lighting */
GLfloat light_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_position[] = { 85.0, 85.0, 30.0, 1.0 };

/* materials */
GLfloat mat_ambient[] = { 0.3, 0.3, 0.3, 1.0 };
GLfloat mat_diffuse[] = { 0.3, 0.3, 0.3, 1.0 };
GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat low_shininess[] = { 2.0 };

/* the spline array */
struct spline *g_Splines;

/* total number of splines */
int g_iNumOfSplines;

/* Set up the field of view. */
const float FOV = 60.0; 

/* Camera initial set up for gluLookAt */
double eyeX = 0.0, eyeY = 0.0, eyeZ = 5.0;
double centerX = 0.0, centerY = 0.0, centerZ = 0.0;
double upX = 0.0, upY = 1.0, upZ = 0.0;
int maxPointCount=0;

/*  Creating a display list for faster computation */
GLuint displayListSpline;

/* Creating level 1 splines along with tangents,normals and binormals */
spline level1Spline;
spline tangents;
spline normals;
spline binormals;
int pointsCount = 0;

// Calculation of the magnitude of the vector
double magnitude_vector(double x,double y,double z)
{ 	
 return sqrt((x*x) + (y*y) + (z*z));
}

int loadSplines(char *argv) 
{
  char *cName = (char *)malloc(128 * sizeof(char));
  FILE *fileList;
  FILE *fileSpline;
  int iType, i = 0, j, iLength;


  /* load the track file */
  fileList = fopen(argv, "r");
  
  if (fileList == NULL) {
    printf ("can't open file\n");
    exit(1);
  }
  
  /* stores the number of splines in a global variable */
  fscanf(fileList, "%d", &g_iNumOfSplines);

  g_Splines = (struct spline *)malloc(g_iNumOfSplines * sizeof(struct spline));

  /* reads through the spline files */
  for (j = 0; j < g_iNumOfSplines; j++) {
    i = 0;
    fscanf(fileList, "%s", cName);
    fileSpline = fopen(cName, "r");

    if (fileSpline == NULL) {
      printf ("can't open file\n");
      exit(1);
    }

    /* gets length for spline file */
    fscanf(fileSpline, "%d %d", &iLength, &iType);

    /* allocate memory for all the points */
    g_Splines[j].points = (struct points *)malloc(iLength * sizeof(struct points));
    g_Splines[j].numControlpoints = iLength;

    /* saves the data to the struct */
    while (fscanf(fileSpline, "%lf %lf %lf", 
	   &g_Splines[j].points[i].x, 
	   &g_Splines[j].points[i].y, 
	   &g_Splines[j].points[i].z) != EOF) {
      i++;
    }
  }

  free(cName);

  return 0;
}

/* Write a screenshot to the specified filename */
void saveScreenshot (char *filename)
{
  int i, j;
  Pic *in = NULL;

  if (filename == NULL)
    return;

  /* Allocate a picture buffer */
  in = pic_alloc(640, 480, 3, NULL);

  printf("File to save to: %s\n", filename);

  for (i=479; i>=0; i--) {
    glReadPixels(0, 479-i, 640, 1, GL_RGB, GL_UNSIGNED_BYTE,
                 &in->pix[i*in->nx*in->bpp]);
  }

  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);
}

/* Calculation of tangents, normals, binormals which are a orthogonal set up */ 
void calculateTangents()
{
 	// Calculation of tangents by taking the derivative of the u scalar.

 	tangents.numControlpoints=g_Splines[0].numControlpoints*20-40;
 	tangents.points= new points[tangents.numControlpoints];
 	maxPointCount=tangents.numControlpoints;

 	// Standard tension parameter set to 0.5
	double s = 0.5;

	// Catmull-Rom Spline matrix
	double basisMatrix[4][4] = {
		
		{ -s, (2 - s), (s - 2), s },
		{ (2 * s), (s - 3), (3 - (2 * s)), -s },
		{ -s, 0, s, 0 },
		{ 0, 1, 0, 0 }
								};

	// Constraint matrix = basis matrix * control matrix

	double constraint[4][3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// Total points loaded into the gsplines object
	int numpoints = g_Splines->numControlpoints;
	pointsCount=0;
	// Since we are taking 4 points at a time, need to remove 3 from the total count
	for (int m = 0; m < numpoints - 3; m++) 
	{
		// Create control matrix
		double controlMatrix[4][3] =
		{
			{ g_Splines[0].points[m].x , g_Splines[0].points[m].y, g_Splines[0].points[m].z },
			{ g_Splines[0].points[m + 1].x, g_Splines[0].points[m + 1].y, g_Splines[0].points[m+1].z },
			{ g_Splines[0].points[m + 2].x, g_Splines[0].points[m + 2].y, g_Splines[0].points[m+2].z },
			{ g_Splines[0].points[m + 3].x, g_Splines[0].points[m + 3].y, g_Splines[0].points[m+3].z }
		};
		
		// Clear constraint matrix after computation of each set
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 3; j++)
			 {
				constraint[i][j] = 0;
			 }
		}
		// Multiplying basis matrix and control matrix to obtain constraint matrix 
		for (int i = 0; i < 4; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				for (int k = 0; k < 4; k++)
				 {
					constraint[i][j] += basisMatrix[i][k] * controlMatrix[k][j];
				 }
			}
		}
		
		// Output co-ordinates of the Catmull-Rom Spline matrix
		double x = 0, y = 0, z = 0;
		
		// Using brute force for thr scalar u with increments of 0.05 
		for (double u = 0.05; u <= 1; u += 0.05) 
		{ 
			// Final matrix multipication and results obtained with scalar in consideration(3*u^2  2*u 1 0)
			x = (3*u*u*constraint[0][0]) + (2*u*constraint[1][0]) + (1 * constraint[2][0]);
			y = (3*u*u*constraint[0][1]) + (2*u*constraint[1][1]) + (1 * constraint[2][1]);
			z = (3*u*u*constraint[0][2]) + (2*u*constraint[1][2]) + (1 * constraint[2][2]);

			
			double mag=magnitude_vector(x,y,z);

			
			if(mag==0)
			{
				tangents.points[pointsCount].x =0.0;
				tangents.points[pointsCount].y =0.0;
				tangents.points[pointsCount].z =0.0;

			}
			else{
			// Storing globally as well after display list has a copy
			tangents.points[pointsCount].x = x/(1.0*mag);
			tangents.points[pointsCount].y = y/(1.0*mag);
			tangents.points[pointsCount].z = z/(1.0*mag);
			}
			pointsCount++;	
		}
		}
	// Reset points count at each stage for new spline
	pointsCount = 0;

	
}

// Calculation for 0th points
void InitialEstimatedNormals()
{
	normals.numControlpoints = tangents.numControlpoints;
	normals.points = new points[normals.numControlpoints];


	// Taking a arbitrary points(0,1,0) and (0,0,0) along with the tangent matrix points on the spline
	// Cross product of the vectors would yield a normal vector which is orthogonal to both. 
	// Hence since we 3 vertices easy to calculate the normal

	double normal_mag=magnitude_vector(0.0*(tangents.points[0].z-tangents.points[0].y),-1.0*(tangents.points[0].z-tangents.points[0].x),0.0*(tangents.points[0].y-tangents.points[0].x));


	if(normal_mag==0)
	{
	normals.points[0].x=0.0;
	normals.points[0].y=0.0;
	normals.points[0].z=0.0;	
	}
	else
	{
		//Normalised normal vector obtained after dividing the magnitude of the vector
	normals.points[0].x=0.0*(tangents.points[0].z-tangents.points[0].y)/(normal_mag);
	normals.points[0].y=-1.0*(tangents.points[0].z-tangents.points[0].x)/(normal_mag);
	normals.points[0].z=0.0*(tangents.points[0].y-tangents.points[0].x)/(normal_mag);

	}
		
}

// Calculation of the binormal at the 0th points using the 0th points normal.
void InitialEstimatedBinormals()
{
	binormals.numControlpoints=tangents.numControlpoints;
	binormals.points=new points[binormals.numControlpoints];

	// Binormals are obtained by the cross product of the normal and the tangent vectors at each points on the spline
	// Cross product values need to be normalised as well
	
	double binorm_x =(tangents.points[0].y * normals.points[0].z)-(tangents.points[0].z * normals.points[0].y);
	double binorm_y =(tangents.points[0].z * normals.points[0].x)-(tangents.points[0].x * normals.points[0].z);
	double binorm_z =(tangents.points[0].x * normals.points[0].y)-(tangents.points[0].y * normals.points[0].x);

	// Calculation of the normal.
	double binorm_mag=magnitude_vector(binorm_x,binorm_y,binorm_z);
	
	if(binorm_mag==0)
	{
	binormals.points[0].x=0.0;
	binormals.points[0].y=0.0;
	binormals.points[0].z=0.0;	
	}
	else
	{
	binormals.points[0].x= -1.0*binorm_x/(1.0*binorm_mag) ;
	binormals.points[0].y= -1.0*binorm_y/(1.0*binorm_mag);
	binormals.points[0].z= -1.0*binorm_z/(1.0*binorm_mag);
	}

}

void NormalsBinormals()
{
	InitialEstimatedNormals();
	InitialEstimatedBinormals();
	

	// Calculation of normals and the binormals from the 0th points (Normal=Binormal cross Tangent)

	for (int i = 1; i < level1Spline.numControlpoints; i++) {

	
		double norm_x = (binormals.points[i - 1].y * tangents.points[i].z) - (binormals.points[i - 1].z * tangents.points[i].y);
		double norm_y = (binormals.points[i - 1].z * tangents.points[i].x) - (binormals.points[i - 1].x * tangents.points[i].z);
		double norm_z = (binormals.points[i - 1].x * tangents.points[i].y) - (binormals.points[i - 1].y * tangents.points[i].x);

		// Normalised values by taking magnitude of the vector

		double normalMagnitude = magnitude_vector((norm_x),(norm_y),(norm_z));
		
		if(normalMagnitude==0)
		{
			norm_x=0.0;
			norm_y=0.0;
			norm_z=0.0;
			
		}
		else
		{
		// normalize the normal vector
		norm_x = 1.0*norm_x /(normalMagnitude);
		norm_y = 1.0*norm_y / normalMagnitude;
		norm_z = 1.0*norm_z / normalMagnitude;
		}
		
		// Final normals
		normals.points[i].x = norm_x;
		normals.points[i].y = norm_y;
		normals.points[i].z = norm_z;
	
		// Tangent cross normal= binormal vector

		double binorm_x = (tangents.points[i].y * normals.points[i].z ) - (tangents.points[i].z * normals.points[i].y);
		double binorm_y = (tangents.points[i].z * normals.points[i].x) - (tangents.points[i].x * normals.points[i].z );
		double binorm_z = (tangents.points[i].x * normals.points[i].y) - (tangents.points[i].y * normals.points[i].x);

		// Calculate the magnitude of the binormal vector for normalization
		double binormalMagnitude = magnitude_vector( (binorm_x) , (binorm_y) , (binorm_z));

		if(binormalMagnitude==0)
		{
			binorm_x=0.0;
			binorm_y=0.0;
			binorm_z=0.0;
		}
		else{
		// normalize the binormal vector
	    binorm_x = binorm_x / binormalMagnitude;
		binorm_y = binorm_y / binormalMagnitude;
		binorm_z = binorm_z / binormalMagnitude;
		}
		// Final binormals for vectors
		binormals.points[i].x = binorm_x;
		binormals.points[i].y = binorm_y;
		binormals.points[i].z = binorm_z;

	}
}

void initSpline() 
{
	// Create spline and store it in display list

	level1Spline.numControlpoints = g_Splines[0].numControlpoints * 20-40;
	// Storing all the control points found
	level1Spline.points = new points[level1Spline.numControlpoints];

	// Display the number of points in the input an the total points created for generation of spline.
	cout<< "The number of points found: "<<g_Splines[0].numControlpoints<<endl;
	cout<< "The number of points used for creating spline: "<<level1Spline.numControlpoints<<endl;
	
	// Standard tension parameter set to 0.5
	double s = 0.5;

	// Catmull-Rom Spline matrix
	double basisMatrix[4][4] = {
		
		{ -s, (2 - s), (s - 2), s },
		{ (2 * s), (s - 3), (3 - (2 * s)), -s },
		{ -s, 0, s, 0 },
		{ 0, 1, 0, 0 }
								};

	// Constraint matrix = basis matrix * control matrix

	double constraint[4][3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// Total points loaded into the gsplines object
	int numpoints = g_Splines->numControlpoints;

	// Since we are taking 4 points at a time, need to remove 3 from the total count
	for (int m = 0; m < numpoints - 3; m++) 
	{
		// Create control matrix
		double controlMatrix[4][3] =
		{
			{ g_Splines[0].points[m].x , g_Splines[0].points[m].y, g_Splines[0].points[m].z },
			{ g_Splines[0].points[m + 1].x, g_Splines[0].points[m + 1].y, g_Splines[0].points[m+1].z },
			{ g_Splines[0].points[m + 2].x, g_Splines[0].points[m + 2].y, g_Splines[0].points[m+2].z },
			{ g_Splines[0].points[m + 3].x, g_Splines[0].points[m + 3].y, g_Splines[0].points[m+3].z }
		};
		
		// Clear constraint matrix after computation of each set
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 3; j++)
			 {
				constraint[i][j] = 0;
			 }
		}
		// Multiplying basis matrix and control matrix to obtain constraint matrix 
		for (int i = 0; i < 4; i++) 
		{
			for (int j = 0; j < 3; j++) 
			{
				for (int k = 0; k < 4; k++)
				 {
					constraint[i][j] += basisMatrix[i][k] * controlMatrix[k][j];
				 }
			}
		}
		
		// Output co-ordinates of the Catmull-Rom Spline matrix
		double x = 0, y = 0, z = 0;
		
		// Using brute force for thr scalar u with increments of 0.05 
		for (double u = 0.05; u <= 1; u += 0.05) 
		{ 
			// Final matrix multipication and results obtained
			x = (u*u*u*constraint[0][0]) + (u*u*constraint[1][0]) + (u * constraint[2][0])+ constraint[3][0];
			y = (u*u*u*constraint[0][1]) + (u*u*constraint[1][1]) + (u * constraint[2][1])+ constraint[3][1];
			z = (u*u*u*constraint[0][2]) + (u*u*constraint[1][2]) + (u * constraint[2][2])+ constraint[3][2];

			// store the vertices in display list
			glVertex3f(x, y, z);
			pointsCount++;

			// Storing globally as well after display list has a copy
			level1Spline.points[pointsCount].x = x;
			level1Spline.points[pointsCount].y = y;
			level1Spline.points[pointsCount].z = z;
			
		}
	
	}
	// Reset points count at each stage for new spline
	pointsCount = 0;
}

void initLevel2Ground() 
{
	groundImage_1024 = jpeg_read("background.jpg", NULL);

	// Placeholder creation
	glGenTextures(1, &texture_ground);

	// Activate textures
	glBindTexture(GL_TEXTURE_2D, texture_ground);

	// Repeat Textures by wrapping
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	
	// Linear filter used for magnification and vice versa
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Image loaded into the placeholder
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, groundImage_1024->pix);
}

void createHeightField(Pic * img)
{

   glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
   glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
   glEnable(GL_TEXTURE_GEN_S);
   glEnable(GL_TEXTURE_GEN_T);
   //glEnable(GL_TEXTURE_2D);


  int img_height= img->ny;
  int img_width= img->nx;
  float Scale =0.125* img_height/ 100.0;// scale corresponding to the max height. 

  // Generate our vertices--we should later use a triangle strip, that will make this so much more efficient. We center the heightmap at the origin.
 
 	// Fixed lights and no modulation based on light
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	// Binding the texture with the object
	glBindTexture(GL_TEXTURE_2D, texture_ground);
	
	// Enabling on for texture mapping
	glEnable(GL_TEXTURE_2D);


  for (int i = 0; i < img_height - 1; i++) {
     glBegin(GL_TRIANGLE_STRIP);
        for (int j = 0; j < img_width - 1; j++) {
 
        float intensity_lower_left=PIC_PIXEL(img,j, i, 0);
        float intensity_upper_left=PIC_PIXEL(img,j, i + 1 , 0);

        float intensity_lower_left_c=PIC_PIXEL(img,j, i, 0)/255.0;
        float intensity_upper_left_c=PIC_PIXEL(img,j, i + 1 , 0)/255.0;

        glColor4f(intensity_lower_left_c,intensity_lower_left_c,intensity_lower_left_c,1.0);
        glVertex3f(j,i,GLfloat(Scale * intensity_lower_left));
        
        glColor4f(intensity_upper_left_c,intensity_upper_left_c,intensity_upper_left_c,1.0);
        glVertex3f(j,i+1, GLfloat(Scale * intensity_upper_left));

        }
     glEnd();

    }
    glDisable(GL_TEXTURE_2D);
} 
 
void initSky()
{
	skyImage_1024=jpeg_read("sky.jpg",NULL);

	// Placeholder creation
	glGenTextures(1, &texture_sky);

	// Activate textures
	glBindTexture(GL_TEXTURE_2D, texture_sky);

	// Repeat Textures by wrapping
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	
	// Linear filter used for magnification and vice versa
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Image loaded into the placeholder
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, skyImage_1024->pix);


	// Top Image also needs the same initialisation as all the above

	topImage_1024=jpeg_read("sky.jpg",NULL);

	// Placeholder creation
	glGenTextures(1, &texture_top);

	// Activate textures
	glBindTexture(GL_TEXTURE_2D, texture_top);

	// Repeat Textures by wrapping
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	
	// Linear filter used for magnification and vice versa
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Image loaded into the placeholder
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, topImage_1024->pix);
}

void initLights()
{
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_ambient);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
}

void initMaterials()
{
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, low_shininess);
}

void initRailTexture() 
{
	imageRail = jpeg_read("metalTexture.jpg", NULL);

	// create placeholder for texture
	glGenTextures(1, &railTexture);

	// make texture active
	glBindTexture(GL_TEXTURE_2D, railTexture);

	// texture parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	// use linear filter both for magnification and minification
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// load image data stored at pointer "imageRail" into currently active texture
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, imageRail->pix);
}

void initWoodTexture()
{

	woodTexture = jpeg_read("wood.jpg", NULL);

	// create placeholder for texture
	glGenTextures(1, &texture_wood);

	// make texture active
	glBindTexture(GL_TEXTURE_2D, texture_wood);

	// texture parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	// use linear filter both for magnification and minification
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// load image data stored at pointer "imageRail" into currently active texture
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, woodTexture->pix);
}	


int pointindex=0;
clock_t t;
void cameraSetup()
{
	if(pointindex==maxPointCount-50)
	pointindex=10;
	clock_t newT = clock() - t;
	pointindex++;
	//float sec = ((float)newT) / CLOCKS_PER_SEC;

	eyeX = level1Spline.points[pointindex].x+normals.points[pointindex].x/4 ;
	eyeY = level1Spline.points[pointindex].y+normals.points[pointindex].y/4 ;
	eyeZ = level1Spline.points[pointindex].z+0.25;


	centerX=eyeX+tangents.points[pointindex].x;
	centerY=eyeY+tangents.points[pointindex].y;
	centerZ=eyeZ+tangents.points[pointindex].z;

	if(binormals.points[pointindex].x==0 && binormals.points[pointindex].y==0 && binormals.points[pointindex].z==0)
	{
	upX=0.0;
	upY=-1.0;
	upZ=0.0;
	}
	else
	{	
	upX = -binormals.points[pointindex].x;
	upY = -binormals.points[pointindex].y;
	upZ = -binormals.points[pointindex].z;
	}

	t=clock();
}

void displaySpline() 
{
	// Rendering the Spline white in the background
	glColor3f(1.0, 1.0, 1.0);

	if(displayrail==0)
	{

	/* Draw line through the vertices obtained.*/
	
	glBegin(GL_LINES);
	for (int i = 0; i < level1Spline.numControlpoints-1; i++) {
			
			
			glVertex3f(level1Spline.points[i].x, level1Spline.points[i].y, level1Spline.points[i].z);
		
			glVertex3f(level1Spline.points[i+1].x, level1Spline.points[i+1].y, level1Spline.points[i+1].z);
		}
	glEnd();

	/* Create cross sections */
	glBegin(GL_LINES);

	for (int i = 0; i < normals.numControlpoints; i++) {
		glVertex3f(level1Spline.points[i].x, level1Spline.points[i].y, level1Spline.points[i].z);
		glVertex3f(level1Spline.points[i].x + normals.points[i].x / 2,
			level1Spline.points[i].y + normals.points[i].y / 2,
			level1Spline.points[i].z + normals.points[i].z / 2);
	}
	
	/* Create support beams */
	
	glBegin(GL_LINES);
	for (int i = 0; i < normals.numControlpoints; i++) {
		if (i % 20 == 0) {
			glVertex3f(level1Spline.points[i].x, level1Spline.points[i].y, level1Spline.points[i].z);
			glVertex3f(level1Spline.points[i].x, level1Spline.points[i].y, -1.0);

		}
	}
	glEnd();
	
	}
	else
	{	// Creating a display list for direct rendering on the GPU for repeated calls.
		glCallList(displayListSpline);
	}
	
	// Flush the display list for refresh of the list. 
	glFlush();
}

/* Display module for Ground */
void displayGround() 
{

	// Fixed lights and no modulation based on light
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	// Binding the texture with the object
	glBindTexture(GL_TEXTURE_2D, texture_ground);
	
	// Enabling on for texture mapping
	glEnable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);

	glTexCoord2f(0.0, 0.0);
	glVertex3f(-30.0, -30.0, -1.0);
	glTexCoord2f(0.0, 50.0);
	glVertex3f(-30.0, 60.0, -1.0);
	glTexCoord2f(50.0, 0.0);
	glVertex3f(60.0, 60.0, -1.0);
	glTexCoord2f(50.0, 50.0);
	glVertex3f(60.0, -30.0, -1.0);

	glEnd();

	glDisable(GL_TEXTURE_2D);
}

void displaySkybox() 
{
	// Fixed lights and no modulation based on light
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	// Binding the texture with the object
	glBindTexture(GL_TEXTURE_2D, texture_sky);
	// Enabling on for texture mapping
	glEnable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);

	// Creation of box with 4 sides having sky text with essential locations 
	// x variations:-30,60
	// y variations:-30,60
	// z variations:-1,60
	// Side 1 : 
	glTexCoord2f(1.0, 1.0);
	glVertex3f(-30.0, -30.0, -1.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(60.0, -30.0, -1.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(60.0, -30.0, 50.0);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(-30.0, -30.0, 50.0);

	// Side 2 :  
	glTexCoord2f(1.0, 1.0);
	glVertex3f(-30.0, -30.0, -1.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(-30.0, 60.0, -1.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(-30.0, 60.0, 50.0);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(-30.0, -30.0, 50.0);

	// Side 3 : 
	glTexCoord2f(1.0, 1.0);
	glVertex3f(-30.0, 60.0, -1.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(60.0, 60.0, -1.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(60.0, 60.0, 50.0);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(-30.0, 60.0, 50.0);

	// Side 4 :
	glTexCoord2f(1.0, 1.0);
	glVertex3f(60.0, 60.0, -1.0);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(60.0, -30.0, -1.0);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(60.0, -30.0, 50.0);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(60.0, 60.0, 50.0);

	glEnd();

	glDisable(GL_TEXTURE_2D);

	// Binding the texture with the object
	glBindTexture(GL_TEXTURE_2D, texture_top);
	// Enabling on for texture mapping
	glEnable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);

	// top of skybox
	// Important to map z accurately(range of z) for perfect top box.
	glTexCoord2f(1.0, 1.0);
	glVertex3f(60.0, 60.0, 49.9);
	glTexCoord2f(0.0, 1.0);
	glVertex3f(60.0, -30.0, 49.9);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(-30.0, -30.0, 49.9);
	glTexCoord2f(1.0, 0.0);
	glVertex3f(-30.0, 60.0, 49.9);

	glEnd();
	glDisable(GL_TEXTURE_2D);
}

void keyboardFunc(unsigned char key, int x, int y) 
{
  switch (key) {
    
    case 27: // ESC key
        exit(0); // exit the program
        break;

    case 'x':
      // take a screenshot
        saveScreenshot("screenshot.jpg");
        break;

    case 'y':
    // Issue with CTRL and ALT keys dealt with y key.
        g_ControlState = SCALE;
        break;

    // Start the roller coaster by camera movement. 
    case 'a':
    
    	// Ensure the fact that we are still on the coaster while animating.
    	g_vLandRotate[0]=0.0;
    	g_vLandRotate[1]=0.0;
    	g_vLandRotate[2]=0.0;
    	g_vLandTranslate[0]=0.0;
    	g_vLandTranslate[1]=0.0;
    	g_vLandTranslate[2]=0.0;

    	ride=!ride;
    	
    	break;

    case 'l':
    	displayrail=!displayrail;
    	
    	break;

    // Zoom in Zoom out
    case 'r':
    	eyeZ -= 1.0;
    	cout<<"Zoom in"<<endl;
    	break;

    case 'd':
    	centerZ -= 1.0;
    	cout<<"Move center to the right"<<endl;
    	break;

    //	Look left right
    case 'g':
     	centerZ += 1.0;
     	cout<<"Move center to the left"<<endl;
    	break;

    case 'f':
     	eyeZ += 1.0;
     	cout<<"Zoom out"<<endl;
    	break;

    case 'q':
	    // Screenshot capture to commence
	    Screenshots_flag = true;
	    std::cout << "Starting animation!." << std::endl;
	    break;

  }
}

/* converts mouse drags into information about 
rotation/translation/scaling */
void mousedrag(int x, int y)
{
  int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};

  
  switch (g_ControlState)
  {
    case TRANSLATE:  
      if (g_iLeftMouseButton)
      {
        g_vLandTranslate[0] += vMouseDelta[0]*0.5;
        g_vLandTranslate[1] -= vMouseDelta[1]*0.5;
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandTranslate[2] += vMouseDelta[1]*0.5;
      }

      break;
    case ROTATE:
      if (g_iLeftMouseButton)
      {
        g_vLandRotate[0] += vMouseDelta[1];
        g_vLandRotate[1] += vMouseDelta[0];
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandRotate[2] += vMouseDelta[1];
      }
      
      break;
    case SCALE:
      if (g_iLeftMouseButton)
      {
        g_vLandScale[0] *= 1.0+vMouseDelta[0]*0.01;
        g_vLandScale[1] *= 1.0-vMouseDelta[1]*0.01;
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandScale[2] *= 1.0-vMouseDelta[1]*0.01;
      }
   
     
      break;
  }
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mouseidle(int x, int y)
{
  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void mousebutton(int button, int state, int x, int y){


  switch (button)
  {
    case GLUT_LEFT_BUTTON:
      g_iLeftMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_MIDDLE_BUTTON:
      g_iMiddleMouseButton = (state==GLUT_DOWN);
      break;
    case GLUT_RIGHT_BUTTON:
      g_iRightMouseButton = (state==GLUT_DOWN);
      break;
  }
 
  switch(glutGetModifiers())
  {
    case GLUT_ACTIVE_SHIFT:
      g_ControlState = TRANSLATE;
      
      break;
    case GLUT_ACTIVE_ALT:
      g_ControlState = SCALE;
      // Issue with MAC, hence bind with a key (y)
      break;
    default:
      g_ControlState = ROTATE;
      break;
  }

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}

void reshape(int w, int h)
{
    // setup image size
    glViewport(0,0, w, h); 
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // setup camera
    gluPerspective (FOV, float(w)/float(h), 0.01, 1000.0);
    glMatrixMode(GL_MODELVIEW);
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f(0.0, 0.0, 1.0);
	glShadeModel(GL_SMOOTH);
	glLoadIdentity();


	if(ride)
	cameraSetup();
	/* camera view */
	gluLookAt(eyeX, eyeY, eyeZ,
		centerX, centerY, centerZ,
		upX, upY, upZ);

	// Write a basic model to read and apply the spline to points.
	glPushMatrix();
	glTranslatef(1.0*g_vLandTranslate[0],1.0*g_vLandTranslate[1],1.0*g_vLandTranslate[2]);
    // Scaling on mouse commands.
    glScalef(g_vLandScale[0],g_vLandScale[1],g_vLandScale[2]);
    // Rotate using mouse callback
    glRotatef(g_vLandRotate[0], 1.0, 0.0,0.0);
    glRotatef(g_vLandRotate[1], 0.0, 1.0,0.0);
    glRotatef(g_vLandRotate[2], 0.0, 0.0,1.0);
    // Write a basic model to read and apply the spline to points.
    

    displaySkybox();
    displaySpline();
    displayGround();

  	/*
    glPushMatrix();
    glTranslatef(0.0,0.0,0.0);
    glScalef(0.02,0.02,0.05);
    createHeightField(g_pHeightData);
    glPopMatrix();
    */
   
    glPopMatrix();

	/* needed for double buffering*/
	glutSwapBuffers();
}

 /* Creation of idle function to post redisplay every frame */
void doIdle()
{
	// Screen shot files to be created for the animation.
  	if (screenshotCounter < totalframes && Screenshots_flag) {
    std::stringstream ss;
    ss << "../assign2/TrialAnimation/trial-" << screenshotCounter << ".jpg";
    saveScreenshot ((char*)ss.str().data());
    screenshotCounter++;
  }
	/* make the screen update */
	glutPostRedisplay();
}

/* OpenGL init */
void myInit() 
{

	t=clock();
	// Initialisation for the scene with right materials and a single source of light
	initLights();
	initMaterials();

	//g_pHeightData=jpeg_read("SantaMonicaMountains-768.jpg",NULL);

	// Single source of illumination.
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	// Initial conditions for rendering of the different levels
	initLevel2Ground();
	initSky();
	initRailTexture();
	initWoodTexture();
	initSpline();
	calculateTangents();
	NormalsBinormals();

	/* Intialise the display list */
	displayListSpline = glGenLists(1);
	glNewList(displayListSpline, GL_COMPILE);
	
	// modulate texture with lighting
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glBindTexture(GL_TEXTURE_2D, railTexture);
	// turn on texture mapping
	glEnable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);
	double distance_factor = 0.05;
	for (int i = 0; i < level1Spline.numControlpoints; i++)
	{
		double X = level1Spline.points[i].x;
		double Y = level1Spline.points[i].y;
		double Z = level1Spline.points[i].z;
		
		double X2 = level1Spline.points[i + 1].x;
		double Y2 = level1Spline.points[i + 1].y;
		double Z2 = level1Spline.points[i + 1].z;
		
		double X1_rail1 = X - distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y1_rail1 = Y - distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z1_rail1 = Z - distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X2_rail1 = X + distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y2_rail1 = Y + distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z2_rail1 = Z + distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X3_rail1 = X + distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y3_rail1 = Y + distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z3_rail1 = Z + distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X4_rail1 = X - distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y4_rail1 = Y - distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z4_rail1 = Z - distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X5_rail1 = X2 - distance_factor*(normals.points[i+1].x - binormals.points[i + 1].x);
		double Y5_rail1 = Y2 - distance_factor*(normals.points[i+1].y - binormals.points[i + 1].y);
		double Z5_rail1 = Z2 - distance_factor*(normals.points[i+1].z - binormals.points[i + 1].z);

		double X6_rail1 = X2 + distance_factor*(normals.points[i+1].x + binormals.points[i + 1].x);
		double Y6_rail1 = Y2 + distance_factor*(normals.points[i+1].y + binormals.points[i + 1].y);
		double Z6_rail1 = Z2 + distance_factor*(normals.points[i+1].z + binormals.points[i + 1].z);

		double X7_rail1 = X2 + distance_factor*(normals.points[i+1].x - binormals.points[i + 1].x);
		double Y7_rail1 = Y2 + distance_factor*(normals.points[i+1].y - binormals.points[i + 1].y);
		double Z7_rail1 = Z2 + distance_factor*(normals.points[i+1].z - binormals.points[i + 1].z);

		double X8_rail1 = X2 - distance_factor*(normals.points[i+1].x + binormals.points[i + 1].x);
		double Y8_rail1 = Y2 - distance_factor*(normals.points[i+1].y + binormals.points[i + 1].y);
		double Z8_rail1 = Z2 - distance_factor*(normals.points[i+1].z + binormals.points[i + 1].z);

		//Rail1

		glNormal3f(-normals.points[i].x, -normals.points[i].y, -normals.points[i].z);
		// Right face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail1, Y1_rail1, Z1_rail1);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X2_rail1, Y2_rail1, Z2_rail1);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X6_rail1, Y6_rail1, Z6_rail1);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail1, Y5_rail1, Z5_rail1);

		glNormal3f(-binormals.points[i].x, -binormals.points[i].y, -binormals.points[i].z);
		// Top face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X2_rail1, Y2_rail1, Z2_rail1);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X3_rail1, Y3_rail1, Z3_rail1);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X7_rail1, Y7_rail1, Z7_rail1);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X6_rail1, Y6_rail1, Z6_rail1);

		glNormal3f(normals.points[i].x, normals.points[i].y, normals.points[i].z);
		// Left face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X3_rail1, Y3_rail1, Z3_rail1);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail1, Y4_rail1, Z4_rail1);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail1, Y8_rail1, Z8_rail1);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X7_rail1, Y7_rail1, Z7_rail1);

		glNormal3f(binormals.points[i].x, binormals.points[i].y, binormals.points[i].z);
		// Bottom face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail1, Y1_rail1, Z1_rail1);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail1, Y4_rail1, Z4_rail1);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail1, Y8_rail1, Z8_rail1);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail1, Y5_rail1, Z5_rail1);


		//Rail 2 using the coordinates.
		double X1_rail2 = X - distance_factor*(-8*normals.points[i].x - binormals.points[i].x);
		double Y1_rail2 = Y - distance_factor*(-8*normals.points[i].y - binormals.points[i].y);
		double Z1_rail2 = Z - distance_factor*(-8*normals.points[i].z - binormals.points[i].z);

		double X2_rail2 = X + distance_factor*(10*normals.points[i].x + binormals.points[i].x);
		double Y2_rail2 = Y + distance_factor*(10*normals.points[i].y + binormals.points[i].y);
		double Z2_rail2 = Z + distance_factor*(10*normals.points[i].z + binormals.points[i].z);

		double X3_rail2 = X + distance_factor*(10*normals.points[i].x - binormals.points[i].x);
		double Y3_rail2 = Y + distance_factor*(10*normals.points[i].y - binormals.points[i].y);
		double Z3_rail2 = Z + distance_factor*(10*normals.points[i].z - binormals.points[i].z);

		double X4_rail2 = X - distance_factor*(-8*normals.points[i].x + binormals.points[i].x);
		double Y4_rail2 = Y - distance_factor*(-8*normals.points[i].y + binormals.points[i].y);
		double Z4_rail2 = Z - distance_factor*(-8*normals.points[i].z + binormals.points[i].z);

		double X5_rail2 = X2 - distance_factor*(-8*normals.points[i+1].x - binormals.points[i + 1].x);
		double Y5_rail2 = Y2 - distance_factor*(-8*normals.points[i+1].y - binormals.points[i + 1].y);
		double Z5_rail2 = Z2 - distance_factor*(-8*normals.points[i+1].z - binormals.points[i + 1].z);

		double X6_rail2 = X2 + distance_factor*(10*normals.points[i+1].x + binormals.points[i + 1].x);
		double Y6_rail2 = Y2 + distance_factor*(10*normals.points[i+1].y + binormals.points[i + 1].y);
		double Z6_rail2 = Z2 + distance_factor*(10*normals.points[i+1].z + binormals.points[i + 1].z);

		double X7_rail2 = X2 + distance_factor*(10*normals.points[i+1].x - binormals.points[i + 1].x);
		double Y7_rail2 = Y2 + distance_factor*(10*normals.points[i+1].y - binormals.points[i + 1].y);
		double Z7_rail2 = Z2 + distance_factor*(10*normals.points[i+1].z - binormals.points[i + 1].z);

		double X8_rail2 = X2 - distance_factor*(-8*normals.points[i+1].x + binormals.points[i + 1].x);
		double Y8_rail2 = Y2 - distance_factor*(-8*normals.points[i+1].y + binormals.points[i + 1].y);
		double Z8_rail2 = Z2 - distance_factor*(-8*normals.points[i+1].z + binormals.points[i + 1].z);

		glNormal3f(-normals.points[i].x, -normals.points[i].y, -normals.points[i].z);
		// Right face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail2, Y1_rail2, Z1_rail2);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X2_rail2, Y2_rail2, Z2_rail2);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X6_rail2, Y6_rail2, Z6_rail2);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail2, Y5_rail2, Z5_rail2);

		glNormal3f(-binormals.points[i].x, -binormals.points[i].y, -binormals.points[i].z);
		// Top face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail2, Y1_rail2, Z2_rail2);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X3_rail2, Y3_rail2, Z3_rail2);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X7_rail2, Y7_rail2, Z7_rail2);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X6_rail2, Y6_rail2, Z6_rail2);

		glNormal3f(normals.points[i].x, normals.points[i].y, normals.points[i].z);
		// Left face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X3_rail2, Y3_rail2, Z3_rail2);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail2, Y4_rail2, Z4_rail2);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail2, Y8_rail2, Z8_rail2);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X7_rail2, Y7_rail2, Z7_rail2);

		glNormal3f(binormals.points[i].x, binormals.points[i].y, binormals.points[i].z);
		// Bottom face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail2, Y1_rail2, Z1_rail2);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail2, Y4_rail2, Z4_rail2);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail2, Y8_rail2, Z8_rail2);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail2, Y5_rail2, Z5_rail2);


	}
	
	glEnd();


	// Guard rail1 and guard rail 2.
	glBegin(GL_QUADS);

	for (int i = 0; i < normals.numControlpoints; i++) {
		if (i % 20 == 0) {

		//guide rail 1

		double Xg = level1Spline.points[i].x;
		double Yg = level1Spline.points[i].y;
		double Zg = level1Spline.points[i].z;
		
		double X2g = level1Spline.points[i].x;
		double Y2g = level1Spline.points[i].y;
		double Z2g = -1.0;

		double X1_rail1g = Xg - distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y1_rail1g = Yg - distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z1_rail1g = Zg - distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X2_rail1g = Xg + distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y2_rail1g = Yg + distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z2_rail1g = Zg + distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X3_rail1g = Xg + distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y3_rail1g = Yg + distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z3_rail1g = Zg + distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X4_rail1g = Xg - distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y4_rail1g = Yg - distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z4_rail1g = Zg - distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X5_rail1g = X2g - distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y5_rail1g = Y2g - distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z5_rail1g = Z2g - distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X6_rail1g = X2g + distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y6_rail1g = Y2g + distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z6_rail1g = Z2g + distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X7_rail1g = X2g + distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y7_rail1g = Y2g + distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z7_rail1g = Z2g + distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X8_rail1g = X2g - distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y8_rail1g = Y2g - distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z8_rail1g = Z2g - distance_factor*(normals.points[i].z + binormals.points[i].z);


		//Rail1g

		glNormal3f(-normals.points[i].x, -normals.points[i].y, -normals.points[i].z);
		// Right face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail1g, Y1_rail1g, Z1_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X2_rail1g, Y2_rail1g, Z2_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X6_rail1g, Y6_rail1g, Z6_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail1g, Y5_rail1g, Z5_rail1g);

		glNormal3f(-binormals.points[i].x, -binormals.points[i].y, -binormals.points[i].z);
		// Top face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X2_rail1g, Y2_rail1g, Z2_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X3_rail1g, Y3_rail1g, Z3_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X7_rail1g, Y7_rail1g, Z7_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X6_rail1g, Y6_rail1g, Z6_rail1g);

		glNormal3f(normals.points[i].x, normals.points[i].y, normals.points[i].z);
		// Left face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X3_rail1g, Y3_rail1g, Z3_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail1g, Y4_rail1g, Z4_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail1g, Y8_rail1g, Z8_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X7_rail1g, Y7_rail1g, Z7_rail1g);

		glNormal3f(binormals.points[i].x, binormals.points[i].y, binormals.points[i].z);
		// Bottom face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail1g, Y1_rail1g, Z1_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail1g, Y4_rail1g, Z4_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail1g, Y8_rail1g, Z8_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail1g, Y5_rail1g, Z5_rail1g);

		// Rail2
		double Xg2 = level1Spline.points[i].x ;
		double Yg2 = level1Spline.points[i].y ;
		double Zg2 = level1Spline.points[i].z ;
		
		double X2g2 = level1Spline.points[i].x ;
		double Y2g2 = level1Spline.points[i].y ;
		double Z2g2 = -1.0;


		double X1_rail2g = Xg2 - distance_factor*(-8*normals.points[i].x - binormals.points[i].x);
		double Y1_rail2g = Yg2 - distance_factor*(-8*normals.points[i].y - binormals.points[i].y);
		double Z1_rail2g = Zg2 - distance_factor*(-8*normals.points[i].z - binormals.points[i].z);

		double X2_rail2g = Xg2 + distance_factor*(10*normals.points[i].x + binormals.points[i].x);
		double Y2_rail2g = Yg2 + distance_factor*(10*normals.points[i].y + binormals.points[i].y);
		double Z2_rail2g = Zg2 + distance_factor*(10*normals.points[i].z + binormals.points[i].z);

		double X3_rail2g = Xg2 + distance_factor*(10*normals.points[i].x - binormals.points[i].x);
		double Y3_rail2g = Yg2 + distance_factor*(10*normals.points[i].y - binormals.points[i].y);
		double Z3_rail2g = Zg2 + distance_factor*(10*normals.points[i].z - binormals.points[i].z);

		double X4_rail2g = Xg2 - distance_factor*(-8*normals.points[i].x + binormals.points[i].x);
		double Y4_rail2g = Yg2 - distance_factor*(-8*normals.points[i].y + binormals.points[i].y);
		double Z4_rail2g = Zg2 - distance_factor*(-8*normals.points[i].z + binormals.points[i].z);

		double X5_rail2g = X2g2 - distance_factor*(-8*normals.points[i].x - binormals.points[i].x);
		double Y5_rail2g = Y2g2 - distance_factor*(-8*normals.points[i].y - binormals.points[i].y);
		double Z5_rail2g = Z2g2 - distance_factor*(-8*normals.points[i].z - binormals.points[i].z);

		double X6_rail2g = X2g2 + distance_factor*(10*normals.points[i].x + binormals.points[i].x);
		double Y6_rail2g = Y2g2 + distance_factor*(10*normals.points[i].y + binormals.points[i].y);
		double Z6_rail2g = Z2g2 + distance_factor*(10*normals.points[i].z + binormals.points[i].z);

		double X7_rail2g = X2g2 + distance_factor*(10*normals.points[i].x - binormals.points[i].x);
		double Y7_rail2g = Y2g2 + distance_factor*(10*normals.points[i].y - binormals.points[i].y);
		double Z7_rail2g = Z2g2 + distance_factor*(10*normals.points[i].z - binormals.points[i].z);

		double X8_rail2g = X2g2 - distance_factor*(-8*normals.points[i].x + binormals.points[i].x);
		double Y8_rail2g = Y2g2 - distance_factor*(-8*normals.points[i].y + binormals.points[i].y);
		double Z8_rail2g = Z2g2 - distance_factor*(-8*normals.points[i].z + binormals.points[i].z);

		glNormal3f(-normals.points[i].x, -normals.points[i].y, -normals.points[i].z);
		// Right face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail2g, Y1_rail2g, Z1_rail2g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X2_rail2g, Y2_rail2g, Z2_rail2g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X6_rail2g, Y6_rail2g, Z6_rail2g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail2g, Y5_rail2g, Z5_rail2g);

		glNormal3f(-binormals.points[i].x, -binormals.points[i].y, -binormals.points[i].z);
		// Top face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail2g, Y1_rail2g, Z2_rail2g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X3_rail2g, Y3_rail2g, Z3_rail2g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X7_rail2g, Y7_rail2g, Z7_rail2g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X6_rail2g, Y6_rail2g, Z6_rail2g);

		glNormal3f(normals.points[i].x, normals.points[i].y, normals.points[i].z);
		// Left face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X3_rail2g, Y3_rail2g, Z3_rail2g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail2g, Y4_rail2g, Z4_rail2g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail2g, Y8_rail2g, Z8_rail2g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X7_rail2g, Y7_rail2g, Z7_rail2g);

		glNormal3f(binormals.points[i].x, binormals.points[i].y, binormals.points[i].z);
		// Bottom face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail2g, Y1_rail2g, Z1_rail2g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail2g, Y4_rail2g, Z4_rail2g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail2g, Y8_rail2g, Z8_rail2g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail2g, Y5_rail2g, Z5_rail2g);

		}
	}

	glEnd();
	glDisable(GL_TEXTURE_2D);

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glBindTexture(GL_TEXTURE_2D, texture_wood);
	// turn on texture mapping
	glEnable(GL_TEXTURE_2D);


	glBegin(GL_QUADS);

	for (int i = 0; i < normals.numControlpoints; i++) {

		double Xg = level1Spline.points[i].x;
		double Yg = level1Spline.points[i].y;
		double Zg = level1Spline.points[i].z;
		
		double X2g = level1Spline.points[i].x + normals.points[i].x / 2.75;
		double Y2g = level1Spline.points[i].y + normals.points[i].y / 2.75;
		double Z2g = level1Spline.points[i].z + normals.points[i].z / 2.75;

		double X1_rail1g = Xg - distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y1_rail1g = Yg - distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z1_rail1g = Zg - distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X2_rail1g = Xg + distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y2_rail1g = Yg + distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z2_rail1g = Zg + distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X3_rail1g = Xg + distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y3_rail1g = Yg + distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z3_rail1g = Zg + distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X4_rail1g = Xg - distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y4_rail1g = Yg - distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z4_rail1g = Zg - distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X5_rail1g = X2g - distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y5_rail1g = Y2g - distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z5_rail1g = Z2g - distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X6_rail1g = X2g + distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y6_rail1g = Y2g + distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z6_rail1g = Z2g + distance_factor*(normals.points[i].z + binormals.points[i].z);

		double X7_rail1g = X2g + distance_factor*(normals.points[i].x - binormals.points[i].x);
		double Y7_rail1g = Y2g + distance_factor*(normals.points[i].y - binormals.points[i].y);
		double Z7_rail1g = Z2g + distance_factor*(normals.points[i].z - binormals.points[i].z);

		double X8_rail1g = X2g - distance_factor*(normals.points[i].x + binormals.points[i].x);
		double Y8_rail1g = Y2g - distance_factor*(normals.points[i].y + binormals.points[i].y);
		double Z8_rail1g = Z2g - distance_factor*(normals.points[i].z + binormals.points[i].z);


		//Rail1g

		glNormal3f(-normals.points[i].x, -normals.points[i].y, -normals.points[i].z);
		// Right face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail1g, Y1_rail1g, Z1_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X2_rail1g, Y2_rail1g, Z2_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X6_rail1g, Y6_rail1g, Z6_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail1g, Y5_rail1g, Z5_rail1g);

		glNormal3f(-binormals.points[i].x, -binormals.points[i].y, -binormals.points[i].z);
		// Top face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X2_rail1g, Y2_rail1g, Z2_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X3_rail1g, Y3_rail1g, Z3_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X7_rail1g, Y7_rail1g, Z7_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X6_rail1g, Y6_rail1g, Z6_rail1g);

		glNormal3f(normals.points[i].x, normals.points[i].y, normals.points[i].z);
		// Left face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X3_rail1g, Y3_rail1g, Z3_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail1g, Y4_rail1g, Z4_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail1g, Y8_rail1g, Z8_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X7_rail1g, Y7_rail1g, Z7_rail1g);

		glNormal3f(binormals.points[i].x, binormals.points[i].y, binormals.points[i].z);
		// Bottom face
		glTexCoord2f(1.0, 0.0);
		glVertex3f(X1_rail1g, Y1_rail1g, Z1_rail1g);
		glTexCoord2f(0.0, 0.0);
		glVertex3f(X4_rail1g, Y4_rail1g, Z4_rail1g);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(X8_rail1g, Y8_rail1g, Z8_rail1g);
		glTexCoord2f(1.0, 1.0);
		glVertex3f(X5_rail1g, Y5_rail1g, Z5_rail1g);

	}

	glEnd();
	glDisable(GL_TEXTURE_2D);
	glEndList();

	/* setup gl view here */
	glClearColor(0.0, 0.0, 0.0, 0.0);
}

int main (int argc, char ** argv)
{
  if (argc<2)
  {  
  printf ("usage: %s <trackfile>\n", argv[0]);
  exit(0);
  }

  loadSplines(argv[1]);

  	/* initialization */
  	glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	/* create window */
	glutInitWindowSize(640, 480);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Assignment 2 - Rollercoaster-AakashShanbhag");

	/* used for double buffering */
	glEnable(GL_DEPTH_TEST);

	/* tells glut to use a particular display function to redraw */
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);

	/* replace with any animate code */
	glutIdleFunc(doIdle);

	/* callback for mouse drags */
  	glutMotionFunc(mousedrag);
  
 	/* callback for idle mouse movement */
  	glutPassiveMotionFunc(mouseidle);
  
  	/* callback for mouse button changes */
  	glutMouseFunc(mousebutton);
 
	/* callback for keyboard */
  	glutKeyboardFunc(keyboardFunc);

	/* enable materials */
	myInit();

	glutMainLoop();

  return 0;
}
