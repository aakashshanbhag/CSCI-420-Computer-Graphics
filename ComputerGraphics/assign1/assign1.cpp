/*
  CSCI 420 Computer Graphics
  Assignment 1: Height Fields
  Author: Aakash Shanbhag
  USC ID:3205699915
*/

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>

using namespace std;

//Global variable description

int g_iMenuId;

// Setting up the screenshot counter and the total number of frames.
const int totalframes = 301;
int screenshotCounter = 0;
bool Screenshots_flag = false;

//OpenGLMatrix matrices;
const float FOV = 60.0; 
GLfloat rotate_factor=2.0;

// Defining the window resolution.
float screenwidth=640.0;
float screenheight=480.0;

// Initilization of the mouse and the keyboard parameters.
int g_vMousePos[2] = {0, 0};
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;
bool stop=0;

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROLSTATE;

CONTROLSTATE g_ControlState = ROTATE;

/* state of the world */
float g_vLandRotate[3] = {0.0, 0.0, 0.0};
float g_vLandTranslate[3] = {0.0, 0.0, 0.0};
float g_vLandScale[3] = {1.0, 1.0, 1.0};


/* see <your pic directory>/pic.h for type Pic */
Pic * g_pHeightData=NULL;
Pic * image2_c=NULL;

// Function to create the heightfield from the img image and the color form the img1 image.
void createHeightField(Pic * img, Pic *img1)
{

  int img_height= img->ny;
  int img_width= img->nx;
  float Scale =0.125* img_height/ 100.0;// scale corresponding to the max height. 

  // Generate our vertices--we should later use a triangle strip, that will make this so much more efficient. We center the heightmap at the origin.
  if(img1==NULL || (img_width!=img1->nx) || (img_height!= img->ny))
  {
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
  } // If img1 exists. Use its intensities as the the color map.
  else{
    for (int i = 0; i < img_height - 1; i++) {
     glBegin(GL_TRIANGLE_STRIP);
        for (int j = 0; j < img_width - 1; j++) {
 
        float intensity_lower_left=PIC_PIXEL(img,j, i, 0);
        float intensity_upper_left=PIC_PIXEL(img,j, i + 1 , 0);

        float intensity_lower_left_c=PIC_PIXEL(img1,j, i, 0)/255.0;
        float intensity_upper_left_c=PIC_PIXEL(img1,j, i + 1 , 0)/255.0;

        glColor4f(intensity_lower_left_c,intensity_lower_left_c,intensity_lower_left_c,1.0);
        glVertex3f(j,i,GLfloat(Scale * intensity_lower_left));
        
        glColor4f(intensity_upper_left_c,intensity_upper_left_c,intensity_upper_left_c,1.0);
        glVertex3f(j,i+1, GLfloat(Scale * intensity_upper_left));

        }
     glEnd();
    }

    }

    
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

void myinit(int argc, char *argv[])
{
  // Tracking if the image is loaded sucessfully
 if (g_pHeightData!=NULL)
  {
    cout<<"Read image Successfully"<<endl;
  }
  else
  {
    cout<<"Failed to read image"<<endl;
  }

  /* setup gl view here */
  glClear(GL_COLOR_BUFFER_BIT |GL_DEPTH_BUFFER_BIT);
 
  // enable depth buffering
  glEnable(GL_DEPTH_TEST);
  
  // interpolate colors during rasterization
  glShadeModel(GL_SMOOTH); 


}
// Function to obtain the heights and the widths of the images to be loaded
float get_height(Pic *image)
{
  float height=image->ny;
  return height;
}

float get_width(Pic * image)
{
  float width=image->nx;
  return width;
}

void display()
{
    // reset transformation
    glLoadIdentity(); 
    //clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
   // glTranslatef(-1.0*get_width(g_pHeightData)/2.0,-1.0*get_height(g_pHeightData)/2.0,-500.0);
    glTranslatef(-1.0*get_width(g_pHeightData)/2.0,-1.0*get_height(g_pHeightData)/2.0,-750.0);
    
    // keeping the transformations local by using the push and pop schemes
    glPushMatrix();
    // Translating on mouse commands
    glTranslatef(1.0*g_vLandTranslate[0],1.0*g_vLandTranslate[1],1.0*g_vLandTranslate[2]);
    // Scaling on mouse commands.
    glScalef(g_vLandScale[0],g_vLandScale[1],g_vLandScale[2]);
    // Rotate using mouse callback
    glRotatef(g_vLandRotate[0], 1.0, 0.0,0.0);
    glRotatef(g_vLandRotate[1], 0.0, 1.0,0.0);
    glRotatef(g_vLandRotate[2], 0.0, 0.0,1.0);
    // Generating heightfield per frame.
    createHeightField(g_pHeightData,image2_c);
    glPopMatrix();

    glutSwapBuffers(); // double buffer flush

}

void menufunc(int value)
{
  switch (value)
  {
    case 0:
      exit(0);
      break;
  }
}

void doIdle()
{
  // Screen shot files to be created for the animation.
  if (screenshotCounter < totalframes && Screenshots_flag) {
    std::stringstream ss;
    ss << "../assign1/TrialAnimation/trial-" << screenshotCounter << ".jpg";
    saveScreenshot ((char*)ss.str().data());
    screenshotCounter++;
  }
  /* make the screen update */
  glutPostRedisplay();
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

void mousebutton(int button, int state, int x, int y)
{

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
    //gluLookAt(0.0,0.0,500.0,0.0,0.0,0.0,0.0,1.0,0.0);
    glMatrixMode(GL_MODELVIEW);
  
}

void keyboardFunc(unsigned char key, int x, int y) {
  switch (key) {
    case 27: // ESC key
      exit(0); // exit the program
      break;

    case 'x':
      // take a screenshot
      saveScreenshot("screenshot.jpg");
      break;

    case 'a':
      // Render points
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      cout<<"Points rendered"<<endl;
      break;

    case 'A':
      // Render points
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      cout<<"Points rendered"<<endl;
      break;

    case 's':
      // Render Wireframe
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      cout<<"Wireframe rendered"<<endl;
      break;

    case 'S':
      // Render Wireframe
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      cout<<"Wireframe rendered"<<endl;
      break;

    case 'd':
      // Render Solid trianges
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      cout<<"Solid trianges rendered"<<endl;
      break;

    case 'D':
      // Render Solid trianges
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      cout<<"Solid trianges rendered"<<endl;
      break;

    case 'q':
      // Screenshot capture to commence
      Screenshots_flag = true;
      std::cout << "Starting animation!." << std::endl;
      break;

    case 'y':
      // Issue with CTRL and ALT keys dealt with y key.
      g_ControlState = SCALE;
    break;

  }
}

int main (int argc, char ** argv)
{
  if (argc<3)
  {  
    printf ("usage: %s heightfield.jpg colormap.jpg \n", argv[0]);
    exit(1);
  }

  g_pHeightData = jpeg_read(argv[1], NULL);
  image2_c=jpeg_read(argv[2],NULL);

  if((g_pHeightData->ny!= image2_c->ny) || (g_pHeightData->nx!= image2_c->nx))
  {cout<<"Sizes of the images do not match and hence color map from second image cannot be generated"<<endl;}
  

  if (!g_pHeightData)
  {
    printf ("error reading %s.\n", argv[1]);
    exit(1);
  }

  glutInit(&argc,argv);
  
  // request double buffer
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

  // set window size
  glutInitWindowSize(screenwidth , screenheight);
    
  // set window position
  glutInitWindowPosition(0, 0);
    
  // creates a window
  glutCreateWindow("HW-1-Aakash_Shanbhag");

  cout << "OpenGL Version: " << glGetString(GL_VERSION) << endl;
  cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
  cout << "Shading Language Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

  
  glutDisplayFunc(display);
  
  /* allow the user to quit using the right mouse button menu */
  g_iMenuId = glutCreateMenu(menufunc);
  glutSetMenu(g_iMenuId);
  glutAddMenuEntry("Quit",0);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  
  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mousedrag);
  
  /* callback for idle mouse movement */
  glutPassiveMotionFunc(mouseidle);
  
  /* callback for mouse button changes */
  glutMouseFunc(mousebutton);
 
  /* tells glut to use a particular display function to redraw */
  glutReshapeFunc(reshape);

  // callback for keyboard
  glutKeyboardFunc(keyboardFunc);
  
  /* do initialization */
  myinit(argc, argv);

  glutMainLoop();
  
  return(0);
}
