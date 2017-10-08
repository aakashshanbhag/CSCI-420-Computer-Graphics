// the headers

#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>


int g_vMousePos[2] = {0, 0};
int g_iLeftMouseButton = 0;    /* 1 if pressed, 0 if not */
int g_iMiddleMouseButton = 0;
int g_iRightMouseButton = 0;

typedef enum { ROTATE, TRANSLATE, SCALE } CONTROLSTATE;

CONTROLSTATE g_ControlState = ROTATE;

float g_vLandRotate[3] = {0.0, 0.0, 0.0};
float g_vLandTranslate[3] = {0.0, 0.0, 0.0};
float g_vLandScale[3] = {1.0, 1.0, 1.0};

// called before main loop
void init() 
{
    glClearColor(0.0, 0.0, 0.0, 0.0);   // set background color
    glEnable(GL_DEPTH_TEST);            // enable depth buffering
    glShadeModel(GL_FLAT);            // interpolate colors during rasterization
}

// display a frame
void display()
{
    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity(); // reset transformation

    // draw a triangle
    glBegin(GL_TRIANGLES);
	//glColor4f(1.0,0.0,0.0,1.0);
        glVertex3f(0.0, 0.0, -10.0);
	//glColor4f(0.0,1.0,0.0,1.0);
        glVertex3f(-5.0, 0.0, -10.0);
	//glColor4f(0.0,0.0,1.0,1.0);
        glVertex3f(0.0, 1.0, -10.0);
    glEnd();

    glutSwapBuffers(); // double buffer flush
}

// called every time window is resized to update projection matrix
void reshape(int w, int h)
{
    // setup image size
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // setup camera
    glFrustum(-0.1, 0.10, -float(h)/(10.0*float(w)), float(h)/(10.0*float(w)), 0.5, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

}

void mousedrag(int x, int y)
{
  int vMouseDelta[2] = {x-g_vMousePos[0], y-g_vMousePos[1]};

  
  switch (g_ControlState)
  {
    case TRANSLATE:  
      if (g_iLeftMouseButton)
      {
        g_vLandTranslate[0] += vMouseDelta[0]*0.01;
        g_vLandTranslate[1] -= vMouseDelta[1]*0.01;
      }
      if (g_iMiddleMouseButton)
      {
        g_vLandTranslate[2] += vMouseDelta[1]*0.01;
      }
      glMatrixMode(GL_MODELVIEW);
      glTranslatef(g_vLandTranslate[0],g_vLandTranslate[1],g_vLandTranslate[2]);
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
//glMatrixMode(GL_MODELVIEW);
  //    glRotatef(10.0,g_vLandRotate[0],g_vLandRotate[1],g_vLandRotate[2]);
     // cout<<"Mouse drag "<<endl;
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
      glMatrixMode(GL_MODELVIEW);
      glScalef(g_vLandScale[0],g_vLandScale[1],g_vLandScale[2]);
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
    case GLUT_ACTIVE_CTRL:
      g_ControlState = TRANSLATE;
      break;
    case GLUT_ACTIVE_SHIFT:
      g_ControlState = SCALE;
      break;
    default:
      g_ControlState = ROTATE;
      break;
  }

  g_vMousePos[0] = x;
  g_vMousePos[1] = y;
}





// entry point
int main(int argc, char **argv)
{
    
    // initialize GLUT
    glutInit(&argc, argv);
    
    // request double buffer
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    
    // set window size
    glutInitWindowSize(500, 500);
    
    // set window position
    glutInitWindowPosition(0, 0);
    
    // creates a window
    glutCreateWindow("Ahahaha!");


    glutMotionFunc(mousedrag);
  /* callback for idle mouse movement */
  glutPassiveMotionFunc(mouseidle);
  /* callback for mouse button changes */
  glutMouseFunc(mousebutton);

    // initialize states
    init();

    // GLUT callbacks
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);


    // start GLUT program
    glutMainLoop();
    return 0;
}
