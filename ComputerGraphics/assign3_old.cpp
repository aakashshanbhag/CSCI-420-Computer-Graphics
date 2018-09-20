/*
CSCI 420: Computer Graphics.
Assignment 3: Raytracer.
Name: Aakash Shanbhag.
USC ID:3205699915.
*/

//ISSUES ON THE CONVERSION BETWEEN DOUBLE AND UNSIGNED CHAR>>>>> SOFT SHADOWS MISSING 


#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include <iostream>
#include <math.h>
using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define PI 3.14159265
const double MAX_DIST = -1e8;
const double BIAS = 1e-16;

unsigned char buffer[HEIGHT][WIDTH][3];


struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex vertex[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

typedef struct _Raycaster
{
  // By default origin at 0,0,0
  double origin[3]={0,0,0};
  double direction[3];
} Raycaster;

typedef struct Intensity
{
  // Default max intensity (White color background)
  double r=1;
  double g=1;
  double b=1;  
};

// Creating a global intensity map for the frame created
Intensity intensity_map[WIDTH][HEIGHT];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,double r,double g,double b);
void plot_pixel_jpeg(int x,int y,double r,double g,double b);
void plot_pixel(int x,int y,double r,double g,double b);
double dot_product(double a[3],double b[3]);
double area_triangle(double a[3],double b[3],double c[3]);
void raycast();
Intensity ray_trace(Raycaster ray);
Raycaster shadowRay_create(double m[3],Light l);
bool sphere_light_intersect(Sphere s, Raycaster ray,double &t_min, bool &inside);
bool sphere_shadow_intersect(Sphere s, Raycaster ray);
bool triangle_light_intersect(Triangle t, Raycaster ray,double &t_param,double &alpha,double &beta,double &gamma);
bool triangle_shadow_intersect(Triangle t, Raycaster ray);
Intensity sphere_phong(const Sphere s,const Raycaster shadow,const Light l,const double m[3],const bool inside);
Intensity triangle_phong(const Triangle t, const Raycaster shadow, const Light L,const double diffuse_color[3],const double specular_color[3], const double shininess, const double m[3], double n[3]);


double dot_product(double a[3], double b[3]) 
{
  double output = 0;
  for (int i = 0; i < 3; i++)
  {
    output = output + a[i] * b[i];
  }
  return output;
}
// Cross product yields the area of the triangle from 3 vertices.
double area_triangle(double a[2], double b[2], double c[2]) 
{
  double output=0;
  // AB vector
  double ab_x = b[0] - a[0];
  double ab_y = b[1] - a[1];
  // AC vector
  double ac_x = c[0] - a[0];
  double ac_y = c[1] - a[1];
  // Area= 0.5* cross(AB,AC)
  output=0.5*((ab_x * ac_y) - (ac_x * ab_y)) ;
  return output;
}

void draw_scene()
{
  raycast();
  unsigned int x,y;
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      plot_pixel(x,y,intensity_map[x][y].r,intensity_map[x][y].g,intensity_map[x][y].b); 
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, double r, double g, double b)
{
  glColor3f(r,g,b);
  glVertex2f(x, y);
}

void plot_pixel_jpeg(int x,int y,double r,double g,double b)
{
  unsigned char R = int(r * 255);
  unsigned char G = int(g * 255);
  unsigned char B= int(b * 255);

  buffer[HEIGHT-y-1][x][0]=R;
  buffer[HEIGHT-y-1][x][1]=G;
  buffer[HEIGHT-y-1][x][2]=B;
}

void plot_pixel(int x,int y,double r,double g, double b)
{
 plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      
}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.vertex[j].position);
	      parse_doubles(file,"nor:",t.vertex[j].normal);
	      parse_doubles(file,"dif:",t.vertex[j].color_diffuse);
	      parse_doubles(file,"spe:",t.vertex[j].color_specular);
	      parse_shi(file,&t.vertex[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}
// Create all new functions for estimating the creation of all the rays
// and the intersections respectively.

//--------------------------------------------------------------------------------------

// Cast rays and map according to the focal length and field view // Scale by required factor
void raycast()
{
  Raycaster ray;
  //cout<<"Camera origin is "<<ray.origin[0]<<"\t"<<ray.origin[1]<<"\t"<<ray.origin[2]<<endl;

  // Create drawing canvas.
  double aspectRatio=(double)WIDTH / (double)HEIGHT;
  double x_right= aspectRatio * tan(fov* PI / 360);// Since angle is in radians 
  double x_left= -1*x_right;
  double y_top= tan(fov * PI / 360);
  double y_bottom= -1*y_top;
  double z=-1;// focal length fixed at 1 in the neagtive z direction.
  
  double frame_width= x_right-x_left;
  double frame_height= y_top- y_bottom;

  double scale_width = frame_width/WIDTH;
  //cout << "image scale_width : " << scale_width<< endl;
  
  double scale_height = frame_height/HEIGHT;
 // cout << "image scale_height : " << scale_height << endl;

  for (int i = 0; i < HEIGHT; i++){
    for (int j = 0; j < WIDTH; j++)
    {
      double x = x_left + (1.0 * j + 0.5) * scale_width;
      double y = y_bottom + (1.0 * i + 0.5) * scale_height;



      double magnitude = sqrt(x*x + y*y + 1); // Check z*z= 1
      //cout<<"Mag: "<<magnitude<<endl;

      ray.direction[0] = x / magnitude;
      ray.direction[1] = y / magnitude;
      ray.direction[2] = -1 / magnitude;// z=-1 always.

      //cout<<ray.direction[0]<<"\t"<<ray.direction[1]<<"\t"<<ray.direction[2]<<endl;

      double direction_magnitude= sqrt((ray.direction[0]*ray.direction[0])+(ray.direction[1]*ray.direction[1])+(ray.direction[2]*ray.direction[2]));
      // cout<<"Dir mag: "<<direction_magnitude<<endl;
      
      
         intensity_map[j][i] = ray_trace(ray);

    
    }
  }
}

Raycaster shadowRay(double m[3],Light l)
{
  Raycaster shadowray;

  // Shadow object intersection origin matters in this case.
  shadowray.origin[0]=m[0];
  shadowray.origin[1]=m[1];
  shadowray.origin[2]=m[2];
  
  // Direction in terms of the light and the intersection.
  shadowray.direction[0]=l.position[0]-m[0];
  shadowray.direction[1]=l.position[1]-m[1];
  shadowray.direction[2]=l.position[2]-m[2];

  double magnitude= sqrt(shadowray.direction[0]*shadowray.direction[0]+shadowray.direction[1]*shadowray.direction[1]+shadowray.direction[2]*shadowray.direction[2]);

  // Normalising and checking if the unity magnitude is maintained in the entire process.

  shadowray.direction[0]/=magnitude;
  shadowray.direction[1]/=magnitude;
  shadowray.direction[2]/=magnitude;


  return shadowray;
}

Intensity ray_trace(Raycaster ray)
{

 Intensity intensity;
 //cout <<"Default intensity is "<<intensity.r<<"\t"<<intensity.g<<"\t"<<intensity.b<<endl; 

  Sphere s;
  double t_sphere = -1;// parameteric t check  The image frame is at z=-1.
  bool inside = false;

  /* test if the ray intersects a sphere */
  for (int i = 0; i < num_spheres; i++)
   {
    Sphere sphere = spheres[i];
    double t_min;
    if (sphere_light_intersect(sphere,ray, t_min, inside)) {
      //cout<<"t_min"<<t_min<<endl;
      if (t_sphere == -1 || t_min < t_sphere) 
      {
        t_sphere = t_min;
        s = sphere;
      }
    }
  }

  Triangle t;
  double alpha,beta,gamma;
  double temp_alpha,temp_beta,temp_gamma;
  double t_triangle=-1;

  /* test if the ray intersects a triangle */
  for (int i = 0; i < num_triangles; i++) {
    Triangle triangle = triangles[i];
    double t_min;
    if (triangle_light_intersect(triangle,ray, t_min, temp_alpha, temp_beta, temp_gamma)) {
      if (t_triangle== -1 || t_min < t_triangle) {
        t_triangle = t_min;
        t = triangle;
        alpha = temp_alpha;
        beta = temp_beta;
        gamma = temp_gamma;
      }
    }
  }


  /* if a ray intersects both a sphere and a triangle, determine which is closer */
  if (t_sphere != -1 && t_triangle != -1)
   {
    if (t_sphere <= t_triangle)
    {
      t_triangle = -1;
    }
    else 
    {
     t_sphere = -1;
    }
  }
// SPHERE INTERSECTION.

  if (t_sphere != -1)
  {

  cout<<"Rendering SPHERE"<<endl;
  intensity.r = ambient_light[0];
  intensity.g = ambient_light[1];
  intensity.b = ambient_light[2];

  // Intesection point in 3D
  double m[3]={0,0,0};

  m[0]= ray.origin[0] + ray.direction[0] * t_sphere;
  m[1]= ray.origin[1] + ray.direction[1] * t_sphere;
  m[2]= ray.origin[2] + ray.direction[2] * t_sphere;

  //cout<<"M values "<<m[0]<<"\t"<<m[1]<<"\t"<<m[2]<<endl; 

  // For each of the lights in the scene we need to create the phong model 
  for (int i = 0; i < num_lights; i++) 
  {
    Light l = lights[i];
    
    Raycaster shadowray ;
      
    shadowray.origin[0]=m[0];
    shadowray.origin[1]=m[1];
    shadowray.origin[2]=m[2];
  
    // Direction in terms of the light and the intersection.
    shadowray.direction[0]=l.position[0]-m[0];
    shadowray.direction[1]=l.position[1]-m[1];
    shadowray.direction[2]=l.position[2]-m[2];

    double magnitude= sqrt(shadowray.direction[0]*shadowray.direction[0]+shadowray.direction[1]*shadowray.direction[1]+shadowray.direction[2]*shadowray.direction[2]);

   // Normalising and checking if the unity magnitude is maintained in the entire process.

    shadowray.direction[0]/=magnitude;
    shadowray.direction[1]/=magnitude;
    shadowray.direction[2]/=magnitude;

      
      bool block=false;

      for (int j = 0; j < num_spheres; j++)
      {
        Sphere sphere = spheres[j];
        /* determine if shadow ray is blocked by any objects */
        if (sphere_shadow_intersect(sphere,shadowray)) 
        {
          block = true;
          cout<<"SHADOW-S"<<endl;
          break;
        }
      }

    if(!block)
      {
        for (int j = 0; j < num_triangles; j++) {
          Triangle triangle = triangles[j];
          /* determine if shadow ray is blocked by any objects */
          if (triangle_shadow_intersect(triangle,shadowray)) {
            block = true;
             cout<<"SHADOW-T"<<endl;
            break;
          }
        }
      }

      if (!block)
      {
        Intensity sphere_color= sphere_phong( s,shadowray, l,  m, inside); 
        intensity.r += sphere_color.r;
        intensity.g += sphere_color.g;
        intensity.b += sphere_color.b;
      }
    }
  }

  else if(t_triangle != -1)
  {

  cout<<"Rendering triangle"<<endl;
  intensity.r = ambient_light[0];
  intensity.g = ambient_light[1];
  intensity.b = ambient_light[2];

  double m[3]={0,0,0};

  m[0]= ray.origin[0] + ray.direction[0] * t_triangle;
  m[1]= ray.origin[1] + ray.direction[1] * t_triangle;
  m[2]= ray.origin[2] + ray.direction[2] * t_triangle;

 // cout<<"M values "<<m[0]<<"\t"<<m[1]<<"\t"<<m[2]<<endl;

  // PHONG Model for getting the normalised values for each component on the basis of the normals

  /* calculate interpolated diffuse color */
    double diffuse_color[3];
    for (int i = 0; i < 3; i++) {
      diffuse_color[i] = ((alpha * t.vertex[0].color_diffuse[i])
        + (beta * t.vertex[1].color_diffuse[i])
        + (gamma * t.vertex[2].color_diffuse[i]));
    }

    /* calculate interpolated normal of intersection point on triangle */
    double normal[3];
    for (int i = 0; i < 3; i++) {
      normal[i] = ((alpha * t.vertex[0].normal[i])
        + (beta * t.vertex[1].normal[i])
        + (gamma * t.vertex[2].normal[i]));
    }

    /* calculate interpolated specular color */
    double specular_color[3];
    for (int i = 0; i < 3; i++) {
      specular_color[i] = ((alpha * t.vertex[0].color_specular[i])
        + (beta * t.vertex[1].color_specular[i])
        + (gamma * t.vertex[2].color_specular[i]));
    }

    /* calculate interpolated shininess */
    double shininess = (alpha * t.vertex[0].shininess) + (beta  * t.vertex[1].shininess) + (gamma * t.vertex[2].shininess);

    for (int i = 0; i < num_lights; i++) {
      Light l = lights[i];

    Raycaster shadowray ;
      
    shadowray.origin[0]=m[0];
    shadowray.origin[1]=m[1];
    shadowray.origin[2]=m[2];
  
    // Direction in terms of the light and the intersection.
    shadowray.direction[0]=l.position[0]-m[0];
    shadowray.direction[1]=l.position[1]-m[1];
    shadowray.direction[2]=l.position[2]-m[2];

    double magnitude= sqrt(shadowray.direction[0]*shadowray.direction[0]+shadowray.direction[1]*shadowray.direction[1]+shadowray.direction[2]*shadowray.direction[2]);

   // Normalising and checking if the unity magnitude is maintained in the entire process.

    shadowray.direction[0]/=magnitude;
    shadowray.direction[1]/=magnitude;
    shadowray.direction[2]/=magnitude;

      
      bool block=false;

      for (int j = 0; j < num_spheres; j++)
      {
        Sphere sphere = spheres[j];
        /* determine if shadow ray is blocked by any objects */
        if (sphere_shadow_intersect(sphere,shadowray)) 
        {
          block = true;
          cout<<"SHADOW-S"<<endl;
          break;
        }
      }

      if(!block)
      {
        for (int j = 0; j < num_triangles; j++) {
          Triangle triangle = triangles[j];
          /* determine if shadow ray is blocked by any objects */
          if (triangle_shadow_intersect(triangle,shadowray) ) {
            block = true;
             cout<<"SHADOW-T"<<endl;
            break;
          }
        }

      }

      /* if no objects are blocking the shadow ray, calculate color */

      if (!block) {
        Intensity triangle_color = triangle_phong(
          t, shadowray, l, diffuse_color, specular_color, shininess, m, normal);
        
        intensity.r += triangle_color.r;
        intensity.g += triangle_color.g;
        intensity.b += triangle_color.b;
      }

  }
  }

return intensity;

}


bool sphere_light_intersect(Sphere s, Raycaster ray, double &t_min, bool &inside)
{
  // According to the formulae stated in Lecture 8.1 for sphere intersections.

    double x = ray.origin[0] - s.position[0];
    double y = ray.origin[1] - s.position[1];
    double z = ray.origin[2] - s.position[2];

    double a = 1 ;// Should be equal to 1 as it is magnitude of the direction vector. 
    double b = 2 * ((ray.direction[0] * x) + (ray.direction[1] * y) + (ray.direction[2] * z));
    double c = (x * x) + (y * y) + (z * z) - (s.radius * s.radius);

    double discriminant = (b * b) - (4 * a * c);
    if (discriminant < 0) return false;  // No real intersection

    // The intersection based on the discriminant value
    double t0 = (-b - sqrt(discriminant)) / (2 * a); // point along ray which enters/exits the sphere
    double t1 = (-b + sqrt(discriminant)) / (2 * a); // point along ray which exits/enters the sphere

    if (t0 < 0 && t1 < 0) {
      return false;
    }
    // Closer point if both values are positive and eliminate negative value.
    if (t0 < 0) {
      t_min = t1;
      inside = true;
    }
    else if (t1 < 0) {
      t_min = t0;
      inside = true;
    }
    else {
      t_min = min(t0, t1);
    }

    return true;
}

bool sphere_shadow_intersect(Sphere s,Raycaster ray)
{
    double x = ray.origin[0] - s.position[0];
    double y = ray.origin[1] - s.position[1];
    double z = ray.origin[2] - s.position[2];

    double a = 1; // magnitude of unit direction vector.
    double b = 2 * ((ray.direction[0] * x) + (ray.direction[1] * y) + (ray.direction[2] * z));
    double c = (x * x) + (y * y) + (z * z) - (s.radius * s.radius);

    double discriminant = (b * b) - (4 * a * c);
    if (discriminant < 0) {
      return false; // no intersection
    }
    
    double t0 = (-b - sqrt(discriminant)) / 2;
    double t1 = (-b + sqrt(discriminant)) / 2; 

    if (t0 < 0 || t1 < 0) {
      return false;
    }
    else {
      return true;
    }
}

// m is the point where the ray intesects the sphere surface.
Intensity sphere_phong(const Sphere s,const Raycaster shadow,const Light L,const double m[3],const bool inside)
{
  double n[3]={0,0,0};

  // Calculation of the normal at the point of intersection m.

  n[0] = (m[0] - s.position[0]) / s.radius;
  n[1] = (m[1] - s.position[1]) / s.radius;
  n[2] = (m[2] - s.position[2]) / s.radius;

  /* Change the orientation if the ray is inside the sphere */
  if (inside) {
    n[0] *= -1;
    n[1] *= -1;
    n[2] *= -1;
  }

  // For the shadow ray, we need to create a local copy for then calculation of the reflected ray.

  Raycaster l = shadow;// could be blocked when in shadow or direct when unblocked 
  double l_dot_n=dot_product(l.direction,n);
  
  // Clamp between 0 and 1
  
  if(l_dot_n < 0)    {l_dot_n=0;}
  else if(l_dot_n>1) {l_dot_n=1;}
  
  //cout<<"LdotN: "<<l_dot_n<<endl;

  // Calculate the reflected ray r for the specular component = 2*(l_dot_n)n-l.

  Raycaster r= l;
  r.direction[0]= 2 * l_dot_n * n[0]-l.direction[0];
  r.direction[1]= 2 * l_dot_n * n[1]-l.direction[1];
  r.direction[2]= 2 * l_dot_n * n[2]-l.direction[2];

  double r_magnitude=sqrt(r.direction[0]*r.direction[0]+r.direction[1]*r.direction[1]+r.direction[2]*r.direction[2]);

  r.direction[0]/=r_magnitude;
  r.direction[1]/=r_magnitude;
  r.direction[2]/=r_magnitude;


  // We need to normalise v for r.v factor of specular component
  double v[3]={0,0,0};
  double magnitude= sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
  v[0]=-1*m[0]/magnitude;
  v[1]=-1*m[1]/magnitude;
  v[2]=-1*m[2]/magnitude;

  double r_dot_v=dot_product(v,r.direction);
  //cout<<"RdotV ORG: "<<r_dot_v<<endl;
  if(r_dot_v < 0)    {r_dot_v=0;}
  else if(r_dot_v>1) {r_dot_v=1;}

  //cout<<"RdotV: "<<r_dot_v<<endl;
  //cout<<endl;

  // Finally applying the Phong model of illumination for the expected output.
  double R = L.color[0] * ((s.color_diffuse[0] * l_dot_n) + s.color_specular[0] * pow(r_dot_v, s.shininess));

  double G = L.color[1] * ((s.color_diffuse[1] * l_dot_n) + s.color_specular[1] * pow(r_dot_v, s.shininess));

  double B = L.color[2] * ((s.color_diffuse[2] * l_dot_n) + s.color_specular[2] * pow(r_dot_v, s.shininess));
  //cout<<"Light color:"<<L.color[0]<<L.color[1]<<L.color[2]<<endl;
  cout<<endl;
  // Clamp up the channels between 0 and 1 
  if(R < 0)       {R=0;}
  else if(R>1)  {R=1;}

  if(G < 0)      {G=0;}
  else if(G>1) {G=1;}

  if(B< 0)       {B=0;}
  else if(B>1) {B=1;}

  //cout<<"R: "<<R<<" G: "<<G<<" B: "<<B<<endl;
  Intensity finalcolor;
  finalcolor.r=R;
  finalcolor.g=G;
  finalcolor.b=B;

  return finalcolor;
}

bool triangle_light_intersect(Triangle t, Raycaster ray,double &t_param,double &alpha,double &beta,double &gamma)
{
  // Define the edges of the triangles
  double edge1[3]={0,0,0};
  double edge2[3]={0,0,0};

  // Calculation of the edges

  for (int i=0;i<3;i++)
  {
    edge1[i]=t.vertex[1].position[i] - t.vertex[0].position[i];
    edge2[i] =t.vertex[2].position[i] - t.vertex[0].position[i];

  }

  // Normal to the faces of the triangles

  double face[3]={0,0,0};

  // Calculation of the cross product

  face[0]= (edge1[1]*edge2[2])-(edge1[2]*edge2[1]);
  face[1]= (edge1[2]*edge2[0])-(edge1[0]*edge2[2]);
  face[2]= (edge1[0]*edge2[1])-(edge1[1]*edge2[0]);

  // Unit vector of the normal to the face.

  double face_mag= sqrt(face[0]*face[0]+face[1]*face[1]+face[2]*face[2]);

  face[0]/=face_mag;
  face[1]/=face_mag;
  face[2]/=face_mag;

  // Intersection of the ray with the face of the triangle

  // Translation of the origin of the ray to the vertex points

  double trans[3]={0,0,0};
  trans[0]= ray.origin[0] - t.vertex[0].position[0];
  trans[1]= ray.origin[1] - t.vertex[1].position[1];
  trans[2]= ray.origin[2] - t.vertex[2].position[2];

  // Checking if the ray direction and the normal direction is orthogonal

  double d_dot_n=dot_product(ray.direction,face);
  if (d_dot_n==0)
    return false;

  t_param= -1.0 * dot_product(trans,face)/d_dot_n;

  // Check if the intersection is behind the camera origin

  if (t_param <0.00001) return false;


  // Calculate 3D point of intersection

  double point_inter[3]={0,0,0};
  point_inter[0]= ray.origin[0] + (t_param * ray.direction[0]);
  point_inter[1]= ray.origin[1] + (t_param * ray.direction[1]);
  point_inter[2]= ray.origin[2] + (t_param * ray.direction[2]);


  // Check if the point is inside the trainagle by calculation of the barycentric coordinates

  double point_estimate[2], v0_projected[2], v1_projected[2], v2_projected[2];

    for (int i = 0; i < 2; i++) {
      point_estimate[i] = point_inter[i];
      v0_projected[i] = t.vertex[0].position[i];
      v1_projected[i] = t.vertex[1].position[i];
      v2_projected[i] = t.vertex[2].position[i];
    }

    double triArea = area_triangle(v0_projected, v1_projected, v2_projected);

    alpha = area_triangle(point_estimate, v1_projected, v2_projected) / triArea;
    beta = area_triangle(v0_projected, point_estimate, v2_projected) / triArea;
    gamma = area_triangle(v0_projected, v1_projected, point_estimate) / triArea;

    if (alpha * beta >= 0 && beta * gamma >= 0) {
     // point is inside the triangle
      return true;
    }
    else {
      return false;
    } 
}

bool triangle_shadow_intersect(Triangle t,Raycaster shadowRay)
{
  // Define the edges of the triangles
  double edge1[3]={0,0,0};
  double edge2[3]={0,0,0};

  // Calculation of the edges

  for (int i=0;i<3;i++)
  {
    edge1[i]=t.vertex[1].position[i] - t.vertex[0].position[i];
    edge2[i] =t.vertex[2].position[i] - t.vertex[0].position[i];

  }

  // Normal to the faces of the triangles

  double face[3]={0,0,0};

  // Calculation of the cross product

  face[0]= (edge1[1]*edge2[2])-(edge1[2]*edge2[1]);
  face[1]= (edge1[2]*edge2[0])-(edge1[0]*edge2[2]);
  face[2]= (edge1[0]*edge2[1])-(edge1[1]*edge2[0]);

  // Unit vector of the normal to the face.

  double face_mag= sqrt(face[0]*face[0]+face[1]*face[1]+face[2]*face[2]);

  face[0]/=face_mag;
  face[1]/=face_mag;
  face[2]/=face_mag;

  // Intersection of the ray with the face of the triangle

  // Translation of the origin of the ray to the vertex points

  double trans[3]={0,0,0};
  trans[0]= shadowRay.origin[0] - t.vertex[0].position[0];
  trans[1]= shadowRay.origin[1] - t.vertex[1].position[1];
  trans[2]= shadowRay.origin[2] - t.vertex[2].position[2];

  // Checking if the ray direction and the normal direction is orthogonal

  double d_dot_n=dot_product(shadowRay.direction,face);
  if (d_dot_n==0)
    return false;

  double t_param= -1.0 * dot_product(trans,face)/d_dot_n;

  // Check if the intersection is behind the camera origin

  if (t_param <0.001) return false;


  // Calculate 3D point of intersection

  double point_inter[3]={0,0,0};
  point_inter[0]= shadowRay.origin[0] + (t_param * shadowRay.direction[0]);
  point_inter[1]= shadowRay.origin[1] + (t_param * shadowRay.direction[1]);
  point_inter[2]= shadowRay.origin[2] + (t_param * shadowRay.direction[2]);


  // Check if the point is inside the trainagle by calculation of the barycentric coordinates

  double point_estimate[2], v0_projected[2], v1_projected[2], v2_projected[2];

    for (int i = 0; i < 2; i++) {
      point_estimate[i] = point_inter[i];
      v0_projected[i] = t.vertex[0].position[i];
      v1_projected[i] = t.vertex[1].position[i];
      v2_projected[i] = t.vertex[2].position[i];
    }

    double triArea = area_triangle(v0_projected, v1_projected, v2_projected);

    double alpha = area_triangle(point_estimate, v1_projected, v2_projected) / triArea;
    double beta = area_triangle(v0_projected, point_estimate, v2_projected) / triArea;
    double gamma = area_triangle(v0_projected, v1_projected, point_estimate) / triArea;

    if (alpha * beta >= 0 && beta * gamma >= 0) {
     // point is inside the triangle
      return true;
    }
    else {
      return false;
    } 
}

Intensity triangle_phong(const Triangle t, const Raycaster shadow, const Light L,const double diffuse_color[3],const double specular_color[3], const double shininess, const double m[3], double n[3])
{
  Raycaster l = shadow;// could be blocked when in shadow or direct when unblocked 
  double l_dot_n=dot_product(l.direction,n);
  
  // Clamp between 0 and 1
  
  if(l_dot_n < 0)    {l_dot_n=0;}
  else if(l_dot_n>1) {l_dot_n=1;}
  
  //cout<<"LdotN: "<<l_dot_n<<endl;

  // Calculate the reflected ray r for the specular component = 2*(l_dot_n)n-l.

  Raycaster r= l;
  r.direction[0]= 2 * l_dot_n * n[0]-l.direction[0];
  r.direction[1]= 2 * l_dot_n * n[1]-l.direction[1];
  r.direction[2]= 2 * l_dot_n * n[2]-l.direction[2];

  double r_magnitude=sqrt(r.direction[0]*r.direction[0]+r.direction[1]*r.direction[1]+r.direction[2]*r.direction[2]);

  r.direction[0]/=r_magnitude;
  r.direction[1]/=r_magnitude;
  r.direction[2]/=r_magnitude;


  // We need to normalise v for r.v factor of specular component
  double v[3]={0,0,0};
  double magnitude= sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
  v[0]=-1*m[0]/magnitude;
  v[1]=-1*m[1]/magnitude;
  v[2]=-1*m[2]/magnitude;

  double r_dot_v=dot_product(v,r.direction);
 // cout<<"RdotV ORG: "<<r_dot_v<<endl;
  if(r_dot_v < 0)    {r_dot_v=0;}
  else if(r_dot_v>1) {r_dot_v=1;}

  //cout<<"RdotV: "<<r_dot_v<<endl;
  //cout<<endl;

  // Finally applying the Phong model of illumination for the expected output.
  double R = L.color[0] * ((diffuse_color[0] * l_dot_n) + specular_color[0] * pow(r_dot_v, shininess));

  double G = L.color[1] * ((diffuse_color[1] * l_dot_n) + specular_color[1] * pow(r_dot_v, shininess));

  double B = L.color[2] * ((diffuse_color[2] * l_dot_n) + specular_color[2] * pow(r_dot_v, shininess));
  
  //cout<<"Light color:"<<L.color[0]<<L.color[1]<<L.color[2]<<endl;
  //cout<<endl;
  // Clamp up the channels between 0 and 1 
  if(R < 0)       {R=0;}
  else if(R>1)  {R=1;}

  if(G < 0)      {G=0;}
  else if(G>1) {G=1;}

  if(B< 0)       {B=0;}
  else if(B>1) {B=1;}

  //cout<<"R: "<<R<<" G: "<<G<<" B: "<<B<<endl;
  Intensity finalcolor;
  finalcolor.r=R;
  finalcolor.g=G;
  finalcolor.b=B;

  return finalcolor;

}




//--------------------------------------------------------------------------------------
void display()
{

}

void keyboardFunc(unsigned char key, int x, int y) 
{
  switch (key) {
    
    case 27: // ESC key
        exit(0); // exit the program
        break;
        }
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);

  /* callback for keyboard */
    glutKeyboardFunc(keyboardFunc);
  init();
  glutMainLoop();
}
