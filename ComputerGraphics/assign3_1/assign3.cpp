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
#include <cstring>
#include <cmath>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>
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
#define fov 90.0
#define PI 3.14159265
const double MAX_DIST = -1e8;
const double BIAS =1e-10 ;

/* Enable Motion Blur */
bool motionBlur=0;

/* Enable light motion for animation */
bool lightMotion=0;

/* Enable soft shadows */
bool softShadow=0;
int num_areaLightComps=35;

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


struct point_3d
{
  double x;
  double y;
  double z;
};

struct triangle_test
{
  bool inside;
  double alpha;
  double beta;
  double gamma;
};

struct intersection
{
  point_3d m;
  double t_param;
  int index;
  int object_type;
  triangle_test t_test;
};



Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];

double ambient_light[3];

struct point_3d camera;

int num_triangles=0;
int num_spheres=0;
int num_lights=0;


//--------------------------------------------------------------------------------------


/* pixel handling function prototypes */
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

point_3d vector_diff(point_3d v1,point_3d v2);
point_3d scalar_division(point_3d c,double d);
double vector_magnitude(point_3d c);
point_3d vector_normalize(point_3d c);
bool equal_mag(point_3d c, point_3d d);
point_3d create_side(Vertex c,Vertex d);
point_3d cross_product(point_3d c,point_3d d);
point_3d cast_ray(point_3d origin,point_3d direction,double t_param);
triangle_test point_test(Triangle t, point_3d p);
double sphere_intersect(Sphere s,point_3d origin,point_3d direction);
point_3d sphere_normal(point_3d p, int ind);
double triangle_intersect(Triangle t,point_3d origin,point_3d direction);
point_3d triangle_interpolate(Triangle t,triangle_test test,int interpolation_type);
point_3d phong_model(point_3d p, int ind, int obj, triangle_test test,Light light, point_3d camera_origin);
point_3d intensity_map(int x, int y);


//----NEW-----

point_3d vector_diff(point_3d v1,point_3d v2)
{
  point_3d output;

  output.x=v1.x-v2.x;
  output.y=v1.y-v2.y;
  output.z=v1.z-v2.z;

  return output;
}

point_3d scalar_division(point_3d c,double d)
{
  point_3d output;
  output.x=0.0;
  output.y=0.0;
  output.z=0.0;

  if( abs(d)>1e-10)
  {
    output.x=c.x/d;
    output.y=c.y/d;
    output.z=c.z/d;
  }
 return output;
}

double vector_magnitude(point_3d c)
{
  return sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
}
 
point_3d vector_normalize(point_3d c)
{
  point_3d output;
  double mag;

  mag=vector_magnitude(c);
  output=scalar_division(c,mag);

  return output;
}

point_3d cross_product(point_3d c,point_3d d)
{
  point_3d output;

  output.x= c.y * d.z - d.y * c.z;
  output.y= d.x * c.z - c.x * d.z;
  output.z= c.x * d.y - c.y * d.x;

  return output;
} 

double dot_product(point_3d c, point_3d d)
{
  return (c.x * d.x + c.y * d.y + c.z * d.z);
}

point_3d create_side(Vertex c,Vertex d)
{
  point_3d output;

  output.x= c.position[0]  -  d.position[0] ;
  output.y= c.position[1]   - d.position[1] ;
  output.z= c.position[2]   - d.position[2] ;

  return output;
}


point_3d cast_ray(point_3d origin,point_3d direction,double t_param)
{
  point_3d output;

  output.x= origin.x + t_param * direction.x;
  output.y= origin.y + t_param * direction.y;
  output.z= origin.z + t_param * direction.z;

  return output;
}

bool equal_mag(point_3d c, point_3d d)
{
  bool output;
  if ((abs(c.x-d.x)<1e-10) && (abs(c.y-d.y)<1e-10) && (abs(c.z-d.z)<1e-10))
    output=1;
  else output=0;
   
  return output;
}


double sphere_intersect(Sphere s,point_3d origin,point_3d direction)
{
    double t;

    double x1 = origin.x - s.position[0];
    double y1 = origin.y - s.position[1];
    double z1 = origin.z - s.position[2];

    double a = 1 ;// Should be equal to 1 as it is magnitude of the direction vector. 
    double b = 2 * ((direction.x * x1) + (direction.y * y1) + (direction.z * z1));
    double c = (x1 * x1) + (y1 * y1) + (z1 * z1) - (s.radius * s.radius);

    double discriminant = (b * b) - (4 * a * c);
    if (discriminant < 0) t= -1;  // No real intersection

    // The intersection based on the discriminant value
    double t0 = (-1 * b - sqrt(discriminant)) / (2 * a); // point along ray which enters/exits the sphere
    double t1 = (-1 * b + sqrt(discriminant)) / (2 * a); // point along ray which exits/enters the sphere

    if (t0 < 0 && t1 < 0) {
      t=-1;
    }
    // Closer point if both values are positive and eliminate negative value.
    if (t0 < 1e-10) 
    {
      t = t1;
      
    }
    else if (t1 < 1e-10) {
      t = t0;
    }
    else {
      t = min(t0, t1);
    }

    return t;

}

triangle_test point_test(Triangle t, point_3d p)
{

  triangle_test test;
  Vertex P;
  // create the the point for the calculation of the barycentric cordinates
  P.position[0]=p.x;
  P.position[0]=p.y;
  P.position[0]=p.z;

  point_3d PA= create_side(P,t.vertex[0]);
  point_3d PB= create_side(P,t.vertex[1]);
  point_3d PC= create_side(P,t.vertex[2]);

  // Check for the cross products for normals with regards to P

  point_3d crossAB=cross_product(PA,PB);
  point_3d crossBC=cross_product(PB,PC);
  point_3d crossCA=cross_product(PC,PA);

  double alpha_un_norm = vector_magnitude(crossBC);
  double beta_un_norm = vector_magnitude(crossCA);
  double gamma_un_norm= vector_magnitude(crossAB);

  // Normalize the cross products 
  point_3d crossAB_norm=scalar_division(crossAB,gamma_un_norm);  
  point_3d crossBC_norm=scalar_division(crossBC,alpha_un_norm);  
  point_3d crossCA_norm=scalar_division(crossCA,beta_un_norm);  

  if ((alpha_un_norm+beta_un_norm+gamma_un_norm)> 1e-20 )  
  {
      test.alpha=alpha_un_norm/(alpha_un_norm+beta_un_norm+gamma_un_norm);//alpha
      test.beta=beta_un_norm/(alpha_un_norm+beta_un_norm+gamma_un_norm);//beta
      test.gamma=gamma_un_norm/(alpha_un_norm+beta_un_norm+gamma_un_norm);//gamma
  } 


    // Find if the normalized cross products are equal, since they are unit vectors now
    if(equal_mag(crossAB_norm,crossBC_norm) && equal_mag(crossAB_norm,crossCA_norm))
      test.inside=1;
    else test.inside=0;

    return test;  
}

double triangle_intersect(Triangle t,point_3d origin,point_3d direction)
{
  point_3d edge1,edge2;
  point_3d face_normal;
  double t_param=0;
  double dir=0;

  edge1=create_side(t.vertex[0],t.vertex[1]);
  edge2=create_side(t.vertex[0],t.vertex[2]);

  face_normal= cross_product(edge1,edge2);
  face_normal=vector_normalize(face_normal);

   dir = (-1)*(t.vertex[1].position[0] * face_normal.x
    + t.vertex[1].position[1] * face_normal.y
    + t.vertex[1].position[2] * face_normal.z);

  if (abs(dot_product(face_normal,direction))<1e-35) 
    t_param=-1;
  else
    t_param=(-1)*(dot_product(face_normal,origin)+dir)/(dot_product(face_normal,direction));

  return t_param;
}

point_3d sphere_normal(point_3d p, int ind)
{
  point_3d output;
 
  output.x=(p.x-spheres[ind].position[0])/spheres[ind].radius;
  output.y=(p.y-spheres[ind].position[1])/spheres[ind].radius;
  output.z=(p.z-spheres[ind].position[2])/spheres[ind].radius;

  return output;
}

// interpolation_type=0 ----> normals
// interpolation_type=1 ----> diffuse
// interpolation_type=2 ----> specular
point_3d triangle_interpolate(Triangle t,triangle_test test,int interpolation_type)
{
  point_3d p;
  
  if (interpolation_type==0)
  {
    p.x = test.alpha * t.vertex[0].normal[0]
         +test.beta * t.vertex[1].normal[0]
         +test.gamma * t.vertex[2].normal[0];

    p.y = test.alpha * t.vertex[0].normal[1]
         +test.beta * t.vertex[1].normal[1]
         +test.gamma * t.vertex[2].normal[1];

    p.z = test.alpha * t.vertex[0].normal[2]+
         test.beta * t.vertex[1].normal[2]
         +test.gamma * t.vertex[2].normal[2];
  }

  else if (interpolation_type==1)
  {
    p.x = test.alpha * t.vertex[0].color_diffuse[0]
         +test.beta * t.vertex[1].color_diffuse[0]
         +test.gamma * t.vertex[2].color_diffuse[0];

    p.y = test.alpha * t.vertex[0].color_diffuse[1]
         +test.beta * t.vertex[1].color_diffuse[1]
         +test.gamma * t.vertex[2].color_diffuse[1];

    p.z = test.alpha * t.vertex[0].color_diffuse[2]+
         test.beta * t.vertex[1].color_diffuse[2]
         +test.gamma * t.vertex[2].color_diffuse[2];
  }
  else if (interpolation_type==1)
  {
    p.x = test.alpha * t.vertex[0].color_specular[0]
         +test.beta * t.vertex[1].color_specular[0]
         +test.gamma * t.vertex[2].color_specular[0];

    p.y = test.alpha * t.vertex[0].color_specular[1]
         +test.beta * t.vertex[1].color_specular[1]
         +test.gamma * t.vertex[2].color_specular[1];

    p.z = test.alpha * t.vertex[0].color_specular[2]+
         test.beta * t.vertex[1].color_specular[2]
         +test.gamma * t.vertex[2].color_specular[2];
  }

  return p;
}
// p is the point of intersection.
point_3d phong_model(point_3d p, int ind, int obj, triangle_test test,Light light, point_3d camera_origin)
{
  point_3d n,l,v,r,kd,ks;
  point_3d intensity;

  double l_dot_n,r_dot_v,shininess_factor;


  // Convert to a point type
    l.x=light.position[0];
    l.y=light.position[1];
    l.z=light.position[2];

    l=vector_normalize(vector_diff(l,p));
    v=vector_normalize(vector_diff(camera_origin,p));

    // Sphere

    if(obj==1)
    {
      n=sphere_normal(p,ind);

      kd.x= spheres[ind].color_diffuse[0];
      kd.y= spheres[ind].color_diffuse[1];
      kd.z= spheres[ind].color_diffuse[2];


      ks.x=spheres[ind].color_specular[0];  
      ks.y=spheres[ind].color_specular[1];  
      ks.z=spheres[ind].color_specular[2];
      
      shininess_factor=spheres[ind].shininess; 

    }
    // Triangle
    else if(obj ==2)
    {

      n= vector_normalize(triangle_interpolate(triangles[ind],test,0));// using the interpolation_type for this case
      kd= triangle_interpolate(triangles[ind],test,1);
      ks= triangle_interpolate(triangles[ind],test,2);

      shininess_factor= test.alpha *  triangles[ind].vertex[0].shininess + test.beta *  triangles[ind].vertex[1].shininess +test.gamma *  triangles[ind].vertex[2].shininess; 

    }

  l_dot_n=dot_product(l,n);
  if (l_dot_n<0)
    l_dot_n=0;  
  else if (l_dot_n>1.f) 
    l_dot_n=1.f;
  
  // Find the reflection vector with l and n as r=2(l.n)n-l
  r.x=2*l_dot_n*n.x-l.x;
  r.y=2*l_dot_n*n.y-l.y;
  r.z=2*l_dot_n*n.z-l.z;

  // Compute dot product r.v, clamp it to 0-1 
  r_dot_v=dot_product(r,v);
  if (l_dot_n<0) 
    l_dot_n=0;  
  else if (l_dot_n>1.f) 
    l_dot_n=1.f;

  // Compute Intensity using the Phong equation
  intensity.x=light.color[0]*((kd.x)*l_dot_n+((ks.x)*pow((r_dot_v),(shininess_factor)))); // r
  intensity.y=light.color[1]*((kd.y)*l_dot_n+((ks.y)*pow((r_dot_v),(shininess_factor)))); // g 
  intensity.z=light.color[2]*((kd.z)*l_dot_n+((ks.z)*pow((r_dot_v),(shininess_factor)))); // b
  

  return intensity;
}

intersection final_intersection(point_3d p1,point_3d p2, int shadow)
{
  
  point_3d org,dir,a,b;
  intersection final;
  triangle_test test;

  double t_param,t_min=0.0,t_max,t_sphere,t_triangle;
  int ind=-1,obj=-1;

  if (shadow==0)
  {
    org=p1;// camera origin
  }
  else if(shadow==1)
  {
    org= p2;// light origin
  }

  dir.x= p1.x-p2.x;
  dir.y= p1.y-p2.y;
  dir.z= p1.z-p2.z;
  dir=vector_normalize(dir);

  // Sphere

  for (int i=0;i<num_spheres;i++)
  {
    t_sphere=sphere_intersect(spheres[i],org,dir); // Find t
    if (t_min==0 && t_sphere>1e-10) // Check for positive t 
    {
      t_min=t_sphere;
      ind=i;
      obj=1;
    }
    else if (t_sphere<=t_min && t_sphere>1e-10) // Minimum positive t of all spheres 
    {
      t_min=t_sphere;
      ind=i;
      obj=1;
    }
  }

  // Triangle
  for (int i=0;i<num_triangles;i++)
  {
    t_triangle=triangle_intersect(triangles[i],org,dir); // Find t
    a=cast_ray(org,dir,t_triangle); // Find the point corresponding to tT
    test=point_test(triangles[i],a); // Check if it lies in the triangle
    
    if (test.inside==1) // If inside the triangle
    { 
      if (t_min==0 && t_triangle>1e-5) // Check for positive t
      {
        t_min=t_triangle;
        ind=i;
        obj=2;
        if (shadow==0)
          b=a;
        final.t_test.alpha=test.alpha;  
        final.t_test.beta=test.beta;    
        final.t_test.gamma=test.gamma;   
      }
      else if (t_triangle<t_min && t_triangle>1e-5) // Minimum positive t of all triangles and spheres 
      {
        t_min=t_triangle;
        ind=i;
        obj=2;
        if (shadow==0)
          b=a;
        final.t_test.alpha=test.alpha;  
        final.t_test.beta=test.beta;    
        final.t_test.gamma=test.gamma;  
      }
    }
  }
  // SHADOW rendering 

  if (shadow==1) 
  {
    // Find t on the ray that corresponds to the light position
    if (dir.x!=0)
    {
      t_max=(p1.x-org.x)/dir.x;
    }
    else if (dir.y!=0)
    {
      t_max=(p1.y-org.y)/dir.y;
    }
    else if (dir.x!=0)
    {
      t_max=(p1.z-org.z)/dir.z;
    }
    else t_max=0;
  
    // The point should lie on the ray before the light
    // thus, t of the Point should be less than the t of the light  
    if (t_min>=t_max)
    {
      obj=-1;
      ind=-1;
    }
  }
  // If its not a shadow ray and the intersection is with circle
  // find the point of intersection, wasnt found so far for circle
  // If its a triangle, the point is already found
  else if ((t_min>=0) && (obj==1))
    b=cast_ray(org,dir,t_min); 
  
  final.m=b;
  final.t_param=t_min;
  final.index=ind;
  final.object_type=obj;
  
  return final;
}

point_3d intensity_map(int x, int y)
{
  double aspectRatio=(double)WIDTH / (double)HEIGHT;
  double x_right= aspectRatio * tan(fov* PI / 360);// Since angle is in radians 
  double x_left= -1*x_right;
  double y_top= tan(fov * PI / 360);
  double y_bottom= -1*y_top;
  double z=-1;// focal length fixed at 1 in the neagtive z direction.
  
  double frame_width= x_right-x_left;
  double frame_height= y_top- y_bottom;
  intersection object,shadow,soft_shadow;

  point_3d a,b,dir,light,light_shadow;
  point_3d img_color,temp,temp_n;

  int control=1;
    if (softShadow)
    control = num_areaLightComps+1;  

  

  point_3d black;// Setting up the background
  black.x=0.0;
  black.y=0.0;
  black.z=0.0;

  // frame to image mapping.

  a.x=(((double)x/(double)WIDTH)*2*x_right)-x_right;
  a.y=(((double)y/(double)HEIGHT)*2*y_top)-y_top;
  a.z=-1; 

  object=final_intersection(a,camera,0);

    // If it intersects any object only, enter the loop to find phong color
  if (object.index!=-1) 
  { 
    // Point of intersection
    b=object.m; 
    
    // Add the global ambient light first
    img_color.x+=ambient_light[0];
    img_color.y+=ambient_light[1];
    img_color.z+=ambient_light[2];
    
    // Find the contribution of each light source
    for (int h=0;h<num_lights;h+=control)
    {
      // Convert to a point type
      light.x=lights[h].position[0];
      light.y=lights[h].position[1];  
      light.z=lights[h].position[2];
    
      // Shoot shadow rays to the light from point and find intersection with the objects
      shadow=final_intersection(light,b,1);
    
      // If it doesnt intersect any object only, enter to find contribution of the light
      if ((shadow.index==-1))
      {
        // Soft Shadow Computation: shoot rays around the light and average
        if (softShadow) 
        {
          temp=black;
          
          for (int j=h;j<(h+(num_areaLightComps)+1);j++)
          {
            // Convert to a point type
            light_shadow.x=lights[j].position[0];
            light_shadow.y=lights[j].position[1];  
            light_shadow.z=lights[j].position[2];
          
            // Find intersection with Area light  
            soft_shadow=final_intersection(light_shadow,b,1);
          
            if (soft_shadow.index==-1)
            {
              // Calculate phong color for each of the lights
              temp_n=phong_model(b,object.index,object.object_type,object.t_test,lights[j],camera);
              
              // Sum the results
              temp.x+=temp_n.x;
              temp.y+=temp_n.y;
              temp.z+=temp_n.z;
            }
          }
          
          // Average the area light
          if (num_areaLightComps!=0)
            temp=scalar_division(temp,num_areaLightComps+1);
        }       
        else
          // Find the phong model color
          temp=phong_model(b,object.index,object.object_type,object.t_test,lights[h],camera);
        
        // Add to the pixel color
        img_color.x+=temp.x;
        img_color.y+=temp.y;
        img_color.z+=temp.z;
      }
    }

    // If pixel color goes above 1.f, clamp it to 1.f
    if (img_color.x > 1) img_color.x=1.f;
    if (img_color.y > 1) img_color.y=1.f;
    if (img_color.z > 1) img_color.z=1.f;
  }
  else img_color=black;
  
  return img_color;
}

void final_color()
{
  unsigned int x,y;
  point_3d img_color;
  for(x=0; x<WIDTH; x++)
  {
    for(y=0;y < HEIGHT;y++)
    {
      img_color=intensity_map(x,y);
      plot_pixel_jpeg(x,y,abs(img_color.x)*255,abs(img_color.y)*255,abs(img_color.z)*255);
    }
  }
}




//=============================Starter code =====================
// change this draw scene according to newer functions. 
void draw_scene()
{
  //raycast();
  unsigned int x,y;
  glPointSize(2.0);
 
  //simple output

   glBegin(GL_POINTS);
  for(x=0; x<WIDTH; x++)
  {
    for(y=0;y < HEIGHT;y++)
    {  
      plot_pixel_display(x,y,buffer[HEIGHT-y-1][x][0],buffer[HEIGHT-y-1][x][1],buffer[HEIGHT-y-1][x][2]);
    }
  }
  glEnd();
    
    if(!motionBlur && !lightMotion)
    glFlush();

  printf("Done!\n");

  if(!motionBlur && !lightMotion)
  fflush(stdout);

}

/* Plot to the GL window */
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

/* Put the values in a buffer */
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

/* Choose to draw to window or save to jpeg */
void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
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

  if (motionBlur)
  {
    for (int i=HEIGHT-1; i>=0; i--) 
    {
      glReadPixels(0, HEIGHT-1-i, WIDTH, 1, GL_RGB, GL_UNSIGNED_BYTE,
                 &in->pix[i*in->nx*in->bpp]);
    }
  }

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


//--------------------------------------------------------------------------------------
void display()
{
  static int once=0;
  static int count=0;
  static int n=15;
  static double blur=0;
  static int k=0;
  string fnameMB; 
  string fname; 

  // If light is moving, every frame is refreshed
  if (lightMotion)
    glClear(GL_COLOR_BUFFER_BIT);// | GL_DEPTH_BUFFER_BIT);

  // Camera Position for the scene
  camera.x=0.0;
  camera.y=0.0;
  camera.z=0.0;

  // Reset the world
  glLoadIdentity();

  // First time writing to color Buffer and lightMotion computation
  if((lightMotion && once<100) || (!once && !lightMotion))
  {
    // First time filling the color buffer by ray Tracing 
    final_color();
    
    // If motionBlur=1, dont draw it here
    if (!motionBlur)
      draw_scene();

    // Save lightMotion screenshots
    if (mode==MODE_JPEG & !motionBlur & lightMotion)
    {
      k++;
      
      if (k>0 && k<100) 
      {
        stringstream ss;
          ss << k;
          ss >> fnameMB;
  
        if (k<10)
          fnameMB.insert(0,"00");
        else if (k<100)
          fnameMB.insert(0,"0");
        
        fname="lM_"+fnameMB+".jpg";
        
        strcpy(filename, fname.c_str());
        
        save_jpg();
      }
    }
    else if (mode==MODE_JPEG && !motionBlur && !lightMotion)
      save_jpg();
    
    // Move the lights in the scene by an offset
    if (lightMotion && !motionBlur)
    {
      for (int l=0;l<num_lights;l++)
      {
        lights[l].position[0]+=(once*0.1-2);
        lights[l].position[0]+=(once*0.1-2);
      }
    }
  }
  once++;
  
  // Swap buffers and Redisplay for lightMotion 
  if (lightMotion && once<100) 
  {
    glutSwapBuffers();
    glutPostRedisplay();
  }   

  // Motion Blur computation
  if (motionBlur)
  {
    // variable to keep a track of iterations
    count++;

    // Measure of the movement/blur
    blur+=2;
    
    // Move the scene
    glTranslatef(blur,0,0);
    
    // Redraw the scene
    draw_scene();

    // Load into the accumulation buffer and keep accumulating a bunch of frames
    if (once==1)
      glAccum(GL_LOAD,0.5/n); //0.5/n
    else 
      glAccum(GL_ACCUM,0.5/n);

    // Stop the blurring when count reaches 2*n
    if (count<2*n) 
    { 
      // Copy the accumulation buffer contents
      glAccum(GL_RETURN,1.0f);

      // Swap buffers and Redisplay
      glutSwapBuffers();
      glutPostRedisplay();

      // Save motionBlur images
      if (mode==MODE_JPEG)
      {
        k++;

        if (k>0 && k<100) 
        {
          stringstream ss;
            ss << k;
            ss >> fnameMB;
  
          if (k<10)
            fnameMB.insert(0,"00");
          else if (k<100)
            fnameMB.insert(0,"0");
          
          fname="mB_"+fnameMB;
          fname+= ".jpg";
  
          strcpy(filename, fname.c_str());

          save_jpg();
        }
        glClear(GL_COLOR_BUFFER_BIT);
      }
    }

    // Bunch of frames to accumulate, the accumulation buffer is reloaded
    // when once=1
    if (once>=n) {
      once=1;
    }
  }


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
   glClear(GL_ACCUM_BUFFER_BIT);

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
 if (argc!=2 && argc!=3 && argc!=6)
  {  
    // For animation, MUST give all flags
    printf ("\nUsage1: %s <scenefile> <motionBlurFlag> <lightMotionFlag> <softShadowFlag> <jpegFlag>\n", argv[0]);
    printf ("\nFlag=1, enable\nFlag=0, disable");
    printf ("\nPlease enable only one of motion Blur, light motion or soft shadows at one time\n");
    printf ("\neg for motion blur: ./assign3 table.scene 0 1 0 0");
    printf ("\neg for light motion: ./assign3 spheres.scene 1 0 0 0");
    printf ("\neg for soft shadows: ./assign3 test2.scene 0 0 1 0\n\n");
    // For simple ray tracing
    printf ("\nUsage2: %s <scenefile> [jpegname]\n\n", argv[0]);
    exit(0);
  }

  if(argc == 6)
  {
    if ( atoi(argv[4])==1 )
      mode = MODE_JPEG;
    else 
      mode = MODE_DISPLAY;
    motionBlur = atoi(argv[2]);
    lightMotion = atoi(argv[3]);
    softShadow = atoi(argv[4]);
    filename = argv[5];
  }
  else if (argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  // Double buffering enabled in case of animations 
  if (motionBlur || lightMotion)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ACCUM); 
  else
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE | GLUT_ACCUM);
 
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();

}
