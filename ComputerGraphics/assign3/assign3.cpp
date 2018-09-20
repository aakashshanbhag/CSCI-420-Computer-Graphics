
/*
CSCI 420: Computer Graphics.
Assignment 3: Raytracer.
Name: Aakash Shanbhag.
USC ID:3205699915.
*/

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

/* Global Defines */
#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 100

/* Filename handle */
char *filename=0;

/* Different display modes */
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

/* Size of the window */
#define WIDTH 640
#define HEIGHT 480

/* Field of view of the cameraera */
#define fov 60.0// change to any angle needed.
#define PI 3.14159265

/* Enable Motion Blur */
bool motionBlur=0;

/* Enable light motion for animation */
bool lightMotion=0;

/* Enable soft shadows */
bool softShadow=0;
int total_lights=35;



/* Color buffer to hold the rgb values each 8-bit */
unsigned char buffer[HEIGHT][WIDTH][3];

/* Program-specific structs for easier access */
struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct point_3d
{
	double x;
	double y;
	double z;
};


struct triangle_test
{
	bool in; // =0 if outside , =1 if inside
	double alpha;
	double beta;
	double gamma; //alpha,beta,gamma
};


struct intersection
{
	point_3d p; // intersection point_3d
	double t; // final t param
	int tID; // index of (object=triangle/sphere)
	int tObj; // =1 if n sphere, = 2 if triangle
	triangle_test iO; //  barycentric co-ordinates if triangle
};

/* Objects */
typedef struct _Triangle
{
  struct Vertex v[3];
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

/* Global variables */
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];

double ambient_light[3];

struct point_3d camera;

int num_triangles=0;
int num_spheres=0;
int num_lights=0;


// Plotting helper functions
// Plotting pixels in different cases as unsigned char.
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// Vector creation and manipulation helper functions.
point_3d vector_diff(point_3d v1,point_3d v2);
point_3d scalar_division(point_3d c,double d);
double vector_magnitude(point_3d c);
point_3d vector_normalize(point_3d c);
point_3d cross_product(point_3d c,point_3d d);
double dot_product(point_3d c, point_3d d);
bool equal_mag(point_3d c,point_3d d);
point_3d cast_ray(point_3d origin,point_3d direction,double t_param);

// Sphere specific helper functions
double sphere_intersect(Sphere s,point_3d origin,point_3d direction);
point_3d sphere_normal(point_3d p, int ind);

// Triangle specifc helper functions 
double triangle_intersect(Triangle t,point_3d origin,point_3d direction);
triangle_test point_test(Triangle triangle,point_3d p);
point_3d triangle_interpolate(Triangle triangle,triangle_test iO,int ID);

// Phong illumination model helper function  
point_3d phong_model(point_3d p,int tID,int tObj,triangle_test iO,Light light,point_3d cameraera);

// Difference between 2 points 
point_3d vector_diff(point_3d v1,point_3d v2)
{
  point_3d output;

  output.x=v1.x-v2.x;
  output.y=v1.y-v2.y;
  output.z=v1.z-v2.z;

  return output;
}

// Scalar division for utilisation in unit vector creation.
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

// Vector magnitude 
double vector_magnitude(point_3d c)
{
  return sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
}

//Unit vector creation
point_3d vector_normalize(point_3d c)
{
  point_3d output;
  double mag;

  mag=vector_magnitude(c);
  output=scalar_division(c,mag);

  return output;
}

// Cross product between two vectors.
point_3d cross_product(point_3d c,point_3d d)
{
  point_3d output;

  output.x= c.y * d.z - d.y * c.z;
  output.y= d.x * c.z - c.x * d.z;
  output.z= c.x * d.y - c.y * d.x;

  return output;
} 

// Dot product between 2 vectors
double dot_product(point_3d c, point_3d d)
{
  return (c.x * d.x + c.y * d.y + c.z * d.z);
}
// Equality check between vectors.
bool equal_mag(point_3d c,point_3d d)
{
	bool output;
	
	if ((abs(c.x-d.x)<1e-10) && (abs(c.y-d.y)<1e-10) && (abs(c.z-d.z)<1e-10))
		output=1;
	else output=0;

	return output;
}

// Equation of ray with origin and direction with t_param being the scalar multiplicative.
point_3d cast_ray(point_3d origin,point_3d direction,double t_param)
{
  point_3d output;

  output.x= origin.x + t_param * direction.x;
  output.y= origin.y + t_param * direction.y;
  output.z= origin.z + t_param * direction.z;

  return output;
}

// Finding the t_param for intersection in case of the sphere.
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

// Creating vector to be utilised for sides in triangles.
point_3d create_side(Vertex v1, Vertex v2)
{
	point_3d AB;
	
	AB.x=v1.position[0]-v2.position[0];
	AB.y=v1.position[1]-v2.position[1];
	AB.z=v1.position[2]-v2.position[2];
	
	return AB;
}

// Check if point lies inside trinagle by the barycentric tests.
triangle_test point_test(Triangle triangle,point_3d p)
{
	triangle_test iO;
	Vertex P;

	point_3d PA,PB,PC;
		
	// Convert point_3d to a vertex	
	P.position[0]=p.x;
	P.position[1]=p.y;
	P.position[2]=p.z;
	
	// Find the sides of the triangle
	PA=create_side(P,triangle.v[0]);
	PB=create_side(P,triangle.v[1]);
	PC=create_side(P,triangle.v[2]);

	// Find the un-normalized cross products of the sides
	point_3d crossAB=cross_product(PA,PB);
  	point_3d crossBC=cross_product(PB,PC);
  	point_3d crossCA=cross_product(PC,PA);
	
	// Magnitudes of the cross products
	
	double alpha_un_norm = vector_magnitude(crossBC);
  	double beta_un_norm = vector_magnitude(crossCA);
  	double gamma_un_norm= vector_magnitude(crossAB);
	
	// Find the barycentric co-ords

	if ((alpha_un_norm+beta_un_norm+gamma_un_norm)> 1e-20 )  
  	{
      iO.alpha=alpha_un_norm/(alpha_un_norm+beta_un_norm+gamma_un_norm);//alpha
      iO.beta=beta_un_norm/(alpha_un_norm+beta_un_norm+gamma_un_norm);//beta
      iO.gamma=gamma_un_norm/(alpha_un_norm+beta_un_norm+gamma_un_norm);//gamma
 	}

	// Normalize the cross products	
	point_3d crossAB_norm=scalar_division(crossAB,gamma_un_norm);  
  	point_3d crossBC_norm=scalar_division(crossBC,alpha_un_norm);  
  	point_3d crossCA_norm=scalar_division(crossCA,beta_un_norm);  	

	// Find if the normalized cross products are equal, since they are unit vectors now
	if(equal_mag(crossAB_norm,crossBC_norm) && equal_mag(crossAB_norm,crossCA_norm))
		iO.in=1;
	else iO.in=0;

	return iO;				
}

// Estimating the t_parameter in case of a triangle 
double triangle_intersect(Triangle t,point_3d origin,point_3d direction)
{
  point_3d edge1,edge2;
  point_3d face_normal;
  double t_param=0;
  double dir=0;

  edge1=create_side(t.v[0],t.v[1]);
  edge2=create_side(t.v[0],t.v[2]);

  face_normal= cross_product(edge1,edge2);
  face_normal=vector_normalize(face_normal);

   dir = (-1)*(t.v[1].position[0] * face_normal.x
    + t.v[1].position[1] * face_normal.y
    + t.v[1].position[2] * face_normal.z);

  if (abs(dot_product(face_normal,direction))<1e-35) 
    t_param=-1;
  else
    t_param=(-1)*(dot_product(face_normal,origin)+dir)/(dot_product(face_normal,direction));

  return t_param;
}

// Calculation of the sphere normal
point_3d sphere_normal(point_3d p, int ind)
{
  point_3d output;
 
  output.x=(p.x-spheres[ind].position[0])/spheres[ind].radius;
  output.y=(p.y-spheres[ind].position[1])/spheres[ind].radius;
  output.z=(p.z-spheres[ind].position[2])/spheres[ind].radius;

  return output;
}

/*	Normalisation of the individual components on the basis of ID
		 ID = 0, triangle_interpolate the normal
			= 1, triangle_interpolate the diffuse component,
			= 2, triangle_interpolate the specular component 
*/
point_3d triangle_interpolate(Triangle triangle,triangle_test iO,int ID)
{
	point_3d P;
	
	if (ID==0)
	{
		P.x = iO.alpha*triangle.v[0].normal[0]
				 +iO.beta*triangle.v[1].normal[0]
				 +iO.gamma*triangle.v[2].normal[0];

		P.y = iO.alpha*triangle.v[0].normal[1]
				 +iO.beta*triangle.v[1].normal[1]
				 +iO.gamma*triangle.v[2].normal[1];

		P.z = iO.alpha*triangle.v[0].normal[2]
				 +iO.beta*triangle.v[1].normal[2]
				 +iO.gamma*triangle.v[2].normal[2];
	}

	else if (ID==1)
	{
		P.x = iO.alpha*triangle.v[0].color_diffuse[0]
				 +iO.beta*triangle.v[1].color_diffuse[0]
				 +iO.gamma*triangle.v[2].color_diffuse[0];

		P.y = iO.alpha*triangle.v[0].color_diffuse[1]
				 +iO.beta*triangle.v[1].color_diffuse[1]
				 +iO.gamma*triangle.v[2].color_diffuse[1];

		P.z = iO.alpha*triangle.v[0].color_diffuse[2]
				 +iO.beta*triangle.v[1].color_diffuse[2]
				 +iO.gamma*triangle.v[2].color_diffuse[2];
	}

	else if (ID==2)
	{
		P.x = iO.alpha*triangle.v[0].color_specular[0]
				 +iO.beta*triangle.v[1].color_specular[0]
				 +iO.gamma*triangle.v[2].color_specular[0];

		P.y = iO.alpha*triangle.v[0].color_specular[1]
				 +iO.beta*triangle.v[1].color_specular[1]
				 +iO.gamma*triangle.v[2].color_specular[1];

		P.z = iO.alpha*triangle.v[0].color_specular[2]
				 +iO.beta*triangle.v[1].color_specular[2]
				 +iO.gamma*triangle.v[2].color_specular[2];
	}

	return P;
}

// Phong model used for shading with all the parameters as discussed in lect 8.1
point_3d phong_model(point_3d p,int tID,int tObj,triangle_test iO,Light light,point_3d cameraera)
{
	point_3d n,l,v,r,kd,ks;
	point_3d intensity; 
	double l_dot_n,r_dot_v,alpha;

	// Convert to a point_3d type
	l.x=light.position[0];
	l.y=light.position[1];
	l.z=light.position[2];

	// Find the direction vector from given point_3d to the given light
	l=vector_normalize(vector_diff(l,p));
	
	// Find the direction vector from the given point_3d to the given camera position
	v=vector_normalize(vector_diff(cameraera,p));

	// Find the normal, kd, ks, shininess based on which object the point_3d intersects
	if (tObj==1) // if sphere
	{
		n=sphere_normal(p,tID);
		
		kd.x=spheres[tID].color_diffuse[0];	
		kd.y=spheres[tID].color_diffuse[1];	
		kd.z=spheres[tID].color_diffuse[2];	
		
		ks.x=spheres[tID].color_specular[0];	
		ks.y=spheres[tID].color_specular[1];	
		ks.z=spheres[tID].color_specular[2];
		
		alpha=spheres[tID].shininess;	
	}

	else if (tObj==2) // if triangle
	{
		n=vector_normalize(triangle_interpolate(triangles[tID],iO,0));
		kd=triangle_interpolate(triangles[tID],iO,1);
		ks=triangle_interpolate(triangles[tID],iO,2);
		
		alpha=iO.alpha*triangles[tID].v[0].shininess
				 +iO.beta*triangles[tID].v[1].shininess
				 +iO.gamma*triangles[tID].v[2].shininess;
	}

	//  Clamp it to 0-1	
	l_dot_n=dot_product(l,n);
	if (l_dot_n<0)
		l_dot_n=0;	
	else if (l_dot_n>1.f) 
		l_dot_n=1.f;
	
	// r=2(l.n)n-l
	r.x=2*l_dot_n*n.x-l.x;
	r.y=2*l_dot_n*n.y-l.y;
	r.z=2*l_dot_n*n.z-l.z;

	//  Clamp it to 0-1	
	r_dot_v=dot_product(r,v);
	if (r_dot_v<0) 
		r_dot_v=0;	
	else if (r_dot_v>1.f) 
		r_dot_v=1.f;

	// Compute Intensity using the phong_model equation
	intensity.x=light.color[0]*((kd.x)*l_dot_n+((ks.x)*pow((r_dot_v),(alpha)))); // r
	intensity.y=light.color[1]*((kd.y)*l_dot_n+((ks.y)*pow((r_dot_v),(alpha)))); // g 
	intensity.z=light.color[2]*((kd.z)*l_dot_n+((ks.z)*pow((r_dot_v),(alpha)))); // b
	return intensity;
}

//Estimate the final intersections based on the location of the point in the specific region.
intersection final_intersection(point_3d p1,point_3d p2,int shadow)
{
	point_3d p,q,org,dir;
	intersection final;
	triangle_test iO;	
	double t,t_min,t_max,t_sphere,t_triangle;
	int tID,tObj;
	
	// Initialize t, index and object
	t_min=0;
	tID=-1;
	tObj=-1;	

	// Change origin 	
	if (shadow==0)		
		org=p1; // p2 is camera
	else if (shadow==1)
		org=p2; // p1 is light
	
	// Find the direction vector of the ray	
	dir.x=p1.x-p2.x;
	dir.y=p1.y-p2.y;
	dir.z=p1.z-p2.z;
	dir=vector_normalize(dir);

	
	// First, with spheres	
	for (int i=0;i<num_spheres;i++)
	{
		t_sphere=sphere_intersect(spheres[i],org,dir); // Find t
		if (t_min==0 && t_sphere>1e-10) // Check for positive t 
		{
			t_min=t_sphere;
			tID=i;
			tObj=1;
		}
		else if (t_sphere<=t_min && t_sphere>1e-10) // Minimum positive t of all spheres 
		{
			t_min=t_sphere;
			tID=i;
			tObj=1;
		}
	}
	// with Triangles
	for (int i=0;i<num_triangles;i++)
	{
		t_triangle=triangle_intersect(triangles[i],org,dir); // Find t
		p=cast_ray(org,dir,t_triangle); // point of intersection
		iO=point_test(triangles[i],p); // Check if it lies in the triangle
		
		if (iO.in==1) // If inside the triangle
		{	
			if (t_min==0 && t_triangle>1e-5) // Check for positive t
			{
				t_min=t_triangle;
				tID=i;
				tObj=2;
				if (shadow==0)
					q=p;
			  final.iO.alpha=iO.alpha;	
			  final.iO.beta=iO.beta;	
			  final.iO.gamma=iO.gamma;	
			}
			else if (t_triangle<t_min && t_triangle>1e-5) // Minimum positive t of all triangles and spheres 
			{
				t_min=t_triangle;
				tID=i;
				tObj=2;
				if (shadow==0)
					q=p;
			  final.iO.alpha=iO.alpha;	
			  final.iO.beta=iO.beta;	
			  final.iO.gamma=iO.gamma;	
			}
		}
	}	

	// Check if the intersection point_3d lies beyond the light in case of a shadow ray
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
	
		// The point_3d should lie on the ray before the light
		// thus, t of the point_3d should be less than the t of the light	
		if (t_min>=t_max)
		{
			tObj=-1;
			tID=-1;
		}
	}
	// If its not a shadow ray and the intersection is with circle
	// find the point_3d of intersection, wasnt found so far for circle
	// If its a triangle, the point_3d is already found
	else if ((t_min>=0) && (tObj==1))
		q=cast_ray(org,dir,t_min); 
	
	final.p=q;
	final.t=t_min;
	final.tID=tID;
	final.tObj=tObj;
	return final;
}

// Corresponding intensity map after intersections have been finalised.
point_3d intensity_map(int x, int y) 
{
		 double aspectRatio=(double)WIDTH / (double)HEIGHT;
		/* Max co-ords in world space based on fov */
		double y_top=tan((double) PI*fov/(2*180));
		double x_right=y_top*(aspectRatio);
		
		point_3d p,q,dir,light,lightS;
		point_3d black,img_intensity,temp,tempN;
		intersection object_region,shadow_region, soft_region;
		int control=1;
		
		// Counter to loop through lights
		if (softShadow)
			control = total_lights+1;	

		// Default color	
		black.x=0.0;
		black.y=0.0;
		black.z=0.0;
		img_intensity=black;

		// Convert pixel co-ordinates to world co-ordinates	
		p.x=(((double)x/(double)WIDTH)*2*x_right)-x_right;
		p.y=(((double)y/(double)HEIGHT)*2*y_top)-y_top;
		p.z=-1;	

		// Shoot rays from pixel point_3d p to intersect objects
		object_region=final_intersection(p,camera,0);
		
		// If it intersects any object only, enter the loop to find phong_model color
		if (object_region.tID!=-1) 
		{	
			// point_3d of intersection
			q=object_region.p; 
			
			// Add the global ambient light first
			img_intensity.x+=ambient_light[0];
			img_intensity.y+=ambient_light[1];
			img_intensity.z+=ambient_light[2];
			
			// Find the contribution of each light source
			for (int h=0;h<num_lights;h+=control)
			{
				// Convert to a point_3d type
				light.x=lights[h].position[0];
				light.y=lights[h].position[1];  
				light.z=lights[h].position[2];
			
				// Shoot shadow rays to the light from point_3d and find intersection with the objects
				shadow_region=final_intersection(light,q,1);
			
				// If it doesnt intersect any object only, enter to find contribution of the light
				if ((shadow_region.tID==-1))
				{
					// Soft Shadow Computation: shoot rays around the light and average
					if (softShadow) 
					{
						temp=black;
						
						for (int j=h;j<(h+(total_lights)+1);j++)
						{
							// Convert to a point_3d type
							lightS.x=lights[j].position[0];
							lightS.y=lights[j].position[1];  
							lightS.z=lights[j].position[2];
						
							// Find intersection with Area light	
							soft_region=final_intersection(lightS,q,1);
						
							if (soft_region.tID==-1)
							{
								// Calculate phong_model color for each of the lights
								tempN=phong_model(q,object_region.tID,object_region.tObj,object_region.iO,lights[j],camera);
								
								// Sum the results
								temp.x+=tempN.x;
								temp.y+=tempN.y;
								temp.z+=tempN.z;
							}
						}
						
						// Average the area light
						if (total_lights!=0)
							temp=scalar_division(temp,total_lights+1);
					}				
					else
						// Find the phong_model model color
						temp=phong_model(q,object_region.tID,object_region.tObj,object_region.iO,lights[h],camera);
					
					// Add to the pixel color
					img_intensity.x+=temp.x;
					img_intensity.y+=temp.y;
					img_intensity.z+=temp.z;
				}
			}

			// If pixel color goes above 1.f, clamp it to 1.f
			if (img_intensity.x > 1) img_intensity.x=1.f;
			if (img_intensity.y > 1) img_intensity.y=1.f;
			if (img_intensity.z > 1) img_intensity.z=1.f;
		}
		else img_intensity=black;
		
		return img_intensity;
}

/* Fill the color buffer with all the computed values for a particular scene */
void final_color()
{
  unsigned int x,y;
	point_3d img_intensity;
 
	// For each pixel, find the color 
	for(x=0; x<WIDTH; x++)
  {
    for(y=0;y < HEIGHT;y++)
    {
			// Find the color
      img_intensity=intensity_map(x,y);

			// Put the color values in the buffer
			plot_pixel_jpeg(x,y,abs(img_intensity.x)*255,abs(img_intensity.y)*255,abs(img_intensity.z)*255);
    }
  }
}

/* With the values in the color buffer, draw them to the GL window */
void draw_scene()
{
  unsigned int x,y;
 glPointSize(2.0);
  
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
 		glFlush(); // Comment this
  printf("Done!\n"); 
	if(!motionBlur && !lightMotion)
		fflush(stdout); // Comment this
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

/* save the color buffer to jpeg */
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
	else
		memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  
	if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

/* Reading the scene files and error checking */
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
  Light l,l1;
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
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
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
		
		/* adding some extra lights to make it an area light */
		l1=l;
		if (softShadow)
		{
			for (int aL=1;aL<(total_lights+1);aL++)
			{
				l1.position[2]+=(aL*0.001-0.02);
	  		lights[num_lights++] = l1;
			}
		}
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}
 
void idle()
{
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


void display()
{
  static int once=0;
  static int count=0;
  static int n=15;
  static double blur=0;
  static int k=0;
  string name_anim;	
  string anim_final;	

	// Refresh lights every frame.
	if (lightMotion)
		glClear(GL_COLOR_BUFFER_BIT);

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
			// Just for the correct nomenclature.
			k++;
			
			if (k>0 && k<100) 
			{
				stringstream ss;
    			ss << k;
    			ss >> name_anim;
	
				if (k<10)
					name_anim.insert(0,"00");
				else if (k<100)
					name_anim.insert(0,"0");
				
				anim_final="light_motion_"+name_anim+".jpg";
				
				strcpy(filename, anim_final.c_str());
				
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
		
		count++;
		blur+=2;
		
		// Move the scene
		glTranslatef(blur,0,0);
		
		// Redraw the scene
		draw_scene();

		// Stack to keep all the frames with the accumulation buffer.
		if (once==1)
			glAccum(GL_LOAD,0.5/n); //0.5/n
		else 
			glAccum(GL_ACCUM,0.5/n);

		// Stop  blurring when count reaches 2*n
		if (count<2*n) 
		{ 
			//Local copy storage
			glAccum(GL_RETURN,1.0f);

			// Swap buffers and Redisplay every frame for refresh
			glutSwapBuffers();
			glutPostRedisplay();

			// Stroage and nomenclature motionBlur images
			if (mode==MODE_JPEG)
			{
				k++;
				cout<<"here"<<endl;
				if (k>0 && k<100) 
				{
					stringstream ss;
    				ss << k;
    				ss >> name_anim;
	
					if (k<10)
						name_anim.insert(0,"00");
					else if (k<100)
						name_anim.insert(0,"0");
					
					anim_final="trial_"+name_anim;
					anim_final+= ".jpg";
	
					strcpy(filename, anim_final.c_str());

					save_jpg();
				}
				glClear(GL_COLOR_BUFFER_BIT);
			}
		}

		// Reload the accumulation buffer after each frame
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


int main (int argc, char ** argv)
{
	
  if (argc!=2 && argc!=3 && argc!=7)
  {  
		printf ("\nUsage1: %s <scenefile>  <softShado`wFlag> <lightMotionFlag> <motionBlurFlag> <jpegname> <jpegFlag>\n", argv[0]);
		printf("\nFor simple ray tracing\n");
    	printf ("\nUsage2: %s <scenefile> [jpegname]\n\n", argv[0]);  
  }

  if(argc == 7)
  {
  	if ( atoi(argv[6])==1 )
			mode = MODE_JPEG;
		else 
			mode = MODE_DISPLAY;
		
		softShadow = atoi(argv[2]);
		lightMotion = atoi(argv[3]);
		motionBlur = atoi(argv[4]);
		
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
  int window = glutCreateWindow("Ray Tracer-AakashShanbhag");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glutKeyboardFunc(keyboardFunc);
  init();
  glutMainLoop();
}