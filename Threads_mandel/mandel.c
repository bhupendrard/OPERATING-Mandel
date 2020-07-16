/*

	Name: Bhupendra Ramdam
	ID:   10013370027

*/
#include "bitmap.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
//creating struct for storing the respective values
struct mystruct
{
	struct bitmap *bm;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int max;
	int begin_height;
	int end_height;
};

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max ,int thread_num);
void *arg_func(void *arg);
void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}



int main( int argc, char *argv[] )
{
	char c;

	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int    max = 1000;
	int    thread_num = 1;

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:o:n:h"))!=-1)
	{
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'n':
				thread_num = atoi(optarg);
				break;
			case 'h':
				show_help();
				exit(1);
				break;

		}
	}

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d thread number=%d outfile=%s\n",xcenter,ycenter,scale,max,thread_num,outfile);


	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,thread_num);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}


/*
this is the void argument function that will be used while doing pthread
Scale the image with the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/
void *arg_func(void *arg)
{
	//we are creating the pointer to the struct inorder to get the value from the struct
	struct mystruct* newstr;
	newstr=(struct mystruct*)arg; //typecasting the void type to the struture type
	int i,j;
	int width = bitmap_width(newstr->bm);
	int height = bitmap_height(newstr->bm);

	int b_heigh= newstr->begin_height;
	int e_heigh= newstr->end_height;


	for(j=b_heigh;j<e_heigh;j++)
	{
		for(i=0;i<width;i++)
		{

			// Determine the point in x,y space for that pixel.
			double x = newstr->xmin + i*(newstr->xmax-newstr->xmin)/width;
			double y = newstr->ymin + j*(newstr->ymax-newstr->ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,newstr->max);

			// Set the pixel in the bitmap.
			bitmap_set(newstr->bm,i,j,iters);
		}
	}
	return 0;//we have void* to return nothing we do return 0.
}


//Compute an entire Mandelbrot image, writing each point to the given bitmap.
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int thread_num)
{
	//creating array of struct
	struct mystruct struct_array[thread_num];
	int i;
	pthread_t tid[thread_num];
	int previous=0;

	int height = bitmap_height(bm);
	// seperate images into band of images
	int new_height= height/thread_num;
	//int remainder = height%new_height;
	//will give result for the first portion of images depending on thread number
	struct_array[0].begin_height=0;
	struct_array[0].end_height=new_height;
	//looping for rest portion if got any
	for(i = 0; i < thread_num; i++)
	{

		struct_array[i].bm = bm;
		struct_array[i].xmin = xmin;
		struct_array[i].xmax = xmax;
		struct_array[i].ymin = ymin;
		struct_array[i].ymax = ymax;
		struct_array[i].max = max;

		if(i!=0)
		{
			//save size of the last band of images so that we can use it to track the upcomming band image
			previous=struct_array[i-1].end_height;
			struct_array[i].begin_height=previous;
			struct_array[i].end_height=previous+new_height;
			//to avoid the lost pixel by the above computation
			// changing the value of end height on the last thread for image completion
			if(i==thread_num-1)
				struct_array[i].end_height=height;

		}
		//creating the thread by accessing the structure
		pthread_create(&tid[i],NULL,arg_func,(void*) &struct_array[i]);

	}
	for(i=0;i< thread_num;i++)
	{
		//wait until the thread is terminated
		pthread_join(tid[i],NULL);
	}

}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max )
	{

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
