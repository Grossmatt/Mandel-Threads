#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <time.h>


int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int num_threads);
void * threadbuilder();

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
        printf("-n <threads> The number of pthreads to run the programs. (default = 0)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}


// Struct to pass into the pthread function as to have access to all the var info.
struct image_vars {
  double xmin_var;
  double xmax_var;
  double ymin_var;
  double ymax_var;
  int width_var;
  int height_var;
  int j_var;
  int max_var;
  struct bitmap *bm;
};



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
        int    num_threads = 1;
        time_t start;
        time_t end;
        double time_var;




        start = time(NULL);
        printf(ctime(&start));


	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:n:m:o:h"))!=-1) {
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
                        case 'n':
                                num_threads = atoi(optarg);
                                break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}



	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n",xcenter,ycenter,scale,max,outfile);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,num_threads);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}



        // Calculating the time differences for the graphs
        end = time(NULL);
        printf(ctime(&end));
        
        time_var = difftime(end, start);
        printf("%.f seconds have passed.\n", time_var);

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int num_threads)
{
	int i,j;
	int rc;
        pthread_t *my_thread;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

        struct image_vars *mybit =  malloc(sizeof(struct image_vars));

        // Assinging all the vars passed into computer_image into their respective vars in the struct.

        mybit->ymin_var = ymin;
        mybit->ymax_var = ymax;
        mybit->xmax_var = xmax;
        mybit->xmin_var = xmin;
        mybit->width_var = width;
        mybit->height_var = height;
        mybit->max_var = max;
        mybit->bm = bm;


        my_thread = malloc(num_threads * (sizeof(pthread_t)));


	// Creating the array of threads to send into the void * function then once completed joining them again.


	for (j = 0; j < height; j = j + num_threads)
        {

          mybit->j_var = j;
      
          for (i = 0; i < num_threads; i++)
          {
            rc = pthread_create(&my_thread[i], NULL, threadbuilder, (void *)mybit);
            if (rc)
            {
              printf("error\n");
              exit(-1);
            }
  
          }



          for (i = 0; i < num_threads; i++)
          {

            rc = pthread_join(my_thread[i],NULL);
            if (rc)
            {
              printf("error\n");
              exit(-1);
            }
          }
	}




}

// Void * function that the pthreads use to build each line of the bitmap image. Passes in the struct to have access to all the vars.

void * threadbuilder (void * mybit)
{
  int i;
  struct image_vars *mybit2 = mybit;

  for (i = 0; i < mybit2->width_var; i++)
  {

    // Determine the point in x,y space for that pixel.
    double x = mybit2->xmin_var + i*(mybit2->xmax_var-mybit2->xmin_var)/mybit2->width_var;
    double y = mybit2->ymin_var + mybit2->j_var*(mybit2->ymax_var-mybit2->ymin_var)/mybit2->height_var;

    // Compute the iterations at that point.
    int iters = iterations_at_point(x,y,mybit2->max_var);

    // Set the pixel in the bitmap.
    bitmap_set(mybit2->bm,i,mybit2->j_var,iters);
  }

  pthread_exit(NULL);

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

	while( (x*x + y*y <= 4) && iter < max ) {

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
