\\This file contains the OpenCL code for much faster calculations of the intensity.
  Separate kernels exist both with and without symmetry and the small angle approximation, as well as those without any points masked.
  The first number in the title indicates if it's summetric (1=symmetric,0=asymmetric).
  The second number in the title indicates if the small-angle approximation is used.\\

#define PYOPENCL_DEFINE_CDOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64: enable

inline float my_dot (float* a, float*b, int length) {
   float c = 0;
   for (int i = 0; i < length; i++) {
      c += a[i] * b[i];
   }
   return c;
}

__kernel void sumint00 (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For no symmetry and no small angle approximation //
   float temp_intensity,temp_intensity_2,QdotR;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float Q[3] = { i*qsize/y_pixels-0.5*qsize, j*qsize/x_pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*y_pixels),2)+pow(j-0.5*x_pixels,2))*qsize/(x_pixels*2*ehc)),2) };
   temp_intensity = 0;
   temp_intensity_2 = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       QdotR = my_dot(Q,R,3);
       temp_intensity += points[p][3]*cos(QdotR);
       temp_intensity_2 += points[p][3]*sin(QdotR);
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2);
}

__kernel void sumint10 (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry but no small angle approximation //
   float temp_intensity;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float Q[3] = { i*qsize/y_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*y_pixels),2)+pow(j-0.5*x_pixels,2))*qsize/(x_pixels*2*ehc)),2) };
   temp_intensity = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,3));
   }
   intensity[n] = pow(temp_intensity,2);
}

__kernel void sumint11 (
   const float qsize, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry and small angle approximation //
   float temp_intensity = 0;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float Q[2] = { i*qsize/y_pixels-0.5*qsize, j*qsize/x_pixels-0.5*qsize};
   for ( int p = 0; p < npts; p++) {
       float R[2] = {points[p][0],points[p][1]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,2));
   }
   intensity[n] = pow(temp_intensity,2);
}

//Pass in separate x and y arrays.//
__kernel void sumint00mask (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels, __constant int* xval, __constant int* yval,
   __global const float4* points, const int npts, __global float* intensity)
{
// For no symmetry and no small angle approximation //
   float temp_intensity,temp_intensity_2,QdotR;
   int n = get_global_id(0);
   int i = xval[n];
   int j = yval[n];
   float Q[3] = { i*qsize/y_pixels-0.5*qsize, j*qsize/x_pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*y_pixels),2)+pow(j-0.5*x_pixels,2))*qsize/(x_pixels*2*ehc)),2) };
   temp_intensity = 0;
   temp_intensity_2 = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       QdotR = my_dot(Q,R,3);
       temp_intensity += points[p][3]*cos(QdotR);
       temp_intensity_2 += points[p][3]*sin(QdotR);
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2);
}


__kernel void sumint10mask (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels, __constant int* xval, __constant int* yval,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry but no small angle approximation //
   float temp_intensity;
   int n = get_global_id(0);
   int i = xval[n];
   int j = yval[n];
   float Q[3] = { i*qsize/y_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*y_pixels),2)+pow(j-0.5*x_pixels,2))*qsize/(x_pixels*2*ehc)),2) };
   temp_intensity = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,3));
   }
   intensity[n] = pow(temp_intensity,2);
}

__kernel void sumint11mask (
   const float qsize, const int x_pixels, const int y_pixels, __constant int* xval, __constant int* yval,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry and small angle approximation //
   float temp_intensity = 0;
   int n = get_global_id(0);
   int i = xval[n];
   int j = yval[n];
   float Q[2] = { i*qsize/y_pixels-0.5*qsize, j*qsize/x_pixels-0.5*qsize};
   for ( int p = 0; p < npts; p++) {
       float R[2] = {points[p][0],points[p][1]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,2));
   }
   intensity[n] = pow(temp_intensity,2);
}
