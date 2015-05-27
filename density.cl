#define PYOPENCL_DEFINE_CDOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64: enable

inline float my_dot (float4 a, float4 b, int length) {
    float c = 0;
    for (int i = 0; i < length; i++) {
        c += a[i] * b[i];
    }
    return c;
}


__kernel void d1sphere (
   const float radius_1, const float rho_1,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   density[n] = 0;
   if (my_dot(points[n],points[n],3) < pow(radius_1,2)) {
      density[n] = rho_1;
   }
}

__kernel void d2cylinder (
   const float radius_1, const float rho_1,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   density[n] = 0;
   if (my_dot(points[n],points[n],2) < pow(radius_1,2)) {
      density[n] = rho_1;
   }
}


__kernel void d3coreshell (
   const float radius_1, const float radius_2, const float rho_1, const float rho_2,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   density[n] = 0;
   float dist = pow(points[n][0],2)+pow(points[n][1],2);
   if (dist < pow(radius_2,2)) {
      density[n] = rho_1;
   }
   else if (dist < pow(radius_1,2)) {
      density[n] = rho_2;
   }
}


__kernel void d4gaussian (
   const float radius_2,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   float dist = pow(points[n][0],2)+pow(points[n][1],2);
   density[n] = exp(-dist/pow(radius_2,2));
}


__kernel void d5choppedcone (
   const float radius_1, const float radius_2, const float rho_1 const float z_dim,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   density[n] = 0;
   float dist = pow(points[n][0],2)+pow(points[n][1],2);
   if (sqrt(dist) < points[n][2]*(radius_2-radius_1)/z_dim+(radius_1+radius_2)*0.5) {
      density[n] = rho_1;
   }
}

__kernel void d6hexprism (
   const float radius_1, const float rho_1,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   density[n] = rho_1;
   float sqrt3over2 = sqrt(3.0)*0.5;
   float coords = points[0:3]/radius_1 //???//
   if (pow(coords[1],2) > 0.75) or (coords[1]+(coords[0]-1)*sqrt3over2 > 0) or (coords[1]+(coords[0]+1)*sqrt3over2 < 0) or (coords[1]-(coords[0]-1)*sqrt3over < 0) or (coords[1]-(coords[0]+1)*sqrt3over2 > 0) {  //???//
      density[n] = 0;
   }
}

__kernel void d7rectprism (
   const float radius_2, const float rho_1,
   __global const float4* points, __global float* density)  //TODO: only uses column 1 of points, so make only pass in column 1.//
{
   int n = get_global_id(0);
   density[n] = 0;
   if (points[n][0] < radius_2) {
      density[n] = rho_1;
   }
}

__kernel void d11doubleslit (
   const float radius_1, const float radius_2, const float rho_1,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   density[n] = 0;
   if (((-radius_1*0.5 < points[n][0]) and (points[n][0] < -radius_2*0.5)) or ((radius_2*0.5 < points[n][0]) and (points[n][0] < radius_1*0.5))) {   //TODO: only uses column 1 of points, so make only pass in column 1.//
      density[n] = rho_1;
   }
}

__kernel void d13sine (
   const float radius_1, const float radius_2, const float rho_1, const float rho_2, const float z_dim,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   float dist = pow(points[n][0],2)+pow(points[n][1],2);
   density[n] = 0;
   if (sqrt(dist) < (radius_1+radius_2)*0.5 + (radius_1-radius_2)*sin(points[n][2]*rho_2*2*pi/z_dim)*0.5) {  //TODO: Define pi?? sin??//
      density[n] = rho_1;
   }
}

__kernel void d14doublecone (
   const float radius_1, const float radius_2, const float rho_1, const float z_dim,
   __global const float4* points, __global float* density)
{
   int n = get_global_id(0);
   float dist = pow(points[n][0],2)+pow(points[n][1],2);
   density[n] = 0;
   if (sqrt(dist) < radius_2+abs(points[n][2])*(radius_1-radius_2)/z_dim*0.5)) {    //TODO: abs??//
      density[n] = rho_1;
   }
}



