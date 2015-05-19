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
   const float radius_1, const float rho_1, const float radius_2, const float rho_2,
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
   density[n] = exp(-dist/pow(radius_2,2))
}

