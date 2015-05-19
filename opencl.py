#!/usr/bin/python

import numpy as np
import pyopencl as cl

import global_vars as g

class OpenCL:
   def __init__(self):
      '''Initilize OpenCL.'''
      self.context=cl.create_some_context()
      self.queue=cl.CommandQueue(self.context)

   def load_program(self,filename):
      '''Loads OpenCL kernel from file.'''
      f = open(filename,'r')
      fstr = "".join(f.readlines())
      self.program=cl.Program(self.context,fstr).build()
      
   def sumint(self,qsize,ehc,x_pixels,y_pixels,points,sym=0,small=0):
      '''Returns the normalized intensity.'''
      npts=points.shape[0]
      buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(points))  #Copy input arrays into memory.
      out=np.empty(x_pixels*y_pixels,dtype=np.float32)
      out_buffer = cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)   #for the output

      if sym == 0:
         self.program.sumint00(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_points,np.int32(npts),out_buffer)
      elif small == 0:
         self.program.sumint10(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_points,np.int32(npts),out_buffer)
      else:
         self.program.sumint11(self.queue,out.shape,None,np.float32(qsize),np.int32(x_pixels),buf_points,np.int32(y_pixels),np.int32(npts),out_buffer)
      #cl.enqueue_read_buffer(self.queue,out_buffer,out).wait()
      cl.enqueue_copy(self.queue,out,out_buffer)
      return out.reshape(y_pixels,-1)    #converts to square    #TODO: or y_pixels??

   def sumint_mask(self,qsize,ehc,mask,points,sym=0,small=0):
      '''Returns the normalized intensity.'''
      npts=points.shape[0]
      x_pixels,y_pixels = mask.shape
      x_coords,y_coords=np.asarray([[i,j] for i in range(mask.shape[0]) for j in range(mask.shape[1]) if mask[i,j]]).T  #get x,y coordinates of non-zero mask points
      buf_x = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.int32(x_coords))
      buf_y = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.int32(y_coords))
      buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(points))  #Copy input arrays into memory.
      out=np.empty(len(x_coords),dtype=np.float32)
      out_buffer = cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)   #for the output

      sym = 0
      if sym == 0:
         self.program.sumint00mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
      elif small == 0:
         self.program.sumint10mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
      else:
         self.program.sumint11mask(self.queue,out.shape,None,np.float32(qsize),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
      #cl.enqueue_read_buffer(self.queue,out_buffer,out).wait()
      cl.enqueue_copy(self.queue,out,out_buffer)
      #return x_coords,y_coords,out
      out2d=np.zeros_like(mask)
      for i in range(len(x_coords)):
         out2d[x_coords[i],y_coords[i]] = out[i]
      return out2d

   def density(self,random_points):
      '''Returns list of densities.'''
      shape = g.dictionary_SI['shape']
      #buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(random_points))  #Copy input arrays into memory.
      buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(np.append(random_points,np.zeros([len(random_points),1]),1)))  #Copy input arrays into memory.
      out=np.empty(random_points.shape[0],dtype=np.float32)
      out_buffer = cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)   #for the output
      if shape == 1:
         radius_1,rho_1 = g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1']
         self.program.d1sphere(self.queue,out.shape,None,np.float32(radius_1),np.float32(rho_1),buf_points,out_buffer)
      elif shape == 2:
         radius_1,rho_1 = g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1']
         self.program.d2cylinder(self.queue,out.shape,None,np.float32(radius_1),np.float32(rho_1),buf_points,out_buffer)
      elif shape == 3:
         radius_1,rho_1 = g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1']
         radius_2,rho_2 = g.dictionary_SI['radius_2'],g.dictionary_SI['rho_2']
         self.program.d3coreshell(self.queue,out.shape,None,np.float32(radius_1),np.float32(rho_1),np.float32(radius_2),np.float32(rho_2),buf_points,out_buffer)
      elif shape == 4:
         radius_2 = g.dictionary_SI['radius_2']
         self.program.d4gaussian(self.queue,out.shape,None,np.float32(radius_2),buf_points,out_buffer)
      cl.enqueue_copy(self.queue,out,out_buffer)
      return out
         


#no mask; implement mask so x,y coordinates of relevant pixels are passed in.

if __name__ == "__main__":
   instance = OpenCL()
   instance.load_program('sumint.cl')
   intensity = instance.sumint(qsize,ehc,pixels,points)
