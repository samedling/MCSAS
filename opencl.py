#!/usr/bin/python

# This file contains the basic OpenCL initilization commands in a class OpenCL.
# After importing OpenCL from this file, create and instance of the OpenCL class, and then load a program.
# instance = OpenCL()
# instance.load_program(filename)
# It also contains functions to perform various calculations:
# sumint: calculates the intensity at all points
# sumint_mask: calculates the intensity at unmasked points only
# density: calculates density for some of the possible shapes (due to if statements, OpenCL isn't really ideal for this)


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
      
   def sumint(self,qsize,ehc,x_pixels,y_pixels,points,sym=0,small=0,coherence_length=0):
      '''Returns the normalized intensity.'''
      npts=points.shape[0]
      buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(points))  #Copy input arrays into memory.
      out=np.empty(x_pixels*y_pixels,dtype=np.float32)
      out_buffer = cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)   #for the output

      if coherence_length and (coherence_length < g.dictionary_SI['z_dim']):
         self.program.sumint_long(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.float32(coherence_length),np.int32(x_pixels),np.int32(y_pixels),buf_points,np.int32(npts),out_buffer)
      elif sym == 0:
         if small == 0:
            self.program.sumint_asym(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_points,np.int32(npts),out_buffer)
         else:
            self.program.sumint_asym_small(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_points,np.int32(npts),out_buffer)
      else:
         if small == 0:
            self.program.sumint_sym(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_points,np.int32(npts),out_buffer)
         else:
            self.program.sumint_sym_small(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_points,np.int32(npts),out_buffer)
      #cl.enqueue_read_buffer(self.queue,out_buffer,out).wait()
      cl.enqueue_copy(self.queue,out,out_buffer)
      return out.reshape(y_pixels,-1)    #Converts 1D intensity back to 2D.

   def sumint_mask(self,qsize,ehc,mask,points,sym=0,small=0,coherence_length=0):
      '''Returns the normalized intensity.'''
      npts=points.shape[0]
      x_pixels,y_pixels = mask.shape
      x_coords,y_coords=np.asarray([[i,j] for i in range(mask.shape[0]) for j in range(mask.shape[1]) if mask[i,j]]).T  #get x,y coordinates of non-zero mask points
      buf_x = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.int32(x_coords))
      buf_y = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.int32(y_coords))
      buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(points))  #Copy input arrays into memory.
      out=np.empty(len(x_coords),dtype=np.float32)
      out_buffer = cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)   #for the output

      try:
         if coherence_length and coherence_length > g.dictionary_SI['z_dim']:
            self.program.sumint_long_mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.float32(coherence_length),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
         elif sym == 0:
            if small == 0:
               g.dprint("Not using radial symmetry or small angle approximation.")
               self.program.sumint_asym_mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
            else:
               g.dprint("Using small angle approximation but not radial symmetry.")
               self.program.sumint_asym_small_mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
         else:
            if small == 0:
               g.dprint("Using radial symmetry but not small angle approximation.")
               self.program.sumint_sym_mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
            else:
               g.dprint("Using both radial symmetry and small angle approximation.")
               self.program.sumint_sym_small_mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
      except:
         print("Runtime Error. Too many pixels being calculated for OpenCL to work with grid compression.")
         print("Either disable grid compression entirely or decrease number of pixels or increase grid compression.")
         if g.debug:
            print("Number of pixels: {0}.".format(len(x_coords)))
            print("(Number of points: {0}.)".format(npts))
         if g.dictionary_SI['grid_compression'] < 5:
            print("Try grid compression of 5 or 10")
         self.program.sumint_asym_mask(self.queue,out.shape,None,np.float32(qsize),np.float32(ehc),np.int32(x_pixels),np.int32(y_pixels),buf_x,buf_y,buf_points,np.int32(npts),out_buffer)
         return
      #cl.enqueue_read_buffer(self.queue,out_buffer,out).wait()
      cl.enqueue_copy(self.queue,out,out_buffer)
      #return x_coords,y_coords,out
      out2d=np.zeros_like(mask)
      for i in range(len(x_coords)):    #Converts 1D intensity back to 2D.
         out2d[x_coords[i],y_coords[i]] = out[i]
      return out2d

   def density(self,random_points):
      '''Returns list of densities.'''
      shape = g.dictionary_SI['shape']
      #buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(random_points))  #Copy input arrays into memory.
      buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(np.append(random_points,np.zeros([len(random_points),1]),1)))  #Copy input arrays into memory. Pad with zeros to use float4*.
      out=np.empty(random_points.shape[0],dtype=np.float32)
      out_buffer = cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)   #for the output
      #TODO: consider all possible variables as a single input array instead of 5 separate variables
      if shape == 1:
         radius_1,rho_1 = g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1']
         self.program.d1sphere(self.queue,out.shape,None,np.float32(radius_1),np.float32(rho_1),buf_points,out_buffer)
      elif shape == 2:
         radius_1,rho_1 = g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1']
         self.program.d2cylinder(self.queue,out.shape,None,np.float32(radius_1),np.float32(rho_1),buf_points,out_buffer)
      elif shape == 3:
         radius_1,radius_2,rho_1,rho_2 = g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['rho_2']
         self.program.d3coreshell(self.queue,out.shape,None,np.float32(radius_1),np.float32(radius_2),np.float32(rho_1),np.float32(rho_2),buf_points,out_buffer)
      elif shape == 4:
         radius_2 = g.dictionary_SI['radius_2']
         self.program.d4gaussian(self.queue,out.shape,None,np.float32(radius_2),buf_points,out_buffer)
      elif shape == 5:
         radius_1,radius_2,rho_1,z_dim = g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['z_dim']
         self.program.d5choppedcone(self.queue,out.shape,None,np.float32(radius_1),np.float32(radius_2),np.float32(rho_1),np.float32(z_dim),buf_points,out_buffer)
      elif shape == 6:
         radius_1,rho_1 = g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1']
         self.program.d6hexprism(self.queue,out.shape,None,np.float32(radius_1),np.float32(rho_1),buf_points,out_buffer)
      elif shape == 7:
         radius_2,rho_1 = g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1']
         self.program.d6hexprism(self.queue,out.shape,None,np.float32(radius_2),np.float32(rho_1),buf_points,out_buffer)
      elif shape == 11:
         radius_1,radius_2,rho_1 = g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1']
         self.program.d11doubleslit(self.queue,out.shape,None,np.float32(radius_1),np.float32(radius_2),np.float32(rho_1),buf_points,out_buffer)
      elif shape == 13:
         radius_1,radius_2,rho_1,rho_2,z_dim = g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary['rho_2'],g.dictionary_SI['z_dim']
         self.program.d13sine(self.queue,out.shape,None,np.float32(radius_1),np.float32(radius_2),np.float32(rho_1),np.float32(rho_2),np.float32(z_dim),buf_points,out_buffer)
      elif shape == 14:
         radius_1,radius_2,rho_1,z_dim = g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['z_dim']
         self.program.d14doublecone(self.queue,out.shape,None,np.float32(radius_1),np.float32(radius_2),np.float32(rho_1),np.float32(z_dim),buf_points,out_buffer)
      cl.enqueue_copy(self.queue,out,out_buffer)
      return out
         

