import os,sys,pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

import global_vars as g

###################               Plotting the points created              ##############################
def Points_Plot(Points, save_name, show=1):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.elev = g.dictionary_SI['altitude']
    ax.azim = g.dictionary_SI['azimuth']

    X = Points[:,0]
    Y = Points[:,1]
    Z = Points[:,2]
    vmin = np.min(Points[:,3])
    vmax = np.max(Points[:3])
    cat = ax.scatter(X, Y, Z, c=Points[:,3], s=50, vmin=vmin, vmax=vmax)    

    #When plotting all points, do you want the axes to have the same scale? yes = 1, no = 0.
    #This sets the default, if 01input_variables has not been run.

    if g.dictionary_SI['scale'] == 1:
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
        mean_x = X.mean()
        mean_y = Y.mean()
        mean_z = Z.mean()
        ax.set_xlim(mean_x - max_range, mean_x + max_range)
        ax.set_ylim(mean_y - max_range, mean_y + max_range)
        ax.set_zlim(mean_z - max_range, mean_z + max_range)
    ax.set_title(g.dictionary_SI['title'])


    plt.savefig(g.dictionary_SI['path_to_subfolder']+save_name+".png")
    plt.show()






##################           Plotting the Intensity                ########################

#save_name and title must be strings "title"
def Intensity_plot(Intensity, name, title, show):
    maximum = g.dictionary_SI['maximum']
    minimum = g.dictionary_SI['minimum']
    log_scale = g.dictionary_SI['log_scale']
    bound = g.dictionary_SI['bound']
    QSize = g.dictionary_SI['QSize']
    save_name = g.dictionary['save_name']
    ThreeD = g.dictionary['ThreeD']
   #getting the dimentions of the Intensity array
    pixels_row, pixels_col = np.array(Intensity).shape

    #This creates a new array which is bounded.
    if bound == 1:
        def limit(x):
            return max(min(x,maximum), minimum)
        limit = np.vectorize(limit)
        bound_intensity = limit(Intensity)

    if 0:#A different way you could plot the intensity - 8bit bmp
        from PIL import Image
        im=Image.Image()
        #im=Image.fromarray(np.log10(Intensity))
        im=Image.fromarray(np.uint32(128*np.log10(Intensity)/np.log10(np.amax(Intensity))))
        im.save(g.dictionary_SI['path_to_subfolder']+'Intensity.bmp')
    

    if ThreeD == 0: #For a 2D heat map
        fig = plt.figure()
        #ax = fig.add_subplot(111)
        ax = fig.add_subplot(111, aspect='equal')
        #CHANGE COLOUR MAP HERE!!
        #Pictures of the colour map are at: http://matplotlib.org/users/colormaps.html
        #put _r at the end to reverse the direction of the colour.
        #the ones below have been good.
        
        #cmap = cm.get_cmap('gist_gray_r')
        #cmap = cm.get_cmap('seismic_r')
        cmap = cm.get_cmap('hot_r')

        if bound == 1:
            if log_scale == 1: #bound with log scale
                cs = plt.pcolor(np.log10(bound_intensity), cmap = cmap)
            else: #bound, no log scale
                cs = plt.pcolor(bound_intensity, cmap = cmap)
        else:
            if log_scale == 1: #not bound with log scale
                cs = plt.pcolor(np.log10(Intensity), cmap = cmap)
            else: #not bound, no log scale
                cs = plt.pcolor(Intensity, cmap = cmap)
    
        #This plots the colour bar.        
        cb = plt.colorbar(cs, orientation='vertical')

        #plt.xticks(range(0, pixels_row, 25), np.arange(-QSize/2,QSize/2, QSize/5), rotation='vertical')
        #plt.yticks(range(0, pixels_col, 10), np.arange(-QSize/2,QSize/2, QSize/11))

        
    else: #For a 3D plot
    
        #CHANGE THE SCALE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        X, Y =np.meshgrid(np.arange(-QSize/2,QSize/2, QSize/pixels_row), np.arange(-QSize/2,QSize/2, QSize/pixels_col))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        #You can change the angle of elevation and azimith - think of an alt/az  telescope
        #leave commented to get a different projection.
        #ax.elev = 0 #to view side on.
        ax.elev = g.dictionary_SI['altitude']
        ax.azim = g.dictionary_SI['azimuth']
    
        #make 0 to cut off the top and bottom.
        if bound == 1:
            if log_scale == 1: #bound with log scale
                p = ax.plot_surface(X, Y, np.log10(bound_intensity), rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            else: #bound, no log scale
                p = ax.plot_surface(X, Y,bound_intensity, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
                
        else:
            if log_scale == 1: #not bound with log scale
                p = ax.plot_surface(X, Y, np.log10(Intensity), rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            else: #not bound, no log scale
                p = ax.plot_surface(X, Y, Intensity, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        #to change the size of the colourbar. I don't know how to change colours for the 3D plot.
        cb = fig.colorbar(p, shrink=0.5)


    ax.set_title(title)
    plt.savefig(g.dictionary_SI['path_to_subfolder']+name+".png")   #should name be save_name?
    if show:
        plt.show()


def Fit_plot(experimental,fit,residuals,orientation=False):
    '''Shows a multiplot containing residuals, experimental data, and fit results.'''
    maximum = g.dictionary_SI['maximum']
    minimum = g.dictionary_SI['minimum']
    log_scale = g.dictionary_SI['log_scale']
    bound = g.dictionary_SI['bound']
    QSize = g.dictionary_SI['QSize']
    save_name = g.dictionary['save_name']
   #getting the dimentions of the Intensity array
    pixels_row, pixels_col = np.array(residuals).shape

    if not orientation:
        if pixels_row > pixels_col:
            orientation = 'Horizontal'
        else:
            orientation = 'Vertical'

    #This creates new arrays which are bounded.
    if bound:
        #def limit(x):
        #    return max(min(x,maximum), minimum)
        #limit = np.vectorize(limit)
        limit = np.vectorize(lambda x: max(min(x,maximum), minimum))
        exp_plot = limit(experimental)
        fit_plot = limit(fit)
        limit = np.vectorize(lambda x: max(min(x,maximum), minimum/10)) #Using smaller minimum for residuals.
        res_plot = limit(residuals)
    else:
        exp_plot = experimental
        fit_plot = fit
        res_plot = residuals

    if orientation[0] in ('V','v'):    #If vertical
       fig,axes = plt.subplots(nrows=3,ncols=1,figsize=(5,10))
    else:
       #fig,axes = plt.subplots(nrows=1,ncols=3,figsize=(12,3))
       fig,axes = plt.subplots(nrows=1,ncols=3,figsize=(8,3))

    #Pictures of the colour map are at: http://matplotlib.org/users/colormaps.html
    #put _r at the end to reverse the direction of the colour.
    #cmap = cm.get_cmap('gist_gray_r')
    #cmap = cm.get_cmap('seismic_r')
    cmap = cm.get_cmap('hot_r')

    if log_scale == 1:
       im = axes[0].imshow(np.log10(exp_plot),cmap=cmap)
       im = axes[1].imshow(np.log10(fit_plot),cmap=cmap)
       im = axes[2].imshow(np.log10(res_plot),cmap=cmap)
    else:
       im = axes[0].imshow(exp_plot,cmap=cmap)
       im = axes[1].imshow(fit_plot,cmap=cmap)
       im = axes[2].imshow(res_plot,cmap=cmap)

    axes[0].set_title('Exp. Data')
    axes[1].set_title('Calc. Intensity')
    axes[2].set_title('|Difference|')

    for i in range(3):
        axes[i].set_xticklabels([])
        axes[i].set_yticklabels([])

    if orientation[0] in ('V','v'):    #If vertical
       fig.subplots_adjust(right=0.85)  #Moves subplots to make room for color bar.
       cbar_ax = fig.add_axes([0.80,0.15,0.05,0.7])
    else:
       #fig.subplots_adjust(right=0.92)  #Moves subplots to make room for color bar.
       #cbar_ax = fig.add_axes([0.95,0.15,0.02,0.7])
       fig.subplots_adjust(right=0.82)  #Moves subplots to make room for color bar.
       cbar_ax = fig.add_axes([0.90,0.15,0.02,0.7])
    fig.colorbar(im, cax=cbar_ax)        

    #plt.savefig(g.dictionary_SI['path_to_subfolder']+save_name+".png")
    plt.show()





##################              Plotting the Radial Intensity          #########################
#this defines the analytic model as a function, if possible.
def radial_intensity_plot(radial_intensity, name, title, show):
    analytic = g.dictionary_SI['analytic']
    save_name = g.dictionary_SI['save_name']
    maximum = g.dictionary_SI['maximum']
    minimum = g.dictionary_SI['minimum']
    plt.close()
    plt.cla()
    plt.clf()

    if g.dictionary_SI['bound'] ==1:
        normalised = [max(min(x[1],maximum), minimum) for x in radial_intensity]
    else:
        normalised=radial_intensity[:,1]

    if g.dictionary_SI['log_scale'] == 1:
        plt.loglog(radial_intensity[:,0]*10**-9, normalised,'.b')#this plots points
        plt.loglog(radial_intensity[:,0]*10**-9, normalised,'-b')#this plots lines between the points
    else:
        plt.plot(radial_intensity[:,0]*10**-9, normalised,'.b')#this plots points
        plt.plot(radial_intensity[:,0]*10**-9, normalised,'-b')#this plots lines between the points
        
    #this plots the analytic model, if there is one
    if analytic!=0:
        theoretical = theory(radial_intensity[:,0])
        ratio = np.amax(radial_intensity[:,1])/np.amax(theoretical)
        TR = theoretical*ratio
        if g.dictionary_SI['bound'] ==1:
             TR = [max(min(x,maximum), minimum) for x in TR]
        if g.dictionary_SI['log_scale'] == 1:
            plt.loglog(radial_intensity[:,0]*10**-9, TR,'r-')
        else:
            plt.plot(radial_intensity[:,0]*10**-9, TR,'r-')

    x1, x2, y1, y2 = plt.axis()
    plt.axis((x1, x2, y1, 1.5*y2))

    pylab.title(title)
    pylab.xlabel('Q Range (nm^-1)')
    pylab.ylabel('Relative Intensity')
        
    plt.savefig(g.dictionary_SI['path_to_subfolder']+name+".png")
    if show == 1:
        pylab.show()
    


##################           plotting by angle            #################

def angle_plot(data, name, title, show):
    analytic = g.dictionary_SI['analytic']
    save_name = g.dictionary_SI['save_name']
    maximum = g.dictionary_SI['maximum']
    minimum = g.dictionary_SI['minimum']
    plt.close()
    plt.cla()
    plt.clf()

    if g.dictionary_SI['bound'] ==1:
        normalised = [max(min(x[1],maximum), minimum) for x in data]
    else:
        normalised=data[:,1]

    if g.dictionary_SI['log_scale'] == 1:
        plt.semilogy(data[:,0], normalised,'.b')#this plots points
        plt.semilogy(data[:,0], normalised,'-b')#this plots lines between the points
    else:
        plt.plot(data[:,0], normalised,'.b')#this plots points
        plt.plot(data[:,0], normalised,'-b')#this plots lines between the points
    x1, x2, y1, y2 = plt.axis()
    plt.axis((x1, x2, y1, 1.5*y2))
#plt.axis((0, 6.29, y1, 1.5*y2))

    pylab.title(title)
    pylab.xlabel('Angle (radians)')
    pylab.ylabel('Relative Intensity')

    
    plt.savefig(g.dictionary_SI['path_to_subfolder']+name+".png")
    if show == 1:
        pylab.show()



