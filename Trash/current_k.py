import numpy as np;
from mpl_toolkits.mplot3d import axes3d;
import matplotlib.pyplot as plt;
from matplotlib import cm
import argparse as argparse


parser = argparse.ArgumentParser(prog="python main.py",
  description = "Takes Delta(k) and plots phase/amp graphs. Use -help to see more information."); 
parser.add_argument('-m', '--mode', help='1 s, 2 g, 3 x2-y2, 4 xy.', default = 1, action='store', type = int);
parser.add_argument('-p1', '--paramfirst', help='First parameter for plotting function.', default = 1.0, action='store', type = float);
parser.add_argument('-p2', '--paramsecond', help='Second parameter for plotting function.', default = 1.0, action='store', type = float);
parser.add_argument('-p3', '--paramthird', help='Third parameter for plotting function.', default = 1.0, action='store', type = float);
parser.add_argument('-ek', '--electronlevel', help='Electron level for tunnel, or something.', default = 1.0, action='store', type = float);
parser.add_argument('-kx', '--xoffset', help='Left/right k_x difference.', default = 1.0, action='store', type = float);
parser.add_argument('-ky', '--yoffset', help='Left/right k_y difference.', default = 1.0, action='store', type = float);
parser.add_argument('-n', '--numbercells', help='Number of unit cells in the BZ to plot.', default = 1, action='store', type = int);
parser.add_argument('-r', '--resolution', help='Resolution is the number of points along an axis.', default = 25, action='store', type = int);
args = parser.parse_args();

numberOfCells = args.numbercells;
resolution = numberOfCells*args.resolution;

k_x, k_y = np.meshgrid(
   np.linspace(-numberOfCells*2*np.pi, numberOfCells*2*np.pi, resolution),
   np.linspace(-numberOfCells*2*np.pi, numberOfCells*2*np.pi, resolution)
);
amplitudeL = 0;
phaseR = 0;


amplitudeR = 0;
phaseL = 0;

offsetX = args.xoffset;
offsetY = args.yoffset;
if args.mode == 1:

  print "Amplitude/Phase for Delta_s (mode %s). Note that delta_s^2 = 0, because we aren't bothered with k_z. " % args.mode;  
  delta_0 = args.paramfirst;
  delta_1 = args.paramsecond;
  
  
  delta = delta_0 + delta_1 * ( np.cos(k_x) + np.cos(k_y) ); 
  amplitudeL = np.absolute(delta);
  phaseL = np.angle(delta); 
  
  delta = delta_0 + delta_1 * ( np.cos(k_x+offsetX) + np.cos(k_y+offsetY) ); 
  amplitudeR = np.absolute(delta);
  phaseR = np.angle(delta); 
  
  
elif args.mode == 2:
  print "Amplitude/Phase for Delta_g (mode %s) " % args.mode;

  delta_0 = args.paramfirst;
  delta_1 = args.paramsecond;
  delta_2 = args.paramthird;
  
  delta = delta_0*(np.sin(2*k_x)*np.sin(k_y) - np.sin( 2*k_y)*np.sin(k_x));
   
  amplitudeL = np.absolute(delta);
  phaseL = np.angle(delta); 
  
  delta = delta_0*(np.sin(2*(k_x+offsetX))*np.sin(k_y+offsetY) - np.sin(2*(k_y+offsetY))*np.sin(k_x+offsetX)); 
  amplitudeR = np.absolute(delta);
  phaseR = np.angle(delta); 
elif args.mode == 3:
  print "Amplitude/Phase for Delta_{x2_y2} (mode %s) " % args.mode;

  delta_0 = args.paramfirst;
  delta_1 = args.paramsecond;
  delta_2 = args.paramthird;
  
  delta = delta_0 * (np.cos(k_x) - np.cos(k_y))
  
  amplitudeL = np.absolute(delta);
  phaseL = np.angle(delta); 
  
  delta = delta_0 * (np.cos(k_x+offsetX) - np.cos(k_y+offsetY))
  amplitudeR = np.absolute(delta);
  phaseR = np.angle(delta); 

elif args.mode == 4:
  print "Amplitude/Phase for Delta_{xy} (mode %s) " % args.mode;

  delta_0 = args.paramfirst;
  delta_1 = args.paramsecond;
  delta_2 = args.paramthird;
  
  delta = delta_0 * np.sin(k_x) * np.sin(k_y)
  
  amplitudeL = np.absolute(delta);
  phaseL = np.angle(delta); 
  
  delta = delta_0 * np.sin(k_x+offsetX) * np.sin(k_y+offsetY)
  amplitudeR = np.absolute(delta);
  phaseR = np.angle(delta); 

electronLevel = args.electronlevel;
phiOffset = args.paramthird;

current = amplitudeL * amplitudeR/ ( np.sqrt(electronLevel**2 + amplitudeL**2) * np.sqrt(electronLevel**2 + amplitudeR**2) * ( np.sqrt(electronLevel**2 + amplitudeL**2) + np.sqrt(electronLevel**2 + amplitudeR**2)))


distance = 1.4
limit = 2.0

  
fig = plt.figure();
ax = fig.add_subplot(111,projection='3d');
 
ax.plot_surface(k_x, k_y, current, rstride=8, cstride=8, alpha=0.3)
cset = ax.contour(k_x, k_y, current, zdir='z', offset=-distance*np.abs(np.min(current)), cmap=cm.coolwarm) 
cset = ax.contour(k_x, k_y, current, zdir='x', offset=-distance*np.abs(np.min(k_x)), cmap=cm.coolwarm) 
cset = ax.contour(k_x, k_y, current, zdir='y', offset=distance*np.max(k_y), cmap=cm.coolwarm)
 

ax.set_xlabel('k_x')
ax.set_xlim(-limit*np.abs(np.min(k_x)), limit*np.max(k_x))
ax.set_ylabel('k_y')
ax.set_ylim(-limit*np.abs(np.min(k_y)), limit*np.max(k_y))
ax.set_zlabel('Ampl')
ax.set_zlim(-limit*np.abs(np.min(current)), limit*np.max(current))

plt.title("Current_k"); 
plt.show()
