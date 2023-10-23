# Creating animation
import sys
sys.path.insert(0, '/home/sid/athena/vis/python')
import athena_read
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from glob import glob
import natsort
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.gridspec import GridSpec

# N = 128
file_num = 2

files = glob("/home/sid/athena/bin/Blast.out1.*.athdf")
files = natsort.natsorted(files)

data = athena_read.athdf(files[file_num],num_ghost=2)
print(np.shape(data['vel2']))

print(data)

v0 = 3e+6  # cm/s

# np.save("x_arr_2048",data['x1v'])

# r = data['x1v']
# theta = data['x2v']
# phi = data['x3v']

# den = data['rho'][0,0,:]

# pres = 0.011*den

v_r_rms_wedge = []
time_wedge = []
for i in range(len(files)):
    data = athena_read.athdf(files[i],num_ghost=2)
    time_wedge.append(data['Time'])
    v_r = data['vel1'][0,0,:]
    v_r2 = v_r**2
    v_r_rms_wedge.append(np.sqrt(np.mean(v_r2)-np.mean(v_r)**(2)))
    # pres = 0.011*den
    # plt.show()
v_r_rms_wedge = np.array(v_r_rms_wedge)*v0
time_wedge = np.array(time_wedge)
plt.plot(time_wedge,v_r_rms_wedge)
# plt.ylim([1e-2,1e5])
plt.ylabel(r"$v_{r_{rms}}$")
plt.xlabel("t")
plt.yscale("log")
# plt.legend()
plt.show() 

# gs = GridSpec(1, 1, wspace=0.7)
# fig = plt.figure()
# ax1 = fig.add_subplot(gs[0,0])
# # ax2 = fig.add_subplot(gs[0,1])

# def animate(i):
#     data = athena_read.athdf(files[i],num_ghost=2)
#     r = data['x1v']
#     den = data['rho'][0,0,:]
#     time = data['Time']
#     # pres = data['press'][0,0,:]
#     ax1.plot(r,den)
#     ax1.set_xlabel("r")
#     ax1.set_ylabel(r"$\rho$")
#     ax1.set_yscale("log")
#     # ax2.plot(r,0.011*den)
#     # ax2.set_xlabel("r")
#     # ax2.set_ylabel("P")
#     fig.suptitle('Time = '+str(time))
#     print(i/len(files)*100)

# anim = animation.FuncAnimation(fig, animate, frames=len(files), interval=100, blit=True) #init_func=init
# anim.save('/home/sid/athena/vis/python/iso_approx.mp4', fps=1)

# j = 1

# files = glob("/home/sid/athena/isothermal_hotjupiter/iso_N="+str(2**(j+5))+"/1D/Blast.out1.*.athdf")
# files = natsort.natsorted(files)
# v0 = 3e+6  # cm/s

# fig = plt.figure()
# vr_rms = []
# t = []
# for i in range(len(files)):
#     data = athena_read.athdf(files[i],num_ghost=2)
#     vr = data['vel1'][0,0,:]
#     time = data['Time']
#     vr_2 = vr**2
#     vr_rms.append(np.sqrt(np.mean(vr_2)-np.mean(vr)**2))
#     t.append(time)
# vr_rms = np.array(vr_rms)*v0
# t = np.array(t)
# plt.plot(t,vr_rms)
# plt.xlabel('time')
# plt.ylabel(r'$v_{r_{rms}}$')
# plt.yscale('log')
# plt.title('Resolution = '+str(2**(j+5)))
# plt.savefig("/home/sid/athena/isothermal_hotjupiter/iso_vr_N="+str(2**(j+5))+"_1D.png")
# plt.show()

# fig = plt.figure()
# vtheta_rms = []
# t = []
# for i in range(len(files)):
#     data = athena_read.athdf(files[i],num_ghost=2)
#     vtheta = data['vel2'][0,0,:]
#     time = data['Time']
#     vtheta_2 = vtheta**2
#     vtheta_rms.append(np.sqrt(np.mean(vtheta_2)-np.mean(vtheta)**2))
#     t.append(time)
# vtheta_rms = np.array(vtheta_rms)*v0
# t = np.array(t)
# plt.plot(t,vtheta_rms)
# plt.xlabel('time')
# plt.ylabel(r'$v_{(\theta)_{rms}}$')
# plt.yscale('log')
# plt.title('Resolution = '+str(2**(j+7)))
# plt.savefig("/home/sid/athena/isothermal_hotjupiter/iso_vtheta_N="+str(2**(j+7))+".png")
# # plt.show()

# fig = plt.figure()
# vphi_rms = []
# t = []
# for i in range(len(files)):
#     data = athena_read.athdf(files[i],num_ghost=2)
#     vphi = data['vel3'][0,0,:]
#     time = data['Time']
#     vphi_2 = vphi**2
#     vphi_rms.append(np.sqrt(np.mean(vphi_2)-np.mean(vphi)**2))
#     t.append(time)
# vphi_rms = np.array(vphi_rms)*v0
# t = np.array(t)
# plt.plot(t,vphi_rms)
# plt.xlabel('time')
# plt.ylabel(r'$v_{(\phi)_{rms}}$')
# # plt.yscale('log')
# plt.title('Resolution = '+str(2**(j+7)))
# plt.savefig("/home/sid/athena/isothermal_hotjupiter/iso_vphi_N="+str(2**(j+7))+".png")
# # plt.show()


# files = glob("/home/sid/athena/bin/disk.out1.*.athdf")
# files = natsort.natsorted(files)

# data = athena_read.athdf(files[0])

# r = data['x1v']
# theta = data['x2v']
# phi = data['x3v']

# den = data['rho']

# w = np.load('w_array.npy')

# data = athena_read.athdf(files[0])
# x = data['x1v']
# y = data['x2v']
# z = w[:,:,0] #np.transpose(data['rho'][0])
# mycmap = plt.get_cmap('viridis') #plt.get_cmap('gist_heat')

# gs = GridSpec(2, 4, width_ratios=[1,0.05,1,0.05], hspace=1, wspace=0.7)
# fig = plt.figure()
# ax1 = fig.add_subplot(gs[0,0], polar=True)
# cbar_ax1 = fig.add_subplot(gs[0,1])
# ax2 = fig.add_subplot(gs[1,0], polar=True)
# cbar_ax2 = fig.add_subplot(gs[1,1])
# ax3 = fig.add_subplot(gs[0,2], polar=True)
# cbar_ax3 = fig.add_subplot(gs[0,3])
# ax4 = fig.add_subplot(gs[1,2], polar=True)
# cbar_ax4 = fig.add_subplot(gs[1,3])

# z = np.transpose(data['rho'][0])
# cs1 = ax1.contourf(y, x, z, cmap=mycmap)
# z = np.transpose(data['vel1'][0])
# cs2 = ax2.contourf(y, x, z, cmap=mycmap)
# z = np.transpose(data['vel2'][0])
# cs3 = ax3.contourf(y, x, z, cmap=mycmap)
# z = w[:,:,0]
# cs4 = ax4.contourf(y, x, z, cmap=mycmap)

# fig.colorbar(cs1,cax=cbar_ax1)
# fig.colorbar(cs2,cax=cbar_ax2)
# fig.colorbar(cs3,cax=cbar_ax3)
# fig.colorbar(cs4,cax=cbar_ax4)

# def animate(i):
#     data = athena_read.athdf(files[i])
#     r = data['x1f'][0:N]
#     theta = data['x2f'][0:64]
#     rho = data['rho'][0,:,:]
#     v_r = data['vel1'][0,:,:]
#     v_theta = data['vel2'][0, :, :]
#     mycmap = plt.get_cmap('viridis')

#     # row, col = np.where(rho>0.001)
#     # col = np.unique(col)
#     # shape = len(col)
#     # if shape != 0:
#     #     cs1 = ax1.contourf(theta, r[col], w[col,:,i], cmap=mycmap)  # np.log10(rho.T)
#     #     fig.colorbar(cs1,cax=cbar_ax1)
#     cs1 = ax1.contourf(theta, r, rho.T, cmap=mycmap)  # np.log10(rho.T)
#     cs2 = ax2.contourf(theta, r, v_r.T, cmap=mycmap)
#     cs3 = ax3.contourf(theta, r, v_theta.T, cmap=mycmap)
#     # cs4 = ax4.contourf(theta, r, w[:,:,i], cmap=mycmap)
#     time = data['Time']
#     ax1.set_title('Density')
#     ax2.set_title('v_r')
#     ax3.set_title('v_theta')
#     ax4.set_title('Omega')
#     fig.suptitle('Time = '+str(time))
#     fig.colorbar(cs1,cax=cbar_ax1)
#     fig.colorbar(cs2,cax=cbar_ax2)
#     fig.colorbar(cs3,cax=cbar_ax3)
#     # fig.colorbar(cs4,cax=cbar_ax4)
#     print(i/len(files)*100)

# anim = animation.FuncAnimation(fig, animate, frames=len(files), interval=100, blit=True) #init_func=init
# anim.save('Animation_new.mp4', fps=1)

# # phi,v_phi = np.loadtxt('v_theta.txt',float,delimiter=',',unpack=True)
# # plt.scatter(phi,v_phi)
# # plt.show()