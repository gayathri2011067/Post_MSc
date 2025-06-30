import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import os, sys
import subprocess

# Paths to data and figures
data_path = "/mnt/0e8cdd72-0fe3-47e3-85d7-75b269878985/The other drive/Post_MSc/refined/data_files"
fig_path = "/mnt/0e8cdd72-0fe3-47e3-85d7-75b269878985/The other drive/Post_MSc/refined/figures"

trial_numbers = []
labels = []

# Read trial numbers and labels from command line arguments
t_num = 1
while True:
    try:
        trial_numbers.append(int(sys.argv[t_num]))
        labels.append(sys.argv[t_num + 1])
        t_num += 2
    except:
        break

# Create folder name for trials
folder_name = f'trials_({"_".join([str(i) for i in trial_numbers])})'

# Check if directory exists
if os.path.exists(f'{fig_path}/{folder_name}'):
    print(f"{folder_name} directory already exists")
    permission = input("Do you want to overwrite the existing directory? (y/n): ")
    if permission.lower() == 'n':
        print("Exiting the program")
        exit()
    elif permission.lower() == 'y':
        print("Overwriting the existing directory")

# Update figure path with folder name
fig_path = f'{fig_path}/{folder_name}'
os.makedirs(fig_path, exist_ok=True)

# Initialize lists to store data
z_lists = []
filename = 'z_values.txt'

# Load z values from the files
for trial_num in trial_numbers:
    file_path = f"{data_path}/trial_{trial_num}/{filename}"
    z = np.loadtxt(file_path)
    z_lists.append(np.copy(z))

# Initialize lists for storing magnetic field components and times
Brs = []
Bphis = []
B_strengths = []
B_average = []
times = []
filename1 = 'Br_final.txt'
filename2 = 'B_phi_final.txt'
filename3 = 'time.txt'
# Load magnetic field and time data
filename = 'pitch_angle.txt'

pitch_angles = []
plt.figure(figsize=(8, 6))  # Set figure size
colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']   # Darker color scheme for trials

# Loop through trials and plot pitch angles
for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    file_path = f"{data_path}/trial_{trial_num}/{filename}"  # File containing pitch angle data
    pitch_list = np.loadtxt(file_path)  # Load pitch angle data
    pitch_angles.append(np.copy(pitch_list))  # Save for later analysis
    
    # Plot pitch angle with consistent style
    plt.plot(z_lists[i], pitch_list[-1], linestyle='-', linewidth=2, color=colors[i % len(colors)])

plt.xlim(0, 0.8)
plt.ylim(0, -45)  # Set x-axis limits
#increase font size of scale markings
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r'z($\times$ 0.5 kpc)', fontsize=17)  # x-axis label
plt.ylabel(r'Pitch Angle (degrees)', fontsize=17)  # y-axis label
# Add grid for better visualization
plt.grid(True, alpha=0.3)

# Create legend with one line per model
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color=colors[i % len(colors)], linewidth=2, linestyle='-', label=rf'{labels[i]}')
    for i in range(len(trial_numbers))
]
plt.legend(handles=legend_elements, loc='upper left', fontsize=17)  # Add legend

# Save the plot
plt.savefig(f'{fig_path}/pitchang_vs_z_final.png', bbox_inches='tight')
plt.close()



pitch_angles_trimmed = []
plt.figure(figsize=(8, 6))  # Set figure size
colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']   # Darker color scheme for trials

# Loop through trials and plot pitch angles
for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    file_path = f"{data_path}/trial_{trial_num}/{filename}"  # File containing pitch angle data
    pitch_list = np.loadtxt(file_path)  # Load pitch angle data
    pitch_angles_trimmed.append(np.copy(pitch_list))  # Save for later analysis
    
    # Remove first three and last three grid points
    z_trimmed = z_lists[i][4:-4]
    pitch_trimmed = pitch_list[:,4:-4]
    
    # Plot trimmed pitch angle with consistent style
    plt.plot(z_trimmed, pitch_trimmed[-1], linestyle='-', linewidth=2, color=colors[i % len(colors)])

# plt.xlim(0, 1)
# plt.ylim(0, -45)  # Set x-axis limits
plt.xlabel(r'z($\times$ 0.5 kpc)', fontsize=17)  # x-axis label
plt.ylabel(r'Pitch Angle (degrees)', fontsize=17)  # y-axis label
# Add grid for better visualization
plt.grid(True, alpha=0.3)
plt.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
plt.axvline(x=0, color='grey', linestyle='--', alpha=0.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# Create legend with one line per model
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color=colors[i % len(colors)], linewidth=2, linestyle='-', label=rf'{labels[i]}')
    for i in range(len(trial_numbers))
]
plt.legend(handles=legend_elements, loc='upper left', fontsize=17)  # Add legend

# Save the plot
plt.savefig(f'{fig_path}/pitchang_vs_z_trimmed.png', bbox_inches='tight')
plt.close()



for i in range(len(trial_numbers)):
    
    trial_num = trial_numbers[i]

    file_path1 = f"{data_path}/trial_{trial_num}/{filename1}"
    file_path2 = f"{data_path}/trial_{trial_num}/{filename2}"
    file_path3 = f"{data_path}/trial_{trial_num}/{filename3}"
    

    Br_list = np.loadtxt(file_path1)
    Bphi_list = np.loadtxt(file_path2)
    time_list = np.loadtxt(file_path3)
   
    
    # Br_list=np.array(lines1,dtype=float)
    # Bphi_list=np.array(lines2,dtype=float)
    B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
    
    B_strengths.append(np.copy(B_strength))
    
    times.append(np.copy(time_list))
    
    Brs.append(np.copy(Br_list))
    Bphis.append(np.copy(Bphi_list))
   

plt.figure(figsize=(8, 6)) 

colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']   # Darker color scheme

for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]

    file_path1 = f"{data_path}/trial_{trial_num}/{filename1}"
    file_path2 = f"{data_path}/trial_{trial_num}/{filename2}"
    file_path3 = f"{data_path}/trial_{trial_num}/{filename3}"
    
    Br_list = np.loadtxt(file_path1)
    Bphi_list = np.loadtxt(file_path2)
    time_list = np.loadtxt(file_path3)
    
    B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
    
    B_strengths.append(np.copy(B_strength))
    times.append(np.copy(time_list))
    Brs.append(np.copy(Br_list))
    Bphis.append(np.copy(Bphi_list))
    
    # Plot Br and Bphi for each model using different colors
    plt.plot(z_lists[i], B_strength[-1], linestyle='-',linewidth=2, color=colors[i % len(colors)])
    # plt.plot(z_lists[i], Bphi_list[-1], linestyle='-',linewidth=2, color=colors[i % len(colors)])
    plt.xlim(-1, 1)
    # plt.ylim(1.70, 1.80)
    
    # plt.yscale('log')
    #increase font size of the labels

plt.xlabel(r'z($\times$ 0.5 kpc)',fontsize=17)
plt.ylabel(r'$B_{strength}/B_{eq}$',fontsize=17)
plt.grid(True, alpha=0.3)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
plt.axvline(x=0, color='grey', linestyle='--', alpha=0.5)
# plt.title(r"Saturated field components vs z")

# Create legend with one line per model
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color=colors[i % len(colors)], linewidth=2, linestyle='-', label=rf'{labels[i]}')
    for i in range(len(trial_numbers))
]
plt.legend(handles=legend_elements, loc='lower right', fontsize=17) 
plt.savefig(f'{fig_path}/Bst_vs_z_final.png')
plt.close()

#     plt.plot(z_lists[i], Br_list[-1], label=f'Br: {labels[i]}')
#     plt.plot(z_lists[i], Bphi_list[-1], label=f'Bphi: {labels[i]}')
#     plt.xlim(-1,1)
#     # plt.ylim(-9,9)
    
# plt.xlabel('z')
# plt.ylabel('Br, Bphi')
# plt.title(f"Br, Bphi vs z at t=({','.join([f'trial{trial_numbers[i]}:{times[i][-1]:.2e}' for j in range(len(trial_numbers))])})")

# plt.legend(loc='lower right', fontsize='small', framealpha=0.5)
# plt.savefig(f'{fig_path}/Br_Bphi_vs_z_final.png')
# plt.close()

# Plot Br and Bphi versus time at the midpoint of z
space_indices = []
plt.figure(figsize=(10, 5))  # Set figure size to rectangular
for i in range(len(trial_numbers)):
    space_idx = (Br_list.shape[1] - 1) // 2
    space_indices.append(space_idx)
    plt.plot(times[i], Brs[i][:, space_idx], label=f'Br: {labels[i]}')
    plt.plot(times[i], Bphis[i][:, space_idx], label=f'Bphi: {labels[i]}')

plt.xlabel(r'$t$')
plt.ylabel(r'$B_r, B_\phi$')
plt.title(r'$B_r, B_\phi \ \mathrm{vs} \ t \ \mathrm{at} \ z=(' +
          ', '.join([f'trial{trial_numbers[i]}:{int(z_lists[i][space_indices[i]])}' for i in range(len(trial_numbers))]) + r')$')
# plt.grid(True)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(f'{fig_path}/Br_Bphi_vs_time.png', bbox_inches='tight')
plt.close()

#'#1E91D6'

plt.figure(figsize=(8, 6))  # Set figure size to rectangular
B_strength_avg = np.mean(B_strengths, axis=2)
line_styles = ['-', '--', '-.', ':']  # Define different line styles
colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']   # Darker color scheme

for i in range(len(trial_numbers)):
    lab = labels[i]
    plt.plot(times[i], B_strength_avg[i], label=f'{lab}', linestyle='-', linewidth=2, color=colors[i % len(colors)])  # Increase line thickness

# plt.xlim(0, )
plt.yscale('log')
plt.ylim(10**-3, 2)
plt.xlabel('time(in Gyr)', fontsize=17)
plt.ylabel(r'${B}_{strength} / B_{eq}$', fontsize=17)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, alpha=0.3)
plt.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
plt.axvline(x=0, color='grey', linestyle='--', alpha=0.5)
# plt.title(f"Effect of different fluxes on magnetic field strength")
# plt.text(0.90, 1.05, r'$B_0=1\mu G$', transform=plt.gca().transAxes, verticalalignment='top')
plt.legend(loc='lower right', labels=[rf'{label}' for label in labels], fontsize=17)
plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()

plt.figure(figsize=(10, 5))  # Set figure size to rectangular
for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    
    plt.plot(times[i], B_strength_avg[i], label=f'{labels[i]}')
plt.xlim(0, 14)
plt.ylim(0, 0.6)
plt.xlabel(r'$t$')
plt.ylabel(r'$\mathrm{B}_{\mathrm{strength, avg}}$')
plt.title(r'$\mathrm{B}_{\mathrm{strength, avg}} \ \mathrm{vs} \ t$')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=17)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# Save or display the plot
plt.savefig(f'{fig_path}/B_strength_avg_vs_time.png', bbox_inches='tight')
# Alternatively, you can use plt.show() to display it in real-time


# Close the plot
plt.close()

# Plot alpha_m versus z for the last time step

# plt.figure(figsize=(8, 6))  # Set figure size to rectangular
# alpha_m_values = []
# filename = 'alpha_m.txt'
# colors = ['#3f8efc', '#940B92']  # Distinct colors for trials

# for i in range(len(trial_numbers)):
#     trial_num = trial_numbers[i]
#     file_path = f"{data_path}/trial_{trial_num}/{filename}"
#     alpha_m_list = np.loadtxt(file_path)
#     alpha_m_values.append(np.copy(alpha_m_list))

#     # Plot alpha_m with solid lines
#     plt.plot(z_lists[i], alpha_m_list[-1], linestyle='-', linewidth=2, color=colors[i])

# plt.axhline(y=0, color='grey', linestyle='--')
# plt.axvline(x=0, color='grey', linestyle='--')
# plt.xlim(-1,1)
# plt.ylim(-30, 30)
# plt.grid(True, alpha=0.3)
# plt.xlabel(r'z($\times$ 0.5 kpc)',fontsize=17)
# plt.ylabel(r'$\alpha_m$ (km s$^{-1}$)',fontsize=17)
# # plt.title(r'$\alpha_m \ \mathrm{vs} \ z$')

# # Create a simplified legend showing only color associations
# from matplotlib.lines import Line2D
# legend_elements = [
#     Line2D([0], [0], color=colors[i], linewidth=2, label=rf'{labels[i]}')
#     for i in range(len(trial_numbers))
# ]
# plt.legend(
#     handles=legend_elements,
#     loc='lower right',  # Place the legend in the lower right corner
#     fontsize=17,
#     framealpha=0.5  # Semi-transparent background for the legend
# )

# plt.savefig(f'{fig_path}/alpha_m_vs_z_final.png', bbox_inches='tight')
# plt.close()

















# plt.figure(figsize=(8,6))  # Set figure size to rectangular
# colors = ['#3f8efc', '#940B92']  # Colors for two trials

# for i in range(len(trial_numbers)):
#     trial_num = trial_numbers[i]

#     file_path1 = f"{data_path}/trial_{trial_num}/{filename1}"
#     file_path2 = f"{data_path}/trial_{trial_num}/{filename2}"
#     file_path3 = f"{data_path}/trial_{trial_num}/{filename3}"

#     Br_list = np.loadtxt(file_path1)
#     Bphi_list = np.loadtxt(file_path2)
#     time_list = np.loadtxt(file_path3)

#     B_strength = np.sqrt(Br_list**2 + Bphi_list**2)

#     B_strengths.append(np.copy(B_strength))
#     times.append(np.copy(time_list))
#     Brs.append(np.copy(Br_list))
#     Bphis.append(np.copy(Bphi_list))

#     # Plot Br as solid line
#     plt.plot(z_lists[i], Br_list[-1], linestyle='-', linewidth=2, color=colors[i], label=f'{labels[i]}')
#     # Plot Bphi as dashed line
#     plt.plot(z_lists[i], Bphi_list[-1], linestyle='--', linewidth=2, color=colors[i])

# plt.xlim(-1, 1)
# plt.xlabel(r'z($\times$ 0.5 kpc)',fontsize=17)
# plt.ylabel(r'$B_r, B_\phi(\mu$ G)',fontsize=17)
# # plt.title(r"Saturated field components vs z")

# # Create legend with one line per model
# from matplotlib.lines import Line2D
# legend_elements = [
#     Line2D([0], [0], color=colors[i % len(colors)], linewidth=2, linestyle='-', label=rf'{labels[i]}')
#     for i in range(len(trial_numbers))
# ]
# plt.legend(handles=legend_elements, loc='lower right', fontsize=17) 
# plt.savefig(f'{fig_path}/44Br_Bphi_vs_z_final.png')
# plt.close()







filename = 'pitch_angle.txt'

pitch_angles = []
plt.figure(figsize=(8, 6)) 
for i in range(len(trial_numbers)):
    
    trial_num = trial_numbers[i]

    file_path = f"{data_path}/trial_{trial_num}/{filename}"
    
    pitch_list = np.loadtxt(file_path)
    
    pitch_angles.append(np.copy(pitch_list))
    
    plt.plot(z_lists[i], pitch_list[-1], label=f'alpha_m: {labels[i]}')
    plt.xlim(-1,1)
    
plt.xlabel('z',fontsize=17)
plt.ylabel('alpha_m',fontsize=17)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title(f"alpha_m vs z at t=({','.join([f'trial{trial_numbers[i]}:{times[i][-1]:.2e}' for i in range(len(trial_numbers))])})")
plt.legend()
plt.savefig(f'{fig_path}/pitch_vs_z_final.png')
plt.close()








# Bstr_values = []
# filename = 'B_strength_avg.txt'

# for i in range(len(trial_numbers)):
#     trial_num = trial_numbers[i]
#     file_path = f"{data_path}/trial_{trial_num}/{filename}"
#     BSTR_list = np.loadtxt(file_path)
#     Bstr_values.append(np.copy(BSTR_list))

#     plt.plot(times[i], BSTR_list[i], label=f'B_strength: {labels[i]}')

# plt.xlabel(r'$t$')
# plt.ylabel(r'$B_avg$')
# plt.title(r'$B_avg \ \mathrm{vs} \ t$')
# # plt.grid(True)
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.savefig(f'{fig_path}/Bavg.png', bbox_inches='tight')
# plt.close()



fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
filename = 'alpha_m.txt'
alpha_m_values = []

# Flatten the 2D array of axes for easier iteration
axs = axs.flatten()

# for i in range(4):  # Assuming only 4 trials
#     trial_num = trial_numbers[i]
#     file_path = f"{data_path}/trial_{trial_num}/{filename}"
#     alpha_m_list = np.loadtxt(file_path)
#     alpha_m_values.append(np.copy(alpha_m_list))
    
#     # Plot on the corresponding subplot
#     axs[i].plot(z_lists[i], alpha_m_list[-1], label=f'alpha_m: {labels[i]}')
#     axs[i].set_xlim(-1,1)
#     axs[i].set_ylim(-10, 10)
#     axs[i].axhline(y=0, color='grey', linestyle='--')
#     axs[i].axvline(x=0, color='grey', linestyle='--')
#     axs[i].set_xlabel(r'$z$')
#     axs[i].set_ylabel(r'$\alpha_m$')
#     axs[i].set_title(f'Trial {trial_num}')
#     axs[i].legend(loc='upper right')
    
# # Adjust the layout so subplots do not overlap
# plt.tight_layout()

# # Save the entire figure as a single image
# plt.savefig(f'{fig_path}/alpha_m_vs_z_subplots.png', bbox_inches='tight')
# plt.close()


# First set of subplots: B_strength vs t at specific z for each trial
fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
axs = axs.flatten()

# for i in range(4):  # Assuming only 4 trials
#     trial_num = trial_numbers[i]
    
#     # Plot B_strength vs time at the specific z location for each trial
#     axs[i].plot(times[i], B_strengths[i][:, space_indices[i]], label=f'{labels[i]}')
#     # axs[i].set_xlim(0, 10)
#     axs[i].set_yscale('log')
#     axs[i].set_xlabel(r'$t$')
#     axs[i].set_ylabel(r'$\mathrm{B}_{\mathrm{strength}}$')
#     axs[i].set_title(f'Trial {trial_num} at z={int(z_lists[i][space_indices[i]])}')
#     axs[i].legend(loc='upper right')

# plt.tight_layout()
# plt.savefig(f'{fig_path}/B_strength_vs_time_subplots.png', bbox_inches='tight')
# plt.close()

# # Second set of subplots: B_strength_avg vs t for each trial
# fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
# axs = axs.flatten()

# B_strength_avg = np.mean(B_strengths, axis=2)

# for i in range(4):  # Assuming only 4 trials
#     trial_num = trial_numbers[i]
    
#     # Plot B_strength_avg vs time for each trial
#     axs[i].plot(times[i], B_strength_avg[i], label=f'{labels[i]}', color='red')
#     axs[i].set_xlim(0, 7)
#     # axs[i].set_ylim(0.01,1)
#     axs[i].set_yscale('log')
#     axs[i].set_xlabel(r'$t \ (\mathrm{computational \ units})$')
#     axs[i].set_ylabel(r'$\log \ \langle B \rangle / B_0$')
#     axs[i].set_title(f'Trial {trial_num}')
#     axs[i].legend(loc='lower right')
    
#     # Add ticklines on all 4 sides of the grid
#     axs[i].tick_params(axis='both', which='both', direction='in', top=True, right=True)

# plt.tight_layout()
# plt.savefig(f'{fig_path}/B_strength_avg_vs_time_subplots.png', bbox_inches='tight')
# plt.close()


# fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
# alpha_m_filename = 'alpha_m.txt'
# turb_vel_filename = 'turb_vel.txt'
# alpha_m_values = []
# turb_vel_values = []
# n1 = 5000
# total_t = 10
# t_values = 2.5  #NOTE: Change this value to plot at different time steps
# t_index = int(t_values * n1 / total_t)
# # Flatten the 2D array of axes for easier iteration
# axs = axs.flatten()

# for i in range(4):  # Assuming only 4 trials
#     trial_num = trial_numbers[i]
    
#     # Load alpha_m data
#     alpha_m_file_path = f"{data_path}/trial_{trial_num}/{alpha_m_filename}"
#     alpha_m_list = np.loadtxt(alpha_m_file_path)
#     alpha_m_values.append(np.copy(alpha_m_list))
    
#     # Load turbulent velocity data
#     turb_vel_file_path = f"{data_path}/trial_{trial_num}/{turb_vel_filename}"
#     with open(turb_vel_file_path, 'r') as f:
#         lines = f.readlines()
#     turb_vel = np.array(lines, dtype=float)
#     turb_vel_values.append(turb_vel)

#     # Plot alpha_m on the corresponding subplot
#     axs[i].plot(z_lists[i], alpha_m_list[t_index], label=f'alpha_m: {labels[i]}')
#     axs[i].axhline(y=0, color='grey', linestyle='--')
#     axs[i].axvline(x=0, color='grey', linestyle='--')

#     # Plot turbulent velocity on the same subplot
#     axs[i].plot(z_lists[i], turb_vel, linestyle='--', color='black', alpha=0.8, label='RMS turb velocity')

#     # Set y-limits based on the subplot index
#     if i == 0:
#         axs[i].set_ylim(-100, 100)
#         axs[i].set_xlim(-1, 1)
#     else:
#         axs[i].set_ylim(-40, 40)
#         axs[i].set_xlim(-1, 1)

#     # Add labels, title, and legend
#     axs[i].set_xlabel(r'$z$')
#     axs[i].set_ylabel(r'$\alpha_m$ / RMS turb velocity')
#     axs[i].set_title(f'Trial {trial_num}')
#     axs[i].legend(loc='upper left')

# # Adjust the layout so subplots do not overlap
# plt.tight_layout()

# # Save the entire figure as a single image
# plt.savefig(f'{fig_path}/alpha_m_and_turb_vel_vs_z_subplots.png', bbox_inches='tight')
# plt.close()


filename = 'alpha_m.txt'

alpha_m_values = []
for i in range(len(trial_numbers)):
    
    trial_num = trial_numbers[i]

    file_path = f"{data_path}/trial_{trial_num}/{filename}"
    
    alpha_m_list = np.loadtxt(file_path)
    
    alpha_m_values.append(np.copy(alpha_m_list))
    
    plt.plot(z_lists[i], alpha_m_list[-1], label=f'{labels[i]}')
    plt.xlim(-1,1)
    
plt.xlabel('z',fontsize=17)
plt.ylabel('alpha_m',fontsize=17)
# plt.title(f"alpha_m vs z at t=({','.join([f'trial{trial_numbers[i]}:{times[i][-1]:.2e}' for i in range(len(trial_numbers))])})")
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(f'{fig_path}/alpha_m_vs_z_final.png')
plt.close()

plt.figure(figsize=(8, 6)) 

colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']  # Darker color scheme

for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]

    file_path1 = f"{data_path}/trial_{trial_num}/{filename1}"
    file_path2 = f"{data_path}/trial_{trial_num}/{filename2}"
    file_path3 = f"{data_path}/trial_{trial_num}/{filename3}"
    
    Br_list = np.loadtxt(file_path1)
    Bphi_list = np.loadtxt(file_path2)
    time_list = np.loadtxt(file_path3)
    
    B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
    
    B_strengths.append(np.copy(B_strength))
    times.append(np.copy(time_list))
    Brs.append(np.copy(Br_list))
    Bphis.append(np.copy(Bphi_list))
    
    # Plot Br and Bphi for each model using different colors
    plt.plot(z_lists[i], Br_list[-1], linestyle='-',linewidth=2, color=colors[i % len(colors)])
    plt.plot(z_lists[i], Bphi_list[-1], linestyle='--',linewidth=2, color=colors[i % len(colors)])
    plt.xlim(-1, 1)
    # plt.yscale('log')
plt.xlabel(r'z($\times$ 0.5 kpc)',fontsize=17)
plt.grid(True, alpha=0.3)
plt.ylabel(r'$B_r, B_\phi(\mu$ G)',fontsize=17)
plt.axhline(y=0, color='grey', linestyle='-', alpha=0.3)
plt.xticks(fontsize=14)
plt.yticks(fontsize=12)
# plt.title(r"Saturated field components vs z")

# Create legend with one line per model
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color=colors[i % len(colors)], linewidth=2, linestyle='-', label=rf'{labels[i]}')
    for i in range(len(trial_numbers))
]
# plt.hlines(0, -1, 1, colors='grey', linestyles='--')
plt.legend(handles=legend_elements, loc='lower right', fontsize=17) 
plt.savefig(f'{fig_path}/Br_Bphi_vs_z_final.png')
plt.close()


import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

fig, axs = plt.subplots(2, 2, figsize=(8, 6), sharex=True, sharey=True)  # 2x2 grid of subplots
colors =['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518'] # Darker, distinct colors

# Loop through trials and plot each in its subplot
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

fig, axs = plt.subplots(2, 2, figsize=(8,6), sharex=True)  # 2x2 grid of subplots, shared x-axis

colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']   # Darker, distinct colors

# Loop through trials and plot each in its subplot
# for i in range(len(trial_numbers)):
#     trial_num = trial_numbers[i]
#     row, col = divmod(i, 2)  # Determine subplot position (row, col)

#     file_path = f"{data_path}/trial_{trial_num}/{filename}"  # File containing alpha_m data
#     alpha_m_list = np.loadtxt(file_path)  # Load alpha_m data
#     alpha_m_values.append(np.copy(alpha_m_list))  # Save for later analysis

#     # Plot alpha_m in the appropriate subplot
#     axs[row, col].plot(
#         z_lists[i], alpha_m_list[-1],
#         linewidth=2, color=colors[i % len(colors)]
#     )
#     axs[row, col].grid(True, alpha=0.3)  # Subtle gridlines
#     axs[row, col].set_title(f'{labels[i]}', fontsize=17)  # Label each subplot with its model
#     axs[row, col].set_xlim(-1, 1)  # Set common x-axis limits
#     axs[row, col].set_ylim(-1, 1)  # Set common y-axis limits

# # Set common labels and title
# fig.text(0.5, 0.02, r'z ($\times$ 0.5 kpc)', ha='center', fontsize=16)  # Common x-axis label
# fig.text(0.02, 0.5, r'$\alpha_m$ (km s$^{-1}$)', va='center', rotation='vertical', fontsize=16)  # Common y-axis label
# # fig.suptitle(r"Profile of $\alpha_m$ vs z at Saturation", fontsize=16)  # Common title

# # Adjust spacing between subplots
# plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])

# # Save the figure
# plt.savefig(f'{fig_path}/alpha_m_vs_z_subplots.png', dpi=300)
# plt.close()

plt.figure(figsize=(8, 6))  # Single plot with defined size

colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']  # Darker, distinct colors

# Loop through trials and plot each on the same plot
for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    file_path = f"{data_path}/trial_{trial_num}/{filename}"  # File containing alpha_m data
    alpha_m_list = np.loadtxt(file_path)  # Load alpha_m data
    alpha_m_values.append(np.copy(alpha_m_list))  # Save for later analysis

    # Plot alpha_m on the same figure
    plt.plot(
        z_lists[i], alpha_m_list[-1],
        label=f'{labels[i]}',  # Label for each trial
        linewidth=2, color=colors[i % len(colors)]
    )

# Add grid, labels, and legend
plt.grid(True, alpha=0.3)  # Subtle gridlines
plt.xlim(-1, 1)  # Set common x-axis limits
plt.ylim(-20, 20)  # Set common y-axis limits
plt.xlabel(r'z ($\times$ 0.5 kpc)', fontsize=17)
plt.ylabel(r'$\alpha_m$ (km s$^{-1}$)', fontsize=17)
# plt.title(r"Profile of $\alpha_m$ vs z at Saturation", fontsize=16)
plt.legend(loc='upper right', fontsize=10)

# Save the figure
plt.tight_layout()
plt.savefig(f'{fig_path}/alpha_m_wrt_trial.png', dpi=300)
plt.close()



#ABCD '#f39237', '#9c3848', '#f39237','#1E91D6'
#CDEF '#f39237','#1E91D6', '#940B92', '#f39237'


if len(trial_numbers) == 2:
    plt.figure(figsize=(8, 6))  # Single plot with defined size
    colors = ['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']   # Distinct colors for two trials

    for i in range(len(trial_numbers)):
        trial_num = trial_numbers[i]
        file_path = f"{data_path}/trial_{trial_num}/{filename}"  # File containing alpha_m data
        alpha_m_list = np.loadtxt(file_path)  # Load alpha_m data
        alpha_m_values.append(np.copy(alpha_m_list))  # Save for later analysis

        # Plot alpha_m on the same figure
        plt.plot(
            z_lists[i], alpha_m_list[-1],
            label=f'Model: {labels[i]}',  # Label for each trial
            linewidth=2, color=colors[i % len(colors)]
        )

    # Add grid, labels, and legend
    plt.grid(True, alpha=0.3)  # Subtle gridlines
    plt.xlim(-1, 1)  # Set common x-axis limits
    plt.ylim(-20, 20)  # Set common y-axis limits
    plt.xlabel(r'z ($\times$ 0.5 kpc)', fontsize=17)
    plt.ylabel(r'$\alpha_m$ (km s$^{-1}$)', fontsize=17)
    plt.title(r"Profile of $\alpha_m$ vs z at Saturation", fontsize=16)
    plt.legend(loc='upper right', fontsize=10)

    # Save the figure with "2" prefix
    plt.tight_layout()
    plt.savefig(f'{fig_path}/2_alpha_m_wrt_trial.png', dpi=300)
    plt.close()
    
    
    
if len(trial_numbers) == 2:
    plt.figure(figsize=(8, 6))  # Single plot for 2 trials
    colors =['#007FFF','#E5385E','#17B169','#FEBE10','#162227','#FF7518']  # Distinct colors for two trials

    for i in range(len(trial_numbers)):
        trial_num = trial_numbers[i]
        file_path1 = f"{data_path}/trial_{trial_num}/{filename1}"
        file_path2 = f"{data_path}/trial_{trial_num}/{filename2}"

        Br_list = np.loadtxt(file_path1)
        Bphi_list = np.loadtxt(file_path2)

        # Plot Br and Bphi
        plt.plot(z_lists[i], Br_list[-1], linestyle='--', linewidth=2, color=colors[i], label=f'Br: {labels[i]}')
        plt.plot(z_lists[i], Bphi_list[-1], linestyle='-', linewidth=2, color=colors[i], label=f'Bphi: {labels[i]}')

    plt.xlim(-1, 1)
    plt.xlabel(r'z ($\times$ 0.5 kpc)', fontsize=17)
    plt.ylabel(r'$B_r, B_\phi (\mu$ G)', fontsize=17)
    plt.grid(True, alpha=0.3)
    plt.axhline(y=0, color='grey', linestyle='--', alpha=0.5)
    plt.axvline(x=0, color='grey', linestyle='--', alpha=0.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=12)
    plt.legend(loc='lower right', fontsize=10)
    plt.savefig(f'{fig_path}/2_Br_Bphi_vs_z_final.png')
    plt.close()



print(np.shape(B_strengths))
print(np.shape(times))
print(np.shape(Brs))
print(np.shape(alpha_m_list))










alpha_m_values = []  # Store all alpha_m data

plt.figure(figsize=(8, 6))  # Create figure

for i in range(len(trial_numbers)):
    trial_num = trial_numbers[i]
    file_path = f"{data_path}/trial_{trial_num}/{filename}"  # File containing alpha_m data
    alpha_m_list = np.loadtxt(file_path)  # Load alpha_m data
    alpha_m_values.append(np.copy(alpha_m_list))  # Save for later analysis

    # Extract time series for alpha_m at index 51
    alpha_m_at_51 = alpha_m_list[:, 79]

    # Plot alpha_m vs. time
    plt.plot(
        times[i], alpha_m_at_51,  
        label=f'{labels[i]}',  # Label for each trial
        linewidth=2, color=colors[i % len(colors)]
    )

# Add grid, labels, and legend
plt.grid(True, alpha=0.3)  # Subtle gridlines
plt.xlabel(r'time(in Gyr)', fontsize=17)
plt.ylabel(r'$\alpha_m$ (km s$^{-1}$) at z $\approx$ 0.5 h', fontsize=17)
# plt.title(r"Evolution of $\alpha_m$ at Index 51", fontsize=16)
plt.legend(loc='lower right', fontsize=14)

# # Save the figure
# plt.tight_layout()
plt.savefig(f'{fig_path}/alpha_m_vs_time_at_51.png', dpi=300)
plt.close()
















