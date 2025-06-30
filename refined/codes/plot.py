import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import os
import subprocess
from scipy.optimize import curve_fit

data_path = "/mnt/0e8cdd72-0fe3-47e3-85d7-75b269878985/The other drive/Post_MSc/refined/run_files"
fig_path = "/mnt/0e8cdd72-0fe3-47e3-85d7-75b269878985/The other drive/Post_MSc/refined/figures"
data_save_path = "/mnt/0e8cdd72-0fe3-47e3-85d7-75b269878985/The other drive/Post_MSc/refined/data_files"

passed_args = argparse.ArgumentParser()

passed_args.add_argument(
    '-t',
    '--trial',
    type=int,
    default=None,
    help='trial number',
)

args = passed_args.parse_args()

if args.trial is None or type(args.trial) != int:
    save_path = Path(data_save_path)
    exisiting_trials = [d for d in save_path.iterdir() if d.is_dir()]
    trial_numbers = []
    for d in exisiting_trials:
        if d.name.startswith('trial_'):
            trial_numbers.append(int(d.name.split('trial_')[-1]))
    # print(f"Existing trials: {trial_numbers}")♥
    last_trial = max(trial_numbers)
    trial_num = last_trial + 1
else:
    trial_num = args.trial
    
# check if the trial directory exists
if os.path.exists(f'{data_save_path}/trial_{trial_num}'):
    print(f"trial_{trial_num} directory already exists")
    permission = input("Do you want to overwrite the existing directory? (y/n): ")
    if permission == 'n' or permission == 'N':
        print("Exiting the program")
        exit()
    elif permission == 'y' or permission == 'Y':
        print("Overwriting the existing directory")
else:
    print(f"trial_{trial_num} directory does not exist. Creating new directory")

fig_path = f'{fig_path}/trial_{trial_num}'
data_save_path = f'{data_save_path}/trial_{trial_num}'

os.makedirs(data_save_path, exist_ok=True)
os.makedirs(fig_path, exist_ok=True)

#import txt file 
# filename = 'eta_fz_values.txt'
# file_path = f"{data_path}/{filename}"

# with open(file_path,'r') as f:
#     lines=f.readlines()
# # print(lines)

# eta=np.array(lines,dtype=float)
# #import txt file

filename = 'z_values.txt'
file_path = f"{data_path}/{filename}"

with open(file_path,'r') as f:
    lines=f.readlines()
# print(lines)

z=np.array(lines,dtype=float)
np.savetxt(f'{data_save_path}/z_values.txt', z)

filename = 'turb_vel.txt'
file_path = f"{data_path}/{filename}"

with open(file_path,'r') as f:
    lines=f.readlines()
# print(lines)

turb_vel=np.array(lines,dtype=float)
np.savetxt(f'{data_save_path}/turb_vel.txt', turb_vel)
# plt.plot(z,eta)
# plt.xlabel('z')
# plt.ylabel('eta')
# plt.title('eta vs z')
# plt.xlim(-15,15)
# plt.savefig(f'{fig_path}/eta_vs_z.png')
# plt.close()
#import txt file
filename = 'alpha_values.txt'

file_path = f"{data_path}/{filename}"

try:
    with open(file_path,'r') as f:
        lines=f.readlines()
    # print(lines)

    alpha=np.array(lines,dtype=float)
    np.savetxt(f'{data_save_path}/alpha_values.txt', alpha)

    plt.plot(z,alpha)
    plt.xlabel('z')
    plt.ylabel('alpha')
    plt.title('alpha vs z')
    # plt.xlim(-15,15)
    plt.savefig(f'{fig_path}/alpha_vs_z.png')
    plt.close()
except:
    print('alpha_values.txt file not found')

#import txt file
filename1 = 'Br_ini.txt'
filename2 = 'B_phi_ini.txt'
file_path1 = f"{data_path}/{filename1}"
file_path2 = f"{data_path}/{filename2}"

try:
    with open(file_path1,'r') as f:
        lines1=f.readlines()
    # print(lines1)
    with open(file_path2,'r') as f:
        lines2=f.readlines()
    # print(lines2)
    Br=np.array(lines1,dtype=float)
    Bphi=np.array(lines2,dtype=float)

    plt.plot(z,Br)
    plt.plot(z,Bphi)
    plt.xlabel('z')
    plt.ylabel('Br, Bphi')
    plt.title('Br, Bphi vs z - initial')
    plt.savefig(f'{fig_path}/Br_Bphi_vs_z_initial.png')
    plt.close()
except:
    print('Br_ini.txt or B_phi_ini.txt file not found')

#import txt file
filename1 = 'Br_final.txt'
filename2 = 'B_phi_final.txt'
filename3 = 'time.txt'
file_path1 = f"{data_path}/{filename1}"
file_path2 = f"{data_path}/{filename2}"
file_path3 = f"{data_path}/{filename3}"


try:
    with open(file_path1,'r') as f:
        lines1=f.readlines()
    # print(lines1)
    with open(file_path2,'r') as f:
        lines2=f.readlines()
    # print(lines2)
    with open(file_path3,'r') as f:
        lines3=f.readlines()
    # print(lines3)

    Br_list = []
    for line in lines1:
        line = line.strip()
        line = line.split()
        curr = np.array(line, dtype=float)
        Br_list.append(curr)
    Br_list = np.array(Br_list)

    Bphi_list = []
    for line in lines2:
        line = line.strip()
        line = line.split()
        curr = np.array(line, dtype=float)
        Bphi_list.append(curr)
    Bphi_list = np.array(Bphi_list)


    # Br_list=np.array(lines1,dtype=float)
    # Bphi_list=np.array(lines2,dtype=float)
    time_list=np.array(lines3,dtype=float)
    # print(time_list)
    np.savetxt(f'{data_save_path}/time.txt', time_list)
    np.savetxt(f'{data_save_path}/Br_final.txt', Br_list)
    np.savetxt(f'{data_save_path}/B_phi_final.txt', Bphi_list)
    # print(Br_list.shape)
    plt.plot(z, Br_list[-1], label='Br')
    plt.plot(z, Bphi_list[-1], label='Bphi')

    # plt.plot(z, Br_list[-4], label='Br')
    # plt.plot(z, Bphi_list[-4], label='Bphi')
    print(Br_list.shape)
    # print(Bphi_list[1])
    # plt.plot(z, Br_list[1], label='Br')
    # plt.plot(z, Bphi_list[1], label='Bphi')
    # plt.plot(z, Br_list[2], label='Br')
    # plt.plot(z, Bphi_list[2], label='Bphi')
    plt.xlabel('z')
    plt.ylabel('Br, Bphi')
    plt.title(f'Br, Bphi vs z at t={time_list[-1]}')
    # plt.xlim(-0.25,0.25)
    plt.axvline(x=1, color='r', linestyle='--')
    plt.axvline(x=-1, color='r', linestyle='--')
    plt.axhline(y=0, color='g', linestyle='--')
    plt.legend()
    plt.savefig(f'{fig_path}/Br_Bphi_vs_z_final.png')
    plt.close()
except:
    print('Br_final.txt or B_phi_final.txt or time.txt file not found')


#plot of Br and Bphi at against time
plt.plot(time_list, Br_list[:,-1], label='Br')
# print(time_list)
plt.plot(time_list, Bphi_list[:,-1], label='Bphi')
plt.xlabel('time')
plt.ylabel('Br, Bphi')
plt.title('Br, Bphi vs time')
plt.legend()
plt.savefig(f'{fig_path}/Br_Bphi_vs_time.png')
plt.close()




# plt.plot(time_list,alpha)
# plt.xlabel('time')
# plt.ylabel('alpha')
# plt.title('alpha vs time')
# plt.xlim(-15,15)
# plt.savefig(f'{fig_path}/alpha_vs_time.png')
# plt.close()
B_strength = np.sqrt(Br_list**2 + Bphi_list**2)
# 320-640, 100-200
# selected_time_list = time_list[320:640]
# selected_B_strength = B_strength[320:640, 51]

np.savetxt(f'{data_save_path}/B_strength.txt', B_strength)

# Plot B_strength at a specific z index





# Select the time indexes between 320 and 640 and their corresponding B_strength


# Take the gradient of log(B_strength) with respect to the selected time range and save it
# grad_B_strength = np.gradient(np.log(selected_B_strength), selected_time_list)
# np.savetxt(f'{data_save_path}/grad_B_strength.txt', grad_B_strength)

# # Plot the growth rates against the selected time range
# plt.plot(selected_time_list, grad_B_strength)
# plt.xlabel('time')
# plt.ylabel('Growth Rate')
# plt.title('Growth Rate vs time')
# plt.savefig(f'{fig_path}/Growth_Rate_vs_time.png')
# plt.close()

# # Calculate the derivative of B_strength with respect to the selected time range
# dB_strength_dt = np.gradient(selected_B_strength, selected_time_list)

# # Calculate the mean of the derivative
# mean_dB_strength_dt = np.mean(dB_strength_dt)

# Save the derivative and mean to files
# np.savetxt(f'{data_save_path}/dB_strength_dt.txt', dB_strength_dt)
# with open(f'{data_save_path}/mean_dB_strength_dt.txt', 'w') as f:
#     f.write(str(mean_dB_strength_dt))

# Plot B_strength vs time with limits
plt.plot(time_list, B_strength[:,51])
plt.xlabel('time')
# plt.xlim(0,10)
# #mark the region of selected time range and the corresponding B_strength by shading
# plt.axvspan(selected_time_list[0], selected_time_list[-1], color='grey', alpha=0.5)
# plt.axhspan(selected_B_strength[0], selected_B_strength[-1], color='grey', alpha=0.5)
#write the value of mean_dB_strength_dt on the plot
# plt.text(0.5, 0.5, f'mean_dB_strength_dt = {mean_dB_strength_dt:.4f}', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes)


# plt.yscale('log')
plt.ylabel('B_strength/B0')
plt.title('B_strength vs time')
plt.yscale('log')


# plt.axhline(y=0.001, color='g', linestyle='--')
# plt.axvline(x=2.5, color='r', linestyle='--')
plt.savefig(f'{fig_path}/B_strength_vs_time.png')
plt.close()





#plot alpha_m vs time


# plt.plot(time_list,alpha)
# plt.xlabel('z')
# plt.ylabel('alpha')
# plt.title('alpha vs z')
# # plt.xlim(-15,15)
# plt.savefig(f'{fig_path}/alpha_vs_time.png')
# plt.close()



# Find growth rate by linear fitting log(B_strength) vs time at mid z 

# mid_z_index = 51 # Select the mid z index
# tolerable_error = 10
# patiance_cells = 50
# max_cells = 500

# log_B_mid = np.log10(B_strength[:, mid_z_index])

# # # find the contant slope region
# first_derivative = np.gradient(log_B_mid, time_list)
# print(first_derivative)
# second_derivative = np.gradient(first_derivative, time_list)
# constant_slope_region = np.argwhere(np.abs(second_derivative) < tolerable_error).flatten()


# # put nan values in the regions of second derivative > 1
# second_derivative[np.abs(second_derivative) > 40] = np.nan
# # plot the derivatives and the regions of constant slope
# plt.plot(time_list, first_derivative, label='First Derivative')
# plt.plot(time_list, second_derivative, label='Second Derivative')
# plt.axhline(y=0, color='grey', linestyle='--', linewidth=1, alpha=0.5)
# plt.xlabel('time')
# plt.ylabel('Second Derivative')
# plt.title('Second Derivative of log(B_strength) vs time at mid z')
# plt.legend()
# plt.savefig(f'{fig_path}/second_derivative_log_B_strength_vs_time_at_mid.png')
# # plt.show()
# plt.close()

# # find continuous regions
# continuous_regions = []
# start_idx = constant_slope_region[0]
# for i in range(1, len(constant_slope_region)):
#     if constant_slope_region[i] - constant_slope_region[i-1] != 1:
#         end_idx = constant_slope_region[i-1] + 1 # to be used as end index in numpy slicing
#         if end_idx - start_idx > patiance_cells and end_idx - start_idx < max_cells:
#             continuous_regions.append((start_idx, end_idx))
#         start_idx = constant_slope_region[i]
#     else:
#         if i == len(constant_slope_region) - 1:
#             end_idx = constant_slope_region[i] + 1
#             if end_idx - start_idx > patiance_cells and end_idx - start_idx < max_cells:
#                 continuous_regions.append((start_idx, end_idx))
            
# slopes = []
# for region in continuous_regions:
#     slope, _ = np.polyfit(time_list[region[0]:region[1]], log_B_mid[region[0]:region[1]], 1)
#     slopes.append(slope)

# # plot log(B_strength) vs time at mid z
# plt.plot(time_list, np.log10(B_strength[:, mid_z_index]))
# # plot the second derivative 

# # highlight the regions of constant slope
# for region in continuous_regions:
#     plt.axvspan(time_list[region[0]], time_list[region[1]-1], color='grey', alpha=0.5)
#     plt.text(time_list[region[0]], 0.5, f'{slopes[continuous_regions.index(region)]:.4f}', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transData)
# plt.xlabel('time')
# plt.ylabel('log(B_strength)')
# plt.title('log(B_strength) vs time at mid z')
# plt.savefig(f'{fig_path}/log_B_strength_vs_time_at_mid.png')
# # plt.show()
# plt.close()
 




# import numpy as np
# import matplotlib.pyplot as plt

# mid_z_index = 51  # Select the mid z index
# tolerance_slope = 0.05  # Slope variation tolerance
# min_consecutive_points = 10  # Minimum number of consecutive points for constant slope
# min_slope_threshold = 0.1  # Minimum absolute slope to consider as exponential
# max_time_cutoff = 4  # Ignore points beyond this time for growth detection

# # Compute log of B_strength at the midplane
# log_B_mid = np.log10(B_strength[:, mid_z_index])
# #if last 10 values of bstrength is same, then it is saturating
# if np.all(log_B_mid[-10:] == log_B_mid[-1]):
#     print("The field is saturating")
# else:
#     print("The field is not saturating")
# # Compute first derivative at every 5th time step
# time_sampled = time_list[::5]  
# first_derivative = np.gradient(log_B_mid[::5], time_sampled)

# # Identify constant slope regions separately for growth and decay
# growth_indices = []
# decay_indices = []
# current_growth = []
# current_decay = []

# for i in range(1, len(first_derivative)):
#     slope = first_derivative[i]
#     prev_slope = first_derivative[i - 1]

#     if np.abs(slope - prev_slope) < tolerance_slope:
#         if slope > min_slope_threshold and time_sampled[i] < max_time_cutoff:
#             current_growth.append(i)
#         elif slope < -min_slope_threshold:
#             current_decay.append(i)
#     else:
#         if len(current_growth) >= min_consecutive_points:
#             growth_indices.extend(current_growth)
#         if len(current_decay) >= min_consecutive_points:
#             decay_indices.extend(current_decay)

#         current_growth = []
#         current_decay = []

# # Ensure the last detected region is included
# if len(current_growth) >= min_consecutive_points:
#     growth_indices.extend(current_growth)
# if len(current_decay) >= min_consecutive_points:
#     decay_indices.extend(current_decay)

# growth_indices = np.array(growth_indices)
# decay_indices = np.array(decay_indices)

# # Compute average slopes
# avg_growth_slope = np.mean(first_derivative[growth_indices]) if len(growth_indices) > 0 else None
# avg_decay_slope = np.mean(first_derivative[decay_indices]) if len(decay_indices) > 0 else None

# print("Average growth slope:", avg_growth_slope)
# print("Average decay slope:", avg_decay_slope)

# # ==== Plot the results ====
# plt.plot(time_list, B_strength[:, mid_z_index], label="B Strength at z={}".format(mid_z_index))

# # Highlight constant slope growth regions in red
# if len(growth_indices) > 0:
#     plt.scatter(time_sampled[growth_indices], B_strength[growth_indices * 5, mid_z_index], 
#                 color='red', label="Exponential Growth", zorder=3)

# # Highlight constant slope decay regions in blue
# if len(decay_indices) > 0:
#     plt.scatter(time_sampled[decay_indices], B_strength[decay_indices * 5, mid_z_index], 
#                 color='blue', label="Exponential Decay", zorder=3)

# plt.xlabel("Time")
# plt.ylabel("Magnetic Field Strength")
# plt.yscale('log')
# plt.legend()
# plt.savefig(f'{fig_path}/B_strength_vs_time_at_mid_z.png')
# plt.show()






# import numpy as np
# import matplotlib.pyplot as plt

# mid_z_index = 51  # Select the mid z index
# tolerance_slope = 0.05  # Slope variation tolerance
# min_consecutive_points = 10  # Minimum number of consecutive points for constant slope
# min_slope_threshold = 0.1  # Minimum absolute slope to consider as exponential

# # Compute log of B_strength at the midplane
# log_B_mid = np.log10(B_strength[:, mid_z_index])

# # Compute first derivative at every 5th time step
# time_sampled = time_list[::5]  
# first_derivative = np.gradient(log_B_mid[::5], time_sampled)

# # Identify the peak (saturation point)
# peak_index = np.argmax(B_strength[:, mid_z_index])

# # Determine whether the field is growing, saturating, or decaying
# is_saturating = (peak_index > 0) and (peak_index < len(B_strength[:, mid_z_index]) - 1)

# growth_indices = []
# decay_indices = []
# current_growth = []
# current_decay = []

# for i in range(1, len(first_derivative)):
#     slope = first_derivative[i]
#     prev_slope = first_derivative[i - 1]
    
#     if np.abs(slope - prev_slope) < tolerance_slope:
#         if slope > min_slope_threshold and not is_saturating:
#             current_growth.append(i)
#         elif slope < -min_slope_threshold and is_saturating:
#             current_decay.append(i)
#     else:
#         if len(current_growth) >= min_consecutive_points:
#             growth_indices.extend(current_growth)
#         if len(current_decay) >= min_consecutive_points:
#             decay_indices.extend(current_decay)
        
#         current_growth = []
#         current_decay = []

# # Ensure last detected regions are included
# if len(current_growth) >= min_consecutive_points and not is_saturating:
#     growth_indices.extend(current_growth)
# if len(current_decay) >= min_consecutive_points and is_saturating:
#     decay_indices.extend(current_decay)

# growth_indices = np.array(growth_indices)
# decay_indices = np.array(decay_indices)

# # Compute average slopes
# avg_growth_slope = np.mean(first_derivative[growth_indices]) if len(growth_indices) > 0 else None
# avg_decay_slope = np.mean(first_derivative[decay_indices]) if len(decay_indices) > 0 else None

# print("Average growth slope:", avg_growth_slope)
# print("Average decay slope:", avg_decay_slope)

# # ==== Plot the results ====
# plt.plot(time_list, B_strength[:, mid_z_index], label="B Strength at z={}".format(mid_z_index))

# # Highlight constant slope growth regions in red (only if no saturation detected)
# if len(growth_indices) > 0:
#     plt.scatter(time_sampled[growth_indices], B_strength[growth_indices * 5, mid_z_index], 
#                 color='red', label="Exponential Growth", zorder=3)

# # Highlight constant slope decay regions in blue (only if no saturation detected)
# if len(decay_indices) > 0 and not is_saturating:
#     plt.scatter(time_sampled[decay_indices], B_strength[decay_indices * 5, mid_z_index], 
#                 color='blue', label="Exponential Decay", zorder=3)

# plt.xlabel("Time")
# plt.ylabel("Magnetic Field Strength")
# plt.yscale('log')
# plt.legend()
# plt.savefig(f'{fig_path}/B_strength_vs_time_at_mid_z.png')
# plt.show()




































# import numpy as np
# import matplotlib.pyplot as plt

# mid_z_index = 51  # Select the mid z index
# tolerance_slope = 0.05  # Slope variation tolerance
# min_consecutive_points = 10  # Minimum number of consecutive points for constant slope
# min_slope_threshold = 0.1  # Minimum absolute slope to consider as exponential

# # Compute log of B_strength at the midplane
# log_B_mid = np.log10(B_strength[:, mid_z_index])

# # Check for saturation
# if np.all(log_B_mid[-10:] == log_B_mid[-1]):
#     print("The field is saturating")
#     max_time_cutoff = 4  # Default cutoff
#     check_growth = True
# else:
#     print("The field is not saturating")
#     max_time_cutoff = time_list[-1]  # Extend cutoff to last time value
#     check_growth = False

# # Compute first derivative at every 5th time step
# time_sampled = time_list[::5]  
# first_derivative = np.gradient(log_B_mid[::5], time_sampled)

# # Identify constant slope regions separately for growth and decay
# growth_indices = []
# decay_indices = []
# current_growth = []
# current_decay = []

# for i in range(1, len(first_derivative)):
#     slope = first_derivative[i]
#     prev_slope = first_derivative[i - 1]

#     if np.abs(slope - prev_slope) < tolerance_slope:
#         if check_growth and slope > min_slope_threshold and time_sampled[i] < max_time_cutoff:
#             current_growth.append(i)
#         elif slope < -min_slope_threshold:
#             current_decay.append(i)
#     else:
#         if check_growth and len(current_growth) >= min_consecutive_points:
#             growth_indices.extend(current_growth)
#         if len(current_decay) >= min_consecutive_points:
#             decay_indices.extend(current_decay)

#         current_growth = []
#         current_decay = []

# # Ensure the last detected region is included
# if check_growth and len(current_growth) >= min_consecutive_points:
#     growth_indices.extend(current_growth)
# if len(current_decay) >= min_consecutive_points:
#     decay_indices.extend(current_decay)

# growth_indices = np.array(growth_indices)
# decay_indices = np.array(decay_indices)

# # Compute average slopes
# avg_growth_slope = np.mean(first_derivative[growth_indices]) if check_growth and len(growth_indices) > 0 else None
# avg_decay_slope = np.mean(first_derivative[decay_indices]) if len(decay_indices) > 0 else None

# print("Average growth slope:", avg_growth_slope)
# print("Average decay slope:", avg_decay_slope)

# # ==== Plot the results ====
# plt.plot(time_list, B_strength[:, mid_z_index], label="B Strength at z={}".format(mid_z_index))

# # Highlight constant slope growth regions in red
# if check_growth and len(growth_indices) > 0:
#     plt.scatter(time_sampled[growth_indices], B_strength[growth_indices * 5, mid_z_index], 
#                 color='red', label="Exponential Growth", zorder=3)

# # Highlight constant slope decay regions in blue
# if len(decay_indices) > 0:
#     plt.scatter(time_sampled[decay_indices], B_strength[decay_indices * 5, mid_z_index], 
#                 color='blue', label="Exponential Decay", zorder=3)

# plt.xlabel("Time")
# plt.ylabel("Magnetic Field Strength")
# plt.yscale('log')
# plt.legend()
# plt.savefig(f'{fig_path}/B_strength_vs_time_at_mid_z.png')
# plt.show()






















import matplotlib.cm as cm


# plot of alpha_m final vs z
filename = 'alpham_final.txt'
file_path = f"{data_path}/{filename}"

try:
    with open(file_path,'r') as f:
        lines=f.readlines()

    alpha_list = []
    for line in lines:
        line = line.strip()
        line = line.split()  
        curr = np.array(line, dtype=float)
        alpha_list.append(curr)
    alpha_m=np.array(alpha_list,dtype=float)
    np.savetxt(f'{data_save_path}/alpha_m.txt', alpha_m)
    n1 = 5000
    total_t = 20
    t_values = np.linspace(0, 20, 20)  # More time points for smoother evolution
    t_indices = [int(t * n1 / total_t) - 1 for t in t_values]

    colors = cm.viridis(np.linspace(0, 1, len(t_values)))  # Dark, high-contrast colormap

    plt.figure(figsize=(8, 6))

    # Plot each time step
    for i, (t_index, t_val) in enumerate(zip(t_indices, t_values)):
        linestyle = '-' if i < len(t_values) // 2 else '--'  # Solid for early times, dashed for later times
        plt.plot(z, alpha_m[t_index],color=colors[i], linestyle=linestyle,linewidth=2, label=f't={t_val:.2f}')

    # Add horizontal and vertical reference lines
    plt.axhline(y=0, color='grey', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(x=0, color='grey', linestyle='--', linewidth=1, alpha=0.5)

    # Labels and formatting
    plt.xlim(-1, 1)
    plt.xlabel(r'z($\times$ 0.5 kpc)', fontsize=17)
    plt.ylabel(r'$\alpha_m$ (km s$^{-1}$)', fontsize=17)
    plt.title(r'$\alpha_m$ evolution', fontsize=17)
    plt.legend(title='Time(Gyr)',fontsize=12)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    # Save figure
    plt.savefig(f'{fig_path}/alpha_m_vs_z_evolution.png')
    plt.close()
except:
    print('alpha_m.txt file not found')





#plot bstrength squared on y and -alpha m on x
# print(np.shape(B_strength))

plt.plot( -alpha_m[-1],B_strength[-1]**2)
plt.ylabel('B_strength_avg')
plt.xlabel('-alpha_m')
plt.title('B_strength_avg vs -alpha_m')
plt.savefig(f'{fig_path}/B_strength_avg_vs_alpha_m.png')
plt.close()






filename = 'alpham_final.txt'
file_path = f"{data_path}/{filename}"
try:
    with open(file_path,'r') as f:
        lines=f.readlines()

    alpha_mm_list = []
    for line in lines:
        line = line.strip()
        line = line.split()  
        curr = np.array(line, dtype=float)
        alpha_mm_list.append(curr)
    alpha_mm=np.array(alpha_mm_list,dtype=float)
    np.savetxt(f'{data_save_path}/aallppmm.txt', alpha_mm)
    # print(alpha_tot)
    plt.plot(z, alpha_mm[-1],label='alpha_m',color='red',linewidth=2)
    plt.xlabel(r'z($\times$ 0.5 kpc)',fontsize=17)
    plt.ylabel(r'$\alpha_m$ (km s$^{-1}$)',fontsize=17)
    plt.vlines([0, 1, -1], -1, 1, color='black', linestyle='--', linewidth=2, alpha=1)
    # plt.title('aallppmm vs z')
    plt.savefig(f'{fig_path}/aallppmm_vs_z.png')
    plt.close()
except:
    print('alpham_final.txt file not found')

#read and save alpha_tot similarly
filename = 'alp_tot.txt'
file_path = f"{data_path}/{filename}"

try:
    with open(file_path,'r') as f:
        lines=f.readlines()

    alpha_tot_list = []
    for line in lines:
        line = line.strip()
        line = line.split()  
        curr = np.array(line, dtype=float)
        alpha_tot_list.append(curr)
    alpha_tot=np.array(alpha_tot_list,dtype=float)
    np.savetxt(f'{data_save_path}/alpha_tot.txt', alpha_tot)
    # print(alpha_tot)
    plt.plot(z, alpha_tot[-1])
    plt.xlabel('z')
    plt.ylabel('alpha_tot')
    plt.title('alpha_tot vs z')
    plt.savefig(f'{fig_path}/alpha_tot_vs_z.png')
    plt.close()
except:
    print('alpha_tot.txt file not found')
    

  
# #  plot alpha and alpha total in same plot against time
plt.plot(time_list,alpha_mm[:,51],label='alpha_m')
plt.plot(time_list,alpha_tot[:,51],label='alpha_total')
plt.xlabel('time')
plt.ylabel('alpha, alpha_tot')
plt.title('alpha, alpha_tot vs time')
plt.legend()
plt.savefig(f'{fig_path}/alpha_alpha_tot_vs_time.png')
plt.close()
   
    

#write a code to find and plot pitch angles |p| = arctan |Br/Bphi|

# Calculate the pitch angle and limit the result to the range [-pi/2, pi/2]
p = np.degrees(np.arctan((Br_list / Bphi_list)))  # Convert the angle to degrees
np.savetxt(f'{data_save_path}/pitch_angle.txt', p)

# Plot the pitch angle vs z
plt.plot(z, p[-1])
plt.xlabel('z')
plt.ylabel('Pitch Angle (degrees)')
plt.title('Pitch Angle vs z')
plt.savefig(f'{fig_path}/pitch_angle_vs_z.png')
plt.close()





















# import numpy as np
# import matplotlib.pyplot as plt

# mid_z_index = 51  # Select the mid z index
# tolerance_slope = 0.05  # Slope variation tolerance
# min_consecutive_points = 10  # Minimum number of consecutive points for constant slope
# min_slope_threshold = 0.1  # Minimum absolute slope to consider as exponential
# max_time_cutoff = 8  # Ignore points beyond this time for growth detection

# # Compute log of B_strength at the midplane
# log_B_mid = np.log10(B_strength[:, mid_z_index])

# # Check for saturation
# is_saturating = np.all(log_B_mid[-10:] == log_B_mid[-1])

# # Compute first derivative at every 5th time step
# time_sampled = time_list[::5]  
# first_derivative = np.gradient(log_B_mid[::5], time_sampled)

# # Identify constant slope regions based on saturation status
# indices_to_detect = []
# current_segment = []

# target_slope_sign = 1 if is_saturating else -1  # Growth (1) if saturating, Decay (-1) if not

# for i in range(1, len(first_derivative)):
#     slope = first_derivative[i]
#     prev_slope = first_derivative[i - 1]

#     if np.abs(slope - prev_slope) < tolerance_slope and slope * target_slope_sign > min_slope_threshold:
#         current_segment.append(i)
#     else:
#         if len(current_segment) >= min_consecutive_points:
#             indices_to_detect.extend(current_segment)
#         current_segment = []

# # Ensure the last detected region is included
# if len(current_segment) >= min_consecutive_points:
#     indices_to_detect.extend(current_segment)

# indices_to_detect = np.array(indices_to_detect)

# # Compute average slope based on the detected region
# avg_detected_slope = np.mean(first_derivative[indices_to_detect]) if len(indices_to_detect) > 0 else None

# # Print results
# if is_saturating:
#     print("The field is saturating.")
#     print("Average growth slope:", avg_detected_slope)
# else:
#     print("The field is not saturating.")
#     print("Average decay slope:", avg_detected_slope)

# # ==== Plot the results ====
# plt.plot(time_list, B_strength[:, mid_z_index], label="B Strength at z={}".format(mid_z_index))

# # Highlight detected regions
# if len(indices_to_detect) > 0:
#     plt.scatter(time_sampled[indices_to_detect], B_strength[indices_to_detect * 5, mid_z_index], 
#                 color='red' if is_saturating else 'blue',
#                 label="Exponential Growth" if is_saturating else "Exponential Decay", 
#                 zorder=3)

# plt.xlabel("Time")
# plt.ylabel("Magnetic Field Strength")
# plt.yscale('log')
# plt.legend()
# plt.savefig(f'{fig_path}/B_strength_vs_time_at_mid_z.png')
# plt.show()










import numpy as np
import matplotlib.pyplot as plt

mid_z_index = 51  # Select the mid z index
tolerance_slope = 0.04  # Slope variation tolerance
min_consecutive_points = 10  # Minimum number of consecutive points for constant slope
min_slope_threshold = 0.004#0.001  # Minimum absolute slope to consider as exponential
tolerance_slope = 0.5  # Slope variation tolerance#for really steep cases

# Compute log of B_strength at the midplane
log_B_mid = np.log10(B_strength[:, mid_z_index])

# Define a tolerance for saturation detection
saturation_tolerance = 1e-5

# Check for saturation: If the last few values are nearly constant
if np.all(np.abs(log_B_mid[-3:] - log_B_mid[-1]) < saturation_tolerance):
    print("The field is saturating")
    max_time_cutoff = 4  # Default cutoff
    detect_growth = True  # Look for both growth and decay
else:
    print("The field is not saturating")
    max_time_cutoff = time_list[-1]  # Extend cutoff to last time value
    detect_growth = False  # Only look for decay

# Compute first derivative at every 5th time step
time_sampled = time_list[::5]  
first_derivative = np.gradient(log_B_mid[::5], time_sampled)

# Identify constant slope regions
growth_indices = []
decay_indices = []
current_decay = []
current_growth = []

for i in range(1, len(first_derivative)):
    slope = first_derivative[i]
    prev_slope = first_derivative[i - 1]

    if np.abs(slope - prev_slope) < tolerance_slope:
        if detect_growth and slope > min_slope_threshold and time_sampled[i] < max_time_cutoff:
            current_growth.append(i)
        elif slope < -min_slope_threshold:
            current_decay.append(i)
    else:
        if len(current_growth) >= min_consecutive_points:
            growth_indices.extend(current_growth)
        if len(current_decay) >= min_consecutive_points:
            decay_indices.extend(current_decay)

        current_growth = []
        current_decay = []

# Ensure the last detected region is included
if len(current_growth) >= min_consecutive_points and detect_growth:
    growth_indices.extend(current_growth)
if len(current_decay) >= min_consecutive_points:
    decay_indices.extend(current_decay)

growth_indices = np.array(growth_indices)
decay_indices = np.array(decay_indices)

# Compute average slopes
if detect_growth:
    avg_growth_slope = np.mean(first_derivative[growth_indices]) if len(growth_indices) > 0 else None
else:
    avg_growth_slope = None

avg_decay_slope = np.mean(first_derivative[decay_indices]) if len(decay_indices) > 0 else None

print("Average growth slope:", avg_growth_slope)
print("Average decay slope:", avg_decay_slope)

# ==== Plot the results ====
plt.plot(time_list, B_strength[:, mid_z_index], label="B Strength at z={}".format(mid_z_index))

# Highlight constant slope decay regions in blue
if len(decay_indices) > 0:
    plt.scatter(time_sampled[decay_indices], B_strength[decay_indices * 5, mid_z_index], 
                color='blue', label="Exponential Decay", zorder=3)

# Highlight constant slope growth regions in red (only if it was analyzed)
if detect_growth and len(growth_indices) > 0:
    plt.scatter(time_sampled[growth_indices], B_strength[growth_indices * 5, mid_z_index], 
                color='red', label="Exponential Growth", zorder=3)

plt.xlabel("Time")
plt.ylabel("Magnetic Field Strength")
plt.yscale('log')
plt.legend()
plt.savefig(f'{fig_path}/B_strength_vs_time_at_mid_z.png')
# plt.show()


#GROWTH OR DECAY, WHATEVER AVERAGE SLOPE, WRITE IT INTO A FILE
with open(f'{data_save_path}/growth_decay_slope.txt', 'w') as f:
    if detect_growth:
        f.write(f"Average growth slope: {avg_growth_slope}\n")
    f.write(f"Average decay slope: {avg_decay_slope}\n")
    
    
import numpy as np
import matplotlib.pyplot as plt

# Load data
filename = 'alpham_final.txt'
file_path = f"{data_path}/{filename}"

try:
    with open(file_path, 'r') as f:
        lines = f.readlines()

    alpha_list = []
    for line in lines:
        line = line.strip().split()  
        curr = np.array(line, dtype=float)
        alpha_list.append(curr)

    alpha_m = np.array(alpha_list, dtype=float)
    np.savetxt(f'{data_save_path}/alpha_m.txt', alpha_m)

    # Define time values
    n1 = 5000  # Number of time steps
    total_t = 20  # Total time in Gyr
    time_values = np.linspace(0, 20, 1)  # Extracted time points

    # Compute mean or representative alpha_m value at each time step
    alpha_mean = np.mean(alpha_m, axis=1)  # Average across z for each time step

    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(time_values, alpha_mean[:len(time_values)], marker='o', linestyle='-', linewidth=2, markersize=6, color='b', label=r'$\langle \alpha_m \rangle$')

    # Labels and formatting
    plt.xlabel('Time (Gyr)', fontsize=17)
    plt.ylabel(r'$\alpha_m$ (km s$^{-1}$)', fontsize=17)
    plt.title(r'$\alpha_m$ Evolution Over Time', fontsize=17)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=14)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    
    # Save figure
    plt.savefig(f'{fig_path}/alpha_m_vs_time.png')
    # plt.show()
#plot the bstrength at z=0 and alpha_m at z=0 against time
    # plt.plot(time_list, B_strength[:,51])



except FileNotFoundError:
    print('alpham_final.txt file not found')

z

plt.plot(time_list, alpha_m[:,51])
plt.xlabel('time')
plt.ylabel('alpha_m')
plt.title('Time evolution of alpha_m at z=0')

plt.savefig(f'{fig_path}/B_strength_alpha_m_vs_time.png')
plt.close()



filename = 'analytic.txt'
file_path = f"{data_path}/{filename}"

with open(file_path,'r') as f:
    lines=f.readlines()
# print(lines)

analyticb=np.array(lines,dtype=float)
np.savetxt(f'{data_save_path}/analytic_values.txt', analyticb)
value=np.mean(analyticb)
# print(value)
plt.plot(z,B_strength[51])
#plot a line with value
plt.axhline(y=value, color='r', linestyle='--')
plt.xlabel('z')
plt.ylabel('B_strength')
plt.title('Analytic vs B_strength')
plt.savefig(f'{fig_path}/analytic_vs_B_strength.png')
plt.close()









save_path = f"{data_save_path}/good_data.txt"
np.savetxt(save_path, B_strength[:, 51], header="B_strength[:, 51]")

# Append growth/decay slopes to the file
with open(save_path, "a") as f:
    f.write("\n--- Summary Statistics ---\n")
    f.write(f"Average growth slope: {avg_growth_slope}\n")
    f.write(f"Average decay slope: {avg_decay_slope}\n")
