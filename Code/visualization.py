# Import necessary modules
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import io

# Parameters
file_out = "chemical_reactions"  # Output file name (without extension)
interval = 50  # Time between frames in milliseconds
save_to_file = True  # False: show the animation on screen, True: save to a file
dpi = 300  # Output video quality (dots per inch)

# Function to load data from a file
def load_data(file_name):
    with open(file_name, "r") as f:
        data_str = f.read()

    frames_data = []
    for frame_data_str in data_str.split("\n\n"):
        if frame_data_str.strip():  # Check if the frame data is not empty
            frame_data = np.loadtxt(io.StringIO(frame_data_str), delimiter=",")
            frames_data.append(frame_data)
    return frames_data

# Load data from the "u" file
f_u = load_data("Chemical_oscillations_u.txt")

# Load data from the "v" file
f_v = load_data("Chemical_oscillations_v.txt")

# Create matrices for RGB color coding
nframes = min(len(f_u), len(f_v))

# Initialize the frames_data list
frames_data = []

# Create an array of zeros with the specified shape
size = len(f_u[0])
zeros_array = 1 - np.zeros((size,size))

# Iterate over each frame
for j in range(nframes):
    # Stack the three arrays along the last axis
    
    max_f_u = np.max(f_u[j])        
    max_f_v = np.max(f_v[j])     
    min_f_u = np.min(f_u[j])        
    min_f_v = np.min(f_v[j])   
    max_f = np.max([max_f_u, max_f_v])
    
    
    if ((max_f_v - min_f_v) <= 0.000001) & ((max_f_v) > 0):
        norm_v = f_v[j] / max_f_v
    elif (max_f_v) > 0:
        norm_v = (f_v[j] - min_f_v) / (max_f_v - min_f_v)
    else:
        norm_v = f_v[j]
        
        
    if ((max_f_u - min_f_u) <= 0.000001) & ((max_f_u) > 0):
        norm_u = f_u[j] / max_f_u
    elif (max_f_u) > 0:
        norm_u = (f_u[j] - min_f_u) / (max_f_u - min_f_u)
    else:
        norm_u = f_u[j]
    
    
    disordered_array = np.stack( ( norm_v, zeros_array, norm_u), axis=-1)
    # Append the stacked array to the frames_data list
    frames_data.append(disordered_array)
    
#alpha = 1
    
#maximo = alpha * np.max(frames_data)

# Normalize the image data to the range [0, 1]
# frames_data_normalized =[1 - (frame_data) / (maximo) for frame_data in frames_data]

# Create the figure and axis objects
fig, ax = plt.subplots()
ax.axis("off")  # Do not show the axes

# IMPORTANTE: AÑADIR _normalizes a frames_data
            
# Display the first frame
im = ax.imshow(frames_data[0])

# Function to update the animation
def update(j_frame, frames_data_normalized, im):
    im.set_data(frames_data[j_frame])
    return im,

# Create the animation
animation = FuncAnimation(fig, update, fargs=(frames_data, im), frames=nframes, blit=True, interval=interval)

# Save or show the animation
if save_to_file:
    animation.save("{}.mp4".format(file_out), dpi=dpi)
'''
else:
    plt.show()
'''
