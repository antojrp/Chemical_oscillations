# Import necessary modules
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import io

# Parameters
file_out = "chemical_reactions"  # Output file name (without extension)
interval = 200  # Time between frames in milliseconds
save_to_file = True  # False: show the animation on screen, True: save to a file
dpi = 150  # Output video quality (dots per inch)

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
frames_data_u = load_data("Chemical_oscillations_u.txt")

# Load data from the "v" file
frames_data_v = load_data("Chemical_oscillations_v.txt")

# Create matrices for RGB color coding
nframes = min(len(frames_data_u), len(frames_data_v))

# Initialize the frames_data list
frames_data = []

# Create an array of zeros with the specified shape
size = len(frames_data_u[0])
zeros_array = np.zeros((size,size))

# Iterate over each frame
for j_frame in range(nframes):
    # Stack the three arrays along the last axis
    disordered_array = np.stack((frames_data_v[j_frame], zeros_array, frames_data_u[j_frame]), axis=-1)
    # Append the stacked array to the frames_data list
    frames_data.append(disordered_array)
    
maximo = np.max(frames_data)
minimo = np.min(frames_data)

# Normalize the image data to the range [0, 1]
frames_data_normalized = [(frame_data - minimo) / (maximo - minimo) for frame_data in frames_data]

# Create the figure and axis objects
fig, ax = plt.subplots()
ax.axis("off")  # Do not show the axes
            
# Display the first frame
im = ax.imshow(frames_data_normalized[0])

# Function to update the animation
def update(j_frame, frames_data_normalized, im):
    im.set_data(frames_data_normalized[j_frame])
    return im,

# Create the animation
animation = FuncAnimation(fig, update, fargs=(frames_data_normalized, im), frames=nframes, blit=True, interval=interval)

# Save or show the animation
if save_to_file:
    animation.save("{}.mp4".format(file_out), dpi=dpi)
else:
    plt.show()

