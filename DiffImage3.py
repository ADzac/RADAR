import numpy as np
from PIL import Image
import matplotlib.pyplot as plt



file = "Before.jpg"
fileN = "BeforeN.jpg"
fileDB = "BeforeDB.jpg"
file_2 = "New.jpg"
file_2N = "NewN.jpg"
file_2DB = "NewDB.jpg"

def find_image_differences(image1_path, image2_path):
    # Open images
    img1 = Image.open(image1_path)
    img2 = Image.open(image2_path)

    # Convert images to RGB mode if not already
    img1 = img1.convert("RGB")
    img2 = img2.convert("RGB")

    # Check if images have the same dimensions
    if img1.size != img2.size:
        raise ValueError("Images must have the same dimensions.")

    # Convert images to numpy arrays to compare
    np_img1 = np.array(img1)
    np_img2 = np.array(img2)

    # Find differences
    differences = np.where(np_img1 != np_img2)

    return differences

def plot_differences(image1_path, image2_path):
    differences = find_image_differences(image1_path, image2_path)
    plt.scatter(differences[1], differences[0], color='r', s=1)
    plt.ylim((0 ,820))
    plt.gca().invert_yaxis()
    plt.show()

def calculate_difference_percentage(image1_path, image2_path):
    differences = find_image_differences(image1_path, image2_path)
    total_pixels = np.prod(np.array(Image.open(image1_path)).size)
    num_differences = len(differences[0])
    return (num_differences / total_pixels) * 100

# Example usage:
image1 = "Realiste.jpg"
image2 = "Realistic.jpg"

#plot_differences(image1, image2)
#print("Difference Percentage for a figure:", calculate_difference_percentage(image1, image2))
'''
plot_differences("RealisteN.jpg", "RealisticN.jpg")
print("Difference Percentage in Linear:", calculate_difference_percentage("RealisteN.jpg", "RealisticN.jpg"))

plot_differences("RealisteDB.jpg", "RealisticDB.jpg")
print("Difference Percentage in dB:", calculate_difference_percentage("RealisteDB.jpg", "RealisticDB.jpg"))
'''

plot_differences(fileN, file_2N)
print("Difference Percentage in Linear:", calculate_difference_percentage(fileN, file_2N))

plot_differences(fileDB,file_2DB)
print("Difference Percentage in dB:", calculate_difference_percentage(fileDB, file_2DB))