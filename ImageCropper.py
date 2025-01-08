from PIL import Image

def crop_image(image_path, left, top, right, bottom, output_path):
    """
    Crop the image based on the specified coordinates and save as JPEG.

    Args:
        image_path (str): Path to the image file.
        left (int): Left coordinate of the crop area.
        top (int): Top coordinate of the crop area.
        right (int): Right coordinate of the crop area.
        bottom (int): Bottom coordinate of the crop area.
        output_path (str): Path to save the cropped image.

    Returns:
        None
    """
    # Open the image
    img = Image.open(image_path)
    
    # Crop the image
    cropped_img = img.crop((left, top, right, bottom))
    
    # Save the cropped image as JPEG
    cropped_img.save(output_path, "JPEG")

# Example usage:
image_path = "Before.jpg"
left = 250  # left coordinate of the crop area
top = 75   # top coordinate of the crop area
right = 807 # right coordinate of the crop area
bottom = 893 # bottom coordinate of the crop area
output_path = "BeforeN.jpg"  # Path to save the cropped image

'''
crop_image(image_path, left, top, right, bottom, output_path)
crop_image("New.jpg", left, top, right, bottom, "NewN.jpg")

crop_image(image_path, 1095, top, 1652, bottom, "BeforeDB.jpg")
crop_image("New.jpg", 1095, top, 1652, bottom, "NewDB.jpg")
'''

import matplotlib.pyplot as plt
from PIL import Image

def view_image(image_path):
    # Open the image using PIL
    img = Image.open(image_path)
    
    plt.figure(facecolor='black')
    # Display the image using matplotlib
    plt.imshow(img)
    plt.title('Cropped Image')
    plt.axis('off')  # Turn off axis numbers and ticks
    plt.show()


from PIL import Image

def crop_and_save(input_image_path, left, top, right, bottom, output_image_path, quality=91 ):
    """
    Crop a region from the input image based on the specified coordinates
    and save the cropped region as a JPEG image with the given quality.
    
    Parameters:
        input_image_path (str): Path to the input image file.
        output_image_path (str): Path to save the output cropped image.
        left (int): Left coordinate of the cropping region.
        top (int): Top coordinate of the cropping region.
        right (int): Right coordinate of the cropping region.
        bottom (int): Bottom coordinate of the cropping region.
        quality (int): JPEG image quality (0-100). Default is 95.
    """
    try:
        # Open the input image
        input_image = Image.open(input_image_path)
        
        # Crop the specified region
        cropped_image = input_image.crop((left, top, right, bottom))
        
        # Save the cropped image as JPEG with specified quality
        cropped_image.save(output_image_path, "JPEG", quality=quality)
        
        print(f"Image cropped and saved as {output_image_path} with quality={quality}")
        
    except Exception as e:
        print(f"Error: {e}")


# Call the function to crop and save the image
crop_and_save("Before.jpg", 330, 97, 1065, 1188, "BeforeN.jpg")

crop_and_save("New.jpg",  330, 97, 1065, 1188, "NewN.jpg")

crop_and_save("Before.jpg", 1458, 97, 2190, 1188, "BeforeDB.jpg")
crop_and_save("New.jpg", 1458, 97, 2190, 1188, "NewDB.jpg")

# Example usage to view the saved cropped images
view_image("BeforeDB.jpg")  # View the first cropped image
'''
view_image("NewN.jpg")     # View the second cropped image
view_image("BeforeDB.jpg")  # View the first cropped image with different coordinates
view_image("NewDB.jpg")     # View the second cropped image with different coordinates
'''
