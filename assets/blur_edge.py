import cv2
import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":
    origin_name = "./cufanggezi0_from_cons.png"
    output_name = "./cufanggezi0_from_cons.blurred.png"
    image = cv2.imread(origin_name, cv2.IMREAD_GRAYSCALE)
    blur_thickness = 5
    gaus_kernel = 9
    blurred_img = cv2.GaussianBlur(image, (gaus_kernel, gaus_kernel), 0)
    mask = np.zeros(image.shape, np.uint8)

    # gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    from copy import deepcopy
    gray = deepcopy(image)
    thresh = cv2.threshold(gray, 60, 255, cv2.THRESH_BINARY)[1]
    # print(f"thresh {thresh}")
    # thresh = thresh[2]
    # exit()
    contours, hierarchy = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL,
                                           cv2.CHAIN_APPROX_SIMPLE)

    cv2.drawContours(mask, contours, -1, (255, ), blur_thickness)
    output = np.where(mask == np.array([255]), blurred_img, image)

    plt.subplot(1, 2, 1)
    plt.imshow(image)
    plt.subplot(1, 2, 2)
    plt.imshow(output)
    plt.show()
    cv2.imwrite(output_name, output)
    print(f"output {output_name}")