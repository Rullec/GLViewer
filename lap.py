import cv2
import numpy as np
from matplotlib import pyplot as plt
import time
from copy import deepcopy


def add_continous_noise(png):
    st = time.time()
    max_png = np.max(png)
    threshold = int(max_png + np.random.rand() * 8)

    all_contour = cv2.Canny(png, 10, 20)
    png[png == 0] = threshold

    all_contour_dilation = cv2.dilate(all_contour, np.ones([4, 4]))

    blurred_img = cv2.GaussianBlur(png, (7, 7), 2)
    mask = all_contour_dilation != 0 
    mask[np.random.random(blurred_img.shape) < 0.4] = 0
    png[mask] = blurred_img[mask]
    png[png == threshold] = 0
    ed = time.time()
    print(f"add_continous_noise cost {ed - st} s")
    return png


if __name__ == "__main__":
    img_path = "assets/cufanggezi0_from_cons.png"
    st = time.time()

    png = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
    plt.subplot(1, 3, 1)
    plt.imshow(png)
    png = add_continous_noise(png)
    plt.subplot(1, 3, 2)
    plt.imshow(png)
    plt.subplot(1, 3, 3)
    plt.imshow(cv2.imread("assets/cufanggezi_nofix/0.png", cv2.IMREAD_GRAYSCALE))

    plt.show()
    exit()
    max_png = np.max(png)
    threshold = int(max_png + np.random.rand() * 20)
    print(f"max png {max_png} thre {threshold}")

    png[png == 0] = threshold
    all_contour = cv2.Canny(png, 10, 20)
    print(f"canny cost {time.time() - st} s")
    # outer_contour = cv2.Canny(png, 100, 200)
    # inner_bd = cv2.bitwise_and(cv2.bitwise_not(outer_contour), all_contour)
    # plt.imshow( cv2.dilate(png, 10, 10))
    # plt.show()

    all_contour_dilation = cv2.dilate(all_contour, np.ones([4, 4]))

    blurred_img = cv2.GaussianBlur(png, (7, 7), 2)
    print(f"dilation and gaussian blur cost {time.time() - st} s")
    new_png = deepcopy(png)
    print(f"")
    new_png[all_contour_dilation != 0] = blurred_img[all_contour_dilation != 0]
    plt.subplot(2, 2, 1)
    png[png == threshold] = 0
    ed5 = time.time()
    print(f"total cost {ed5 - st} s")
    plt.imshow(png)
    plt.title("raw img")
    plt.subplot(2, 2, 2)
    plt.imshow(all_contour_dilation)
    plt.title("all contour")
    plt.subplot(2, 2, 3)
    new_png[new_png == threshold] = 0
    plt.imshow(new_png)
    plt.title("blurred img")

    plt.subplot(2, 2, 4)
    diff = new_png - png
    diff[diff != 0] = 1
    plt.imshow(diff)
    plt.title("diff")
    plt.show()