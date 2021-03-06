import sim_kinect
import cv2
from matplotlib import pyplot as plt
import time
import numpy as np
if __name__ == "__main__":
    # data = np.random.normal(0, 1.0 / 6, size = 1000)
    # plt.hist(data)
    # plt.show()
    kinect = sim_kinect.sim_kinect()
    kinect.Init()
    st = time.time()
    data_lst = []
    kinect.SetBaseline(0.03)
    for i in range(10):
        # -----show baseline difference-----
        # cur_baseline = 0.01 * i
        # kinect.SetBaseline(cur_baseline)
        # cur_focal = kinect.GetFocalLength()
        # res = kinect.LoadAndCalculate("./assets/cufanggezi0_from_cons.png")
        # plt.subplot(3, 4, i + 1)
        # plt.title(f"baseline {cur_baseline} cur_focal {cur_focal}")
        # plt.imshow(res)
        # data_lst.append(res)
        cur_focal = 100 * i
        kinect.SetFocalLength(cur_focal)
        cur_baseline = kinect.GetBaseline()
        img = cv2.imread("./assets/cufanggezi0_from_cons.png", cv2.IMREAD_GRAYSCALE)
        res = kinect.ApplyKinectHoleQuantizationNoise(img)
        plt.subplot(3, 4, i + 1)
        plt.title(f"cur_focal {kinect.GetFocalLength()} baseline {kinect.GetBaseline():.2f} ")
        plt.imshow(res)
        data_lst.append(res)
        


    ed = time.time()
    cost_ms = (ed - st) * 1e3
    print(res.shape, f"total cost {cost_ms} ms, avg {cost_ms / len(data_lst)} ms")

    # plt.imshow(res)
    plt.show()
    plt.savefig("output.png")