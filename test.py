import sim_kinect 
from matplotlib import pyplot as plt
import time
if __name__ == "__main__":
    kinect = sim_kinect.sim_kinect()
    kinect.Init()
    # print(f"[debug] build kinect noise succ")
    # print(f"[debug] init kinect noise succ")
    st = time.time()
    res = kinect.LoadAndCalculate("./assets/cufanggezi0_from_cons.png")
    ed = time.time()


    print(res.shape, f"cost {(ed - st) * 1e3} ms")
    plt.imshow(res)
    plt.show()