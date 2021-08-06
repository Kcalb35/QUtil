import matplotlib.pyplot as plt
import numpy as np

prefix = "build/test/"

if __name__ == '__main__':
    names = ["SAC", "DAC", "ECR", "DBG", "DAG", "DRN"]
    scale = [50, 12, 1, 1, 1, 50]
    for i in range(len(names)):
        plt.close()
        data = np.transpose(np.loadtxt(prefix + names[i] + ".txt"))
        plt.plot(data[0], data[1], label="E1")
        plt.plot(data[0], data[2], label="E2")
        plt.plot(data[0], data[3] / scale[i], label="nac")
        plt.title(names[i])
        plt.show()
