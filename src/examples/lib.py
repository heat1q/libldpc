import ctypes
import threading
import time
import numpy as np

snrs = [0.200, 0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000,
        3.200, 3.400, 3.600, 3.800, 4.000, 4.200, 4.400, 4.600, 4.800, 5.000, 5.200, 5.400, 5.600, 5.800, 6.000]


class sim_results_t(ctypes.Structure):
    _fields_ = [("fer", ctypes.POINTER(ctypes.c_double)),
                ("ber", ctypes.POINTER(ctypes.c_double)),
                ("avg_iter", ctypes.POINTER(ctypes.c_double)),
                ("time", ctypes.POINTER(ctypes.c_double)),
                ("fec", ctypes.POINTER(ctypes.c_ulong)),
                ("frames", ctypes.POINTER(ctypes.c_ulong))]


def test_sim():
    stop = ctypes.c_ubyte(0)
    res = sim_results_t()


    def sim():
        lib = ctypes.cdll.LoadLibrary("../sim_ldpc.so")
        #lib.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_ubyte, ctypes.c_int)
        lib.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_ubyte, ctypes.c_int, sim_results_t)
        lib.simulate("code.txt".encode("utf-8"), "sim.txt".encode("utf-8"), ctypes.c_int(1), ctypes.byref(stop), ctypes.c_int(0), ctypes.byref(res))
        lib.free_results(ctypes.byref(res))


    thread = threading.Thread(target=sim)
    thread.start()
    print("Thread launched with simulation")
    time.sleep(10)
    print("Stopped after 10 seconds")

    print("Results:")
    for i in range(10):
        print(f"FER: {res.fer[i]} - BER: {res.ber[i]} - Avg_iter: {res.avg_iter[i]} - Time: {res.time[i]}s - Frames: {res.frames[i]}")

    stop.value = 1
    thread.join()

    print("Complete")

def test_rank():
    lib = ctypes.cdll.LoadLibrary("../sim_ldpc.so")
    lib.argtypes = (ctypes.c_char_p, )
    rank = lib.calculate_rank("code.txt".encode("utf-8"))
    print(f"Rank of matrix: {rank}")


if __name__ == "__main__":
    test_rank()