import ctypes
import os
import threading
import numpy as np

LIB_PATH = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}/libldpc.so"

class sim_results_t(ctypes.Structure):
    _fields_ = [("fer", ctypes.POINTER(ctypes.c_double)),
                ("ber", ctypes.POINTER(ctypes.c_double)),
                ("avg_iter", ctypes.POINTER(ctypes.c_double)),
                ("time", ctypes.POINTER(ctypes.c_double)),
                ("fec", ctypes.POINTER(ctypes.c_ulong)),
                ("frames", ctypes.POINTER(ctypes.c_ulong))]

class LDPC:
    def __init__(self, code_file: str, sim_file: str):
        self.code_file = code_file
        self.sim_file = sim_file
        self.sim_results_struct = sim_results_t()
        self.sim_stop_flag = ctypes.c_ubyte(1)

        self.results = {}

        with open(self.sim_file, "r") as simf:
            self.snrs = simf.readlines()[4]
            self.snrs = self.snrs.split(":")[-1]
            self.snrs = [float(x) for x in self.snrs.split(",")]

    def simulate(self, num_threads: int, rng_seed: int):
        def sim_thread():
            lib = ctypes.cdll.LoadLibrary(LIB_PATH)

            lib.argtypes = (sim_results_t, ctypes.c_uint64)
            lib.allocate_results(ctypes.byref(self.sim_results_struct), ctypes.c_uint64(len(self.snrs)))

            # set the flag
            self.sim_stop_flag.value = 0

            lib.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_ubyte, ctypes.c_int, sim_results_t)
            lib.simulate(
                self.code_file.encode("utf-8"), 
                self.sim_file.encode("utf-8"), 
                ctypes.c_int(num_threads), 
                ctypes.byref(self.sim_stop_flag), 
                ctypes.c_int(rng_seed), 
                ctypes.byref(self.sim_results_struct)
            )
            
            #free memory of results manually, so that results dont get lost
            lib.argtypes = (sim_results_t, )
            lib.free_results(ctypes.byref(self.sim_results_struct))
        
        th_sim = threading.Thread(target=sim_thread)
        th_sim.start()

    def stop_simulation(self):
        if not self.sim_stop_flag.value:
            # save results before clearing
            self.results = self.get_results()
            self.sim_stop_flag.value = 1

    def get_results(self):
        if not self.sim_stop_flag.value:
            # we are only interested in entries of frames processed > 0
            max_index = np.sum(np.array(self.sim_results_struct.frames[0:len(self.snrs)]) > 0)

            # convert the cstruct to dictionary
            return dict([(x, getattr(self.sim_results_struct, x)[0:max_index]) for (x,_) in self.sim_results_struct._fields_])
        else:
            return self.results

    def calculate_rank(self):
        lib = ctypes.cdll.LoadLibrary(LIB_PATH)
        lib.argtypes = (ctypes.c_char_p, )
        return lib.calculate_rank(self.code_file.encode("utf-8"))
