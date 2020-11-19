import ctypes as ct
import os
import threading
import numpy as np

LIB_PATH = f"{os.path.dirname(os.path.abspath(__file__))}/libldpc.so"

class sim_results_t(ct.Structure):
    _fields_ = [("fer", ct.POINTER(ct.c_double)),
                ("ber", ct.POINTER(ct.c_double)),
                ("avg_iter", ct.POINTER(ct.c_double)),
                ("time", ct.POINTER(ct.c_double)),
                ("fec", ct.POINTER(ct.c_uint64)),
                ("frames", ct.POINTER(ct.c_uint64))]

class decoder_param(ct.Structure):
    _fields_ = [("earlyTerm", ct.c_bool),
                ("iterations", ct.c_uint32),
                ("type", ct.c_char_p)]

class channel_param(ct.Structure):
    _fields_ = [("seed", ct.c_uint64),
                ("double", ct.c_double * 3),
                ("type", ct.c_char_p)]

class simulation_param(ct.Structure):
    _fields_ = [("threads", ct.c_uint32),
                ("maxFrames", ct.c_uint64),
                ("fec", ct.c_uint64),
                ("resultFile", ct.c_char_p)]

class LDPC:
    def __init__(self, pc_file: str, gen_file = "", lib = LIB_PATH):    
        self.pc_file = pc_file
        self.gen_file = gen_file 
        self.n = ct.c_int(0)
        self.m = ct.c_int(0)

        self.lib = ct.cdll.LoadLibrary(lib)
        self.lib.argtypes = (ct.c_char_p, ct.c_char_p, ct.c_int, ct.c_int)
        self.lib.ldpc_setup(pc_file.encode("utf-8"), gen_file.encode("utf-8"), ct.byref(self.n), ct.byref(self.m))

        self.n = self.n.value
        self.m = self.m.value
        self.k = self.n - self.m

        self.sim_stop_flag = ct.c_bool(False)
        vec_double = ct.c_double * 50
        vec_u64 = ct.c_uint64 * 50
        self.sim_results_struct = sim_results_t(
            vec_double(),
            vec_double(),
            vec_double(),
            vec_double(),
            vec_u64(),
            vec_u64(),
        )
        self.results = {}
        self.sim_params = {
            "earlyTerm": True,
            "iterations": 50,
            "decoding": "BP",
            "seed": 0,
            "snr": [],
            "channel": "AWGN",
            "threads": 1,
            "maxFrames": int(10e9),
            "fec": 50
        }


    def encode(self, info_word: np.array) -> np.array:
        """Encode a binary array.

        Args:
            info_word (np.array): Input binary array.

        Raises:
            RuntimeError: No generator matrix is initially provided.

        Returns:
            np.array: Encoded binary codeword.
        """
        if not self.gen_file:
            raise RuntimeError("No generator matrix provided for encoding")

        vec_in = ct.c_uint8 * self.k
        vec_out = ct.c_uint8 * self.n
        in_arr = vec_in(*info_word)
        out_arr = vec_out()

        self.lib.argtypes = (vec_in, vec_out)
        self.lib.encode(ct.byref(in_arr), ct.byref(out_arr))

        return np.array(out_arr[0:self.n])



    def decode(self, llr_in: np.array, early_term=True, iters=50, dec_type="BP") -> np.array:
        """Decode array of input LLRs.

        Args:
            llr_in (np.array): Input LLR, length n (transmitted)
            early_term (bool, optional): Terminate decoding if codeword 
            is valid. Defaults to True.
            iters (int, optional): Number of iterations. Defaults to 50.
            dec_type (str, optional): Type of decoding. See libldpc documentation
            Defaults to "BP".

        Returns:
            np.array: Output LLR, length n (transmitted)
        """
        dec_params = decoder_param(early_term, iters, dec_type.encode("utf-8"))

        vec_double = ct.c_double * self.n
        in_arr = vec_double(*llr_in)
        out_arr = vec_double()

        self.lib.argtypes = (decoder_param, vec_double, vec_double)
        self.lib.restype = ct.c_int
        iter_req = self.lib.decode(dec_params, ct.byref(in_arr), ct.byref(out_arr))

        return np.array(out_arr[0:self.n]), iter_req


    def simulate(self, **args):
        """Start the simulation in threaded mode.

        Args (optional):
            earlyTerm (bool): Terminate decoding if codeword valid
            iterations (int): Number of decoding iterations
            decoding (str): "BP", "BP_MS"
            seed (int): RNG Seed
            snr (list): [MIN, MAX, STEP]
            channel (str): "AWGN", "BSC", "BEC"
            threads (int): Number of parallel threads
            maxFrames (int): Maximum number of frames per SNR value
            fec (int): Number of error frames per SNR value
        """
        snr = ct.c_double * 3
        self.sim_params = {**self.sim_params, **args}
        snr = snr(*self.sim_params["snr"])
        dec_param = decoder_param(self.sim_params["earlyTerm"], self.sim_params["iterations"], self.sim_params["decoding"].encode("utf-8"))
        ch_param = channel_param(self.sim_params["seed"], snr, self.sim_params["channel"].encode("utf-8"))
        sim_param = simulation_param(self.sim_params["threads"], self.sim_params["maxFrames"], self.sim_params["fec"], "".encode("utf-8"))

        def sim_thread():
            self.sim_stop_flag.value = False

            self.lib.argtypes = (decoder_param, channel_param, simulation_param, sim_results_t, ct.c_bool)
            self.lib.simulate(
                dec_param,
                ch_param,
                sim_param,                
                ct.byref(self.sim_results_struct),
                ct.byref(self.sim_stop_flag)
            )
        
        th_sim = threading.Thread(target=sim_thread)
        th_sim.start()

    def stop_simulation(self):
        """Stop the simulation.
        """
        if not self.sim_stop_flag.value:
            # save results before clearing
            self.results = self.get_results()
            self.sim_stop_flag.value = 1

    def get_results(self):
        """Return the simulation results.

        Returns:
            dict: Simulation results
        """
        if not self.sim_stop_flag.value:
            # we are only interested in entries of frames processed > 0
            max_index = np.sum(np.array(self.sim_results_struct.frames[0:50]) > 0)

            # convert the cstruct to dictionary
            return dict([(x, getattr(self.sim_results_struct, x)[0:max_index]) for (x,_) in self.sim_results_struct._fields_])
        else:
            return self.results

    def rank(self):
        """Calculate the rank of the parity-check matrix over GF(2).

        Returns:
            int: Rank of H
        """
        return self.lib.calculate_rank()
