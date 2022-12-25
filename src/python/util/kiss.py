import ctypes
import os


class KissGenerator(object):
    def __init__(self, seed=None, state=None):
        path = os.path.dirname(__file__)
        lib_path = os.path.join(path, '..', '..', '..', 'lib')
        kiss_file = 'kiss.so'

        my_kiss = ctypes.CDLL(os.path.join(lib_path, kiss_file))

        self.state_type = ctypes.c_uint * 4

        self.urcw1 = my_kiss.urcw1_
        self.urcw1.restype = ctypes.c_double

        self.mscw = my_kiss.mscw_
        self.mscw.restype = ctypes.c_int

        self.kiss32 = my_kiss.kiss32_
        self.kiss32.restype = ctypes.c_uint32

        self.__set_kiss32 = my_kiss.set_kiss32_
        self.__set_kiss32.argtypes = [self.state_type]

        self.__get_kiss32 = my_kiss.get_kiss32_
        self.__get_kiss32.argtypes = [self.state_type]

        if state is not None:
            self.set_state(state)
        elif seed is not None:
            self.set_seed(seed)
        else:
            self.set_seed(self.mscw())

        self.urcw1 = my_kiss.urcw1_
        self.urcw1.restype = ctypes.c_double
        my_kiss.urcw1_()
        
    def set_state(self, state):
        c_state = self.state_type(state[0], state[1], state[2], state[3])
        self.__set_kiss32(c_state)

    def set_seed(self, seed):
        state = self.get_state()
        state[3] = seed
        self.set_state(state)

    def get_state(self):
        c_state = self.state_type(0, 0, 0, 0)
        self.__get_kiss32(c_state)
        return [int(i) for i in c_state]
