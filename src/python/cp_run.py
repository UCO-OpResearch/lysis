import cProfile
import os

import lysis


e = lysis.util.Experiment(r'../../data', experiment_code='2022-12-27-1100')
p = {'total_time': 1}
e.initialize_macro_param(p)
macro = lysis.CudaMacroscaleRun(e)

macro.run()
