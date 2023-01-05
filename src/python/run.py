import cProfile
import lysis


e = lysis.util.Experiment(r'../../data', experiment_code='2022-12-27-1100')
p = {'total_time': 10}
e.initialize_macro_param(p)
macro = lysis.MacroscaleRun(e)
cProfile.run("macro.run()")
