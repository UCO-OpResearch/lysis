import cProfile
import os

import lysis


e = lysis.util.Experiment(r'../../data', experiment_code='2022-12-27-1100')
p = {'total_time': 50}
e.initialize_macro_param(p)
macro = lysis.MacroscaleRun(e)
os.makedirs(os.path.join(e.os_path, "macro_pstats"), exist_ok=True)
if os.path.isfile(os.path.join(e.os_path, "macro_pstats", "macro_pstats.sts")):
    last_stats_num = 0
    filename = None
    while True:
        last_stats_num += 1
        filename = "macro_pstats_" + str(last_stats_num).rjust(3, '0') + ".sts"
        if not os.path.isfile(os.path.join(e.os_path, "macro_pstats", filename)):
            break
    os.rename(os.path.join(e.os_path, "macro_pstats", "macro_pstats.sts"),
              os.path.join(e.os_path, "macro_pstats", filename))
cProfile.run("macro.run()", os.path.join(e.os_path, "macro_pstats", "macro_pstats.sts"))
