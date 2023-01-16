import cProfile
import os

from datetime import datetime

import lysis


e = lysis.util.Experiment(r'../../data', experiment_code='2022-12-27-1100')
p = {'total_time': 200}
e.initialize_macro_param(p)
macro = lysis.MacroscaleRun(e)
os.makedirs(os.path.join(e.os_path, "macro_pstats"), exist_ok=True)
now = datetime.now()
timestamp = str.join('', [str(now.year),
                          '_',
                          str(now.month),
                          '_',
                          str(now.day),
                          '_',
                          str(now.hour),
                          str(now.minute),
                          str(now.second)
                          ])
last_stats_num = 0
filename = None
while True:
    last_stats_num += 1
    filename = "macro_pstats_" + str(last_stats_num).rjust(3, '0') + ".sts"
    if not os.path.isfile(os.path.join(e.os_path, "macro_pstats", filename)):
        break
cProfile.run("macro.run()", os.path.join(e.os_path, "macro_pstats", filename))

print(f"cProfile stats saved as {filename}.")
