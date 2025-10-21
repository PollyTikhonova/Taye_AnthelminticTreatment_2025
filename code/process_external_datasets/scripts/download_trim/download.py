import pandas as pd
import ENATool

import sys
args = sys.argv

prj_name = str(args[1]).strip()
study_path = str(args[2]).strip()
output_file = str(args[3]).strip()

table = ENATool.fetch(prj_name, path=study_path)
table.to_csv(output_file)
