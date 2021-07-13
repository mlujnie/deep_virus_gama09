from gama_tools import *

for gama_field in ["gama09E", "gama09F"][1:]:
    save_stars(gama_field, overwrite=False)
    print("Done with ", gama_field)
