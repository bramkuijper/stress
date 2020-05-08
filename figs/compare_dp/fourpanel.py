import stress_panels
import pandas as pd
import numpy as np


dp_summary_file = "summary_dp.csv"
sims_summary_file = "summary_sims.csv"


translate_cols = {
    "t":"time"
    ,"pLeave":"sP2NP_1"
    ,"pArrive":"sNP2P_1"
    ,"pAttack":"p_att"
    ,"hormone":"mean_hormone"
    }


# translate column values
dp_data = dp_data.rename(
    columns=translate_cols
    )

dp_data["aP"] = 1.0
dp_data["ad"] = 1.5
