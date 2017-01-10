import numpy as np
import pandas as pd
import kplr
client = kplr.API()

data = pd.read_csv("garcia.txt")
m = data.period.values < 10
print(data.period.values[m])
print(data.KID.values[m])

stars = data.KID.values[m]
for s in stars[3:]:
    print(s)
    star = client.star(s)
    star.get_light_curves(fetch=True)
