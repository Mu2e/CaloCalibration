import numpy as np

values_start = range(-70, -40, 5)
values_end = range(-30, 50, 5)
#If a special value has to be included, append it to this temporary list

t_inters = []

for start in values_start:
    for end in values_end:
        if end - start >= 25:
            t_inters.append((start, end))
t_inters = np.array(t_inters)