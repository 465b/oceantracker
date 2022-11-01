import numpy as np
import matplotlib.pyplot as plt
import os
import json



bio_labels = ["no modifiers","dryout when stranded","light starvation","splitting & culling","sinking","diel migration"]

path = "/scratch/local1/output/22_07_12_comp_cost_of_bio_v00.py"
files = os.listdir(path)
caseInfos = [item for item in files if 'caseInfo.json' in item]

bio_cases = {}
for item in caseInfos:
    case_file = os.path.join(path,item)

    bio_case = item.split('_')[-2][1:4]
    with open(case_file, 'r') as json_file:
        caseInfo = json.load(json_file)
    model_duration = caseInfo['run_info']['model_run_duration']
    model_duration = int(model_duration.split(':')[1])*60+float(model_duration.split(':')[2]) 
    if bio_case not in bio_cases:
        bio_cases[bio_case] = [model_duration]
    else:
        bio_cases[bio_case].append(model_duration)

mean = [np.mean(bio_cases[key]) for key in bio_cases]
std = [np.std(bio_cases[key]) for key in bio_cases]

plt.bar(x=bio_labels,height=mean,yerr=std)
plt.xticks(rotation=30)
plt.hlines(mean[0],xmin=0,xmax=5)
plt.tight_layout()
plt.savefig('bio_comp_cost.png')

cost = [item-mean[0] for item in mean]
cost[0] = mean[0]

fig, ax = plt.subplots()
# First plot the 'Male' bars for every day.
for ii,item in enumerate(cost):
    if ii == 0:
        ax.bar(x=1, height=item,label=bio_labels[ii])
    else:
        bottom = np.sum(cost[:ii])
        print(bottom)
        ax.bar(x=1, height=item,label=bio_labels[ii],bottom=bottom)
plt.legend()
plt.savefig('bio_summed_cost.png')



print(cost)


    