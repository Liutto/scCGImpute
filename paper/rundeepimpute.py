from deepimpute.multinet import MultiNet
import pandas as pd
import time
from playsound import playsound
start_time = time.time()
data = pd.read_csv('D:/software/RStudio/single_cell/datasets/simulate_splat/simulate_splat.csv', index_col=0) # dimension = (cells x genes)
data = data.transpose()
print('Working on {} cells and {} genes'.format(*data.shape))
model = MultiNet()
model.fit(data)
imputed = model.predict(data)
print('Working on {} cells and {} genes'.format(*imputed.shape))
print("--- %s seconds ---" % (time.time() - start_time))
pd.DataFrame.to_csv(imputed.transpose(), "D:/software/RStudio/single_cell/datasets/simulate_splat/simulate_splat_deepimpute.csv")

playsound('E:/music/未闻花香.m4a')