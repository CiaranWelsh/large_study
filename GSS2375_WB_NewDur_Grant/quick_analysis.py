import pandas


f = r'C:\Users\Ciaran\Documents\LargeStudy\GSS2375_WB_NewDur_Grant\GSS2375_RNA.xlsx'

data = pandas.read_excel(f)

import matplotlib.pyplot as plt
plt.figure()
plt.hist(data['Yield(ug)'],bins=50)
plt.show()