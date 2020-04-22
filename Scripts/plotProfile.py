from matplotlib import pyplot as plt
#from matplotlib import style
import numpy as np

#style.use('ggplot')

y,x = np.loadtxt('profileJob_26435.csv',
                unpack = True,
                delimiter = ',,')

plt.plot(x,y/1000000)

#plt.xticks(np.arange(min(x), max(x)+1, 25.0))
#plt.yticks(np.arange(318000, 525000, 25000))
plt.title('Memory Usage')
plt.ylabel('Memory Usage (GB)')
plt.xlabel('Time, s')
plt.grid()
plt.show()
plt.savefig('mem_26435.png')
