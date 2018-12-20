import csv
import matplotlib.pyplot as mpl

t=[]
dv=[]
udv=[]

with open('parsedfile_velocity.csv') as csvfile:
    header=csv.reader(csvfile,delimiter=',')

    for row in header:
        t.append(float(row[0]))
        dv.append(float(row[1]))
        udv.append(float(row[2]))



fig, ax = mpl.subplots()
ax.plot(t, dv, 'r', label='damped velocity')
ax.plot(t, udv, 'b', label='undamped velocity')

ax.set(xlabel='time (s)', ylabel='velocity (degrees/s)',
       title='comparison for damped and undaped motion')
ax.grid()

mpl.show()
