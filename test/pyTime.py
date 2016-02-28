import datetime
import sys
import threading

startTime = datetime.datetime.now()
print 'linking charged web now!', startTime

endTime = startTime + datetime.timedelta(0,20,0)
print 'will shut down charging web at: ', endTime

threading._sleep(15)

print 'hello world'
print 'current time: ',datetime.datetime.now()
