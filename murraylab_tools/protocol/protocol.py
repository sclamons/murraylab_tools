import ipywidgets as widgets
from IPython.display import display, Javascript
from traitlets import Unicode, validate
import time
import threading

class Timer():

    def __init__(self,limit=180):
        #self.value = Unicode('00:00:00').tag(sync=True)
        self.labwidg = widgets.Label('00:00:00')
        #display(self.labwidg)
        self.limit = limit
        self.stopTimer = False
        self.thread = None
    def timeit(self, limit=180):
        #display(self)
        hours = 0
        mins = 0
        secs = 0
        for i in range(1,(limit*60+1)):
            if(self.stopTimer):
                self.stopTimer = False
                self.thread = None
                break
            if i%60 == 0:
                if i%3600 == 0:
                    secs = 0
                    mins = 0
                    hours += 1
                else:
                    secs = 0
                    mins += 1
            else:
                secs += 1
            time.sleep(.0167)
            self.labwidg.value = '{hour:02}:{minute:02}:{second:02}'.\
                            format(hour=hours,minute=mins,second=secs)
    #def twrap(self,b,limit=180):
    #    self.timeit(limit=limit)
    def threadTimer(self,b):
        if(b.description == "Start"):
            if(self.thread == None):
                self.thread = threading.Thread(target=self.timeit,args=(self.limit,))
                self.thread.start()
                b.description = "Stop"
        elif(b.description == "Stop"):
            if(self.thread != None):
                self.stopTimer=True
                b.description = "Reset"
        else:
            self.labwidg.value = '00:00:00'
            self.Thread = None
            self.
    def stopTime(self,b):
        self.stopTimer=True

#def startTimer(b,timer):
    #, args=(,))
def lay(width,height):
    return widgets.Layout(width=str(width)+"px",height=str(height)+"px")
def display_timer():
    timer = Timer()
    #self.ddlay = widgets.Layout(width='75px',height='30px')
    button = widgets.Button(description="Start", button_style='info',layout=lay(50,30))
    button2 = widgets.Button(description="Stop", button_style='info',layout=lay(50,30))
    #display(progress)
    #thread.start()
    timerhbox = widgets.HBox([timer.labwidg,button,button2])
    #display(button)
    display(timerhbox)
    #display(timer)
    #timer.twrap(None)
    button.on_click(timer.threadTimer)
    button2.on_click(timer.stopTime)
