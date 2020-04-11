import ipywidgets as widgets
from IPython.display import display, Javascript
from traitlets import Unicode, validate
import time
import threading

class Timer():

    def __init__(self,limit=180):
        #self.value = Unicode('00:00:00').tag(sync=True)
        self.labwidg = widgets.Button(description="00:00:00",layout=lay(120,50))
        self.labwidg.on_click(self.threadTimer)
        #self.labwidg.style.button_color = "gree"
        #widgets.Label('00:00:00')
        #display(self.labwidg)
        self.limit = limit
        self.stopTimer = False
        self.thread = None
        self.buttonState = "Start"
    def timeit(self, limit=180):
        #display(self)
        hours = 0
        mins = 0
        secs = 0
        for i in range(1,(limit*60+1)):
            if(self.buttonState == "Start"):
                self.labwidg.button_style="primary"
                #self.labwidg.style.font_weight = "50px"
            elif(self.buttonState == "Stop"):
                self.labwidg.button_style="success"
                #self.labwidg.style.font_weight = "50px"
            elif(self.buttonState == "Reset"):
                self.labwidg.button_style="warning"
                #self.labwidg.style.font_weight = "50px"
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
            self.labwidg.description = '{hour:02}:{minute:02}:{second:02}'.\
                            format(hour=hours,minute=mins,second=secs)
    #def twrap(self,b,limit=180):
    #    self.timeit(limit=limit)
    def threadTimer(self,b):
        #b.description = self.buttonState
        if(self.buttonState == "Start"):
            if(self.thread == None):
                self.thread = threading.Thread(target=self.timeit,args=(self.limit,))
                self.thread.start()
                self.buttonState = "Stop"
        elif(self.buttonState == "Stop"):
            if(self.thread != None):
                self.stopTimer=True
                self.buttonState = "Reset"
        else:
            self.labwidg.description = '00:00:00'
            self.Thread = None
            self.stopTimer=False
            self.buttonState = "Start"
    def stopTime(self,b):
        self.stopTimer=True

#def startTimer(b,timer):
    #, args=(,))
def lay(width,height):
    return widgets.Layout(width=str(width)+"px",height=str(height)+"px")
def display_timer():
    timer = Timer()
    display(timer.labwidg)
