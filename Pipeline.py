# Python script to take variables and write Variables.R

import tkinter as tk

class Application(tk.frame):
	def __init__(self, master=none):
		tk.Frame.__init__(self, master)
		self.grid()
		self.createWidgets()
		
	def createWidgets(self):
		self.quitButton = tk.Button(self, text= 'Quit1', command=self.quit)
		self.quitButton.grid()

app = Application()
app.master.title('Sample Application')
app.mainloop()
