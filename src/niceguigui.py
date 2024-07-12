##
# @internal
# @file gui.py
# @brief Create the graphical user interfaces. 
# There are 3 options:
# -# single_mode
# -# multi_mode
# -# batch_mode

import info
# from process import *
# from plots import *
# from visualisation import *
import os
import oineus
from nicegui import ui

def entry_window():
	# ui.label("Topological Amorphous Material Analysis")
	with ui.tabs().classes('w-full') as tabs:
		one = ui.tab('One')
		two = ui.tab('Two')
	with ui.tab_panels(tabs, value=two).classes('w-full'):
		with ui.tab_panel(one):
			ui.label('First tab')
			# ui.restructured_text(info.license())
		with ui.tab_panel(two):
			ui.label('Second tab')

	ui.run()