# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 06:52:30 2025

@author: eadsr
"""

from pynput import mouse
import time



clicks = []

def on_click(x, y, button, pressed):
    if pressed:
        print(f"Click recorded at {(x, y)} with {button}")
        clicks.append((x, y))
    else:
        # Stop listener if right button released
        if button.name == "right":
            return False  

time.sleep(5)
with mouse.Listener(on_click=on_click) as listener:
    listener.join()

print("Recorded clicks:", clicks)

# import pyautogui
# import time

# #Update, Export, Next Config
# update_button = (-1858, 513)
# export_button = (-1759, 511)
# force_save_btn = (-908, 591)
# next_config_button = (-953, 494)

# pyautogui.click(update_button)
# time.sleep(3)

# pyautogui.typewrite("qtrace_1.csv")
# pyautogui.press("enter")
# time.sleep(1)
# pyautogui.click(force_save_btn)

# print("qtrace_1.csv")
# time.sleep(2)

# for i in range(2, 17):  # 5 files
#     pyautogui.click(next_config_button)
#     time.sleep(1)

#     pyautogui.click(update_button)
#     time.sleep(3)
    
#     pyautogui.click(export_button)
#     pyautogui.typewrite(f"qtrace_{i}.csv")
#     pyautogui.press("enter")
#     time.sleep(1)
#     pyautogui.click(force_save_btn)
    
#     print(f"qtrace_{i}.csv")
#     time.sleep(2)



