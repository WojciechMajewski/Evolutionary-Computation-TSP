### TJWP - Travelling Jehovah's Witness Problem

### The dataset contains 3 columns - x, y coordinates and extra weight of the node, with each row representing a node

### Trasa piesek??? PEPE???


import csv
import os
import pygame

dataset_path = ""

with open('TSPA.csv', newline='') as f:
    reader = csv.reader(f, delimiter=";")
    
    data = [[int(row[0]), int(row[1]), int(row[2])] for row in reader if row]

print(data)


from screeninfo import get_monitors
from pathlib import Path

window_size = [800, 600]
window_starting_pos = ((get_monitors()[0].width - window_size[0]) // 2, 60)
target_framerate = 60



background_color = (10, 120, 80)

def windowInit():
    global clock, window#, pixel_font, pixel_font_small

    os.environ['SDL_VIDEO_WINDOW_POS'] = "%i,%i" % window_starting_pos

    pygame.init()
    pygame.mixer.init()
    pygame.font.init()

    #pixel_font = pygame.font.Font('Assets/Gamer.ttf', 96)
    #pixel_font_small = pygame.font.Font('Assets/Gamer.ttf', 48)

    clock = pygame.time.Clock()
    window = pygame.display.set_mode(window_size)

    pygame.display.set_caption('Chimes')

def windowQuit():
    pygame.font.quit()
    pygame.mixer.quit()
    pygame.quit()

def windowBlit():
    global window
    global background_color
    global file_saved_time
    global missing_track

    # -= BACKGROUND =-
    window.fill(background_color)

    # -= NODES =-
    for node in data:
        color_coef = (255 / 2000)
        pygame.draw.circle(window, ((node[2] * color_coef), (node[2] * color_coef), (node[2] * color_coef)), (node[0] * (window_size[0] / 4000), node[1] * (window_size[1] / 2000)), 2)

    pygame.display.flip()

def mainLoop():

    windowInit()

    done = False
    while not done:
        clock.tick(target_framerate)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                done = True

        windowBlit()

    windowQuit()

mainLoop()






