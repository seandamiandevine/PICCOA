# -*- coding: utf-8 -*-
"""
code to display ethical scenarios
author: David Levari, david.levari@gmail.com
edited: Sean Devine, seandamiandevine@gmail.com
last updated: December 4, 2018

based on studies found in:
Levari, D. E., Gilbert, D. T., Wilson, T. D., Sievers, B., Amodio, D. M.,
& Wheatley, T. (2018). Prevalence-induced concept change in human judgment.
Science, 360(6396), 1465-1467.

Parameters:
    id: Participant ID 
    age: Age of the participant 
    cn: condition (0 = stable prevalence; 1 = decreasing prevalence)
    cb: counterbalance (1 = dots task first; 2 = ethics task first)
    inBattery: is this task in a Battery? 1 for yes, in which case quit at the end. 0 for no 
    debug: use for testing task; quits early 
"""

from psychopy import visual, core, gui, event, monitors
from psychopy.visual import ShapeStim, ImageStim
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np
import pandas as pd
import random
import csv
import os

#Set directory
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

#Define functions
def addOutput(file, output):
    with open(file, 'a') as data:
        writer = csv.writer(data)
        writer.writerow(output)

def runEthics(id, age, cn, cb, inBattery, debug): 
    def stimRandomizer(trials, thisfreq): #need to confirm proportions
        thisAmbig = random.sample(ambig, int(round(trials*.33)))
        signal = random.sample(awful, int(round(trials*thisfreq)))
        noise = random.sample(harmless, int(round(trials*(1-thisfreq-.33))))
        
        thisBlockStim = signal+noise+thisAmbig
        random.shuffle(thisBlockStim)
        return(thisBlockStim)

    #files to be used
    dir = "stimuli/"
    Awful = []
    Ambig = []
    Harmless = []
    with open('normed_ethics.csv') as vignettes: 
        vignettes = csv.reader(vignettes, delimiter = ',')
        line = 0
        for row in vignettes:
            if line > 0: 
                if row[9] == 'awful': Awful.append(row[1])
                elif row[9] == 'meh': Ambig.append(row[1])
                elif row[9] == 'harmless':Harmless.append(row[1])
            line += 1
    awful = list(set(Awful))
    ambig = list(set(Ambig))
    harmless = list(set(Harmless))
    
    #setup datafile
    filename = _thisDir + os.sep + 'data/PICCethics_' + id+'.csv'
    with open(filename, 'w') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerow(['id', 'age', 'condition', 'counterbalance', 'block', 'trialinblock', 'trial', 'tstart', 'tend', 'trialstim', 'response', 'RT', 'trialtype', 'ethicfreq', 'norm_mean', 'norm_sd', 'norm_se'])
    
    #Window setup
    win = visual.Window(
        size=[1920, 1080], fullscr=True, screen=0,
        allowGUI=False, allowStencil=False,
        monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
        blendMode='avg', useFBO=True,
        units='cm')
    
    #set constants 
    numtrials = 24
    numblocks = 10 if not debug else 1
    
    advance = ['space']
    signalKey = "a"
    noiseKey = "l"
    acceptedKeys = ["a","l"]
    
    if cn == '0': 
        ethicsfreq = [.333]*numblocks
    if cn == '1': 
        ethicsfreq = [.333,.333,.333,.333,.25,.166,.083,.0412,.0412,.0412]
    
    trials = []
    for t in range(numblocks): 
        tmp = stimRandomizer(numtrials, ethicsfreq[t])
        #make sure there are are no repeats
        harmless = [x for x in harmless if x not in tmp]
        ambig = [x for x in ambig if x not in tmp]
        awful = [x for x in awful if x not in tmp]
        trials.append(tmp)
    
    #initialize instructions
    instWait = 3 #time to wait before being able to progress instructions
    
    pressSpace = visual.TextStim(win=win, name='i1',
        text="Press SPACE to continue",
        font='Arial',
        pos=(0, -10), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i1 = visual.TextStim(win=win, name='i1',
        text="Welcome to this task! We are interested in studying how people \
make ethical decisions about scientific experiments.",
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
        
    i2 = visual.TextStim(win=win, name='i2',
        text="Many scientific experiments involve some risk for the participants because they can cause psychological distress or physical harm. Universities have to make difficult ethical decisions about whether or not to allow experiments to be conducted.",
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
        
    i3 = visual.TextStim(win=win, name='i3',
        text='Today, you will read about various experiments that could be conducted on human beings. We simply want to know whether you think scientists SHOULD or SHOULD NOT be allowed to conduct each of these experiments.',
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i4 = visual.TextStim(win=win, name='i4',
        text='Because this is an ethical decision, there are no right or wrong answers. We simply want your personal decision for each study.',
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i5 = visual.TextStim(win=win, name='i5',
        text='Here are some things to keep in mind as you make your decisions.\n\n\
1) All of the experiments you will read about will be conducted on adults who have volunteered to take part in exchange for money.\
    \n\n\
2) All of the experiments are part of research on human behavior.\
    \n\n\
3) When scientists must lie to the participants either before or during the experiment, they always tell the participants the truth when the experiment is over.\
    \n\n\
4) Participants are always free to withdraw and can stop participating at any time they wish.',
        font='Arial',
        pos=(0, 0), height=0.9, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
        
    i6 = visual.TextStim(win=win, name='i6',
        text='In the task, you will see descriptions of experiments presented on the screen, one at a time.',
        font='Arial',
        pos=[0, 0], height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i7 = visual.TextStim(win=win, name='i7',
        text='When you read a description of an experiment that you would not allow to be conducted, \
press the "A" key. For all other experiments, press the "L" key.',
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i8 = visual.TextStim(win=win, name='i8',
        text="The experiments will be presented in series, with breaks in between. \
This means that you will read a series of experiments, have a short break, \
and then another series of experiments, until you have seen 10 series.",
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i9 = visual.TextStim(win=win, name='i9',
        text="Some of the series you see may have a lot of unethical experiments, and others may have only a few. There's nothing for you to count or keep track of -- your only task is to approve or reject each experiment.",
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i10 = visual.TextStim(win=win, name='i10',
        text='You will have as much time to answer as you need. That being said, you should do your best to answer quickly and accurately during the study. If you make a mistake and hit the wrong button at any point, just keep going.',
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
        
    i11 = visual.TextStim(win=win, name='i11',
        text='Now you will complete a brief practice round so you can get used to the task.',
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    i12 = visual.TextStim(win, text="You have now completed the practice series. If you have any questions, you can ask the experimenter now. Otherwise, you're ready to begin the study.",
        height=1, color='white', pos=(0, 0), wrapWidth = 30)
    
    instsA = [i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11]
    
    #initialize trial components
    trialClock = core.Clock()
    vignette = visual.ImageStim(
        win=win, name='vignette',
        image='stimuli/prac_217_1.png', mask=None,
        ori=0, pos=(0,0), size=(1250, 750), units = 'pix', 
        color=[1,1,1], colorSpace='rgb', opacity=1,
        flipHoriz=False, flipVert=False,
        texRes=128, interpolate=True, depth=0.0)
    
    fixTime = 0.5 #time that the fixation across appears after a choice
    fix = visual.TextStim(win=win, name='fix',
        text='+',
        font='Arial',
        pos=(0, 0), height=2, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    blockTxt = visual.TextStim(win=win, name='blockTxt',
        text='',
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    endBlocktxt = visual.TextStim(win=win, name='endBlocktxt',
        text="Series complete!\n\nPlease take a short break. We'll start the next series in a moment.",
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    
    end = visual.TextStim(win=win, name='end',
        text="You are finished with this task! Please notify the tester.",
        font='Arial',
        pos=(0, 0), height=1, wrapWidth=30, ori=0,
        color='white', colorSpace='rgb');
    #--------------------------------------------------------Start task-----------------------------------------------------------------------
    #instructions 
    for i in instsA: 
        i.draw()
        win.flip()
        core.wait(instWait) 
        i.draw()
        pressSpace.draw()
        win.flip()
        event.waitKeys(keyList = advance)
        win.flip()
    
    #practice 
    trialClock.reset()
    tstart = trialClock.getTime()
    vignette.draw()
    win.flip()
    RTtime = trialClock.getTime()
    response = event.waitKeys(keyList = acceptedKeys)
    tend = trialClock.getTime()
    RT = tend - RTtime
    addOutput(filename, [id, age, cn, cb, 0, 0, 'practice', tstart, tend, vignette.image, response[0], RT, 'harmless', 'NA', 'NA', 'NA', 'NA'])
    fix.draw()
    win.flip()
    core.wait(fixTime)
    
    #end practice
    i12.draw()
    win.flip()
    core.wait(instWait)
    i12.draw()
    pressSpace.draw()
    win.flip()
    event.waitKeys(keyList = advance)
    win.flip()
    
    #trials
    totalTrial = 0
    trialClock.reset()
    for b in range(numblocks):
        blockTxt.text = 'Block {}'.format(str(b+1))
        blockTxt.draw()
        win.flip()
        core.wait(instWait)
        blockTxt.draw()
        pressSpace.draw()
        win.flip()
        event.waitKeys(keyList = advance)
        win.flip()
        core.wait(fixTime)
        for t in range(numtrials):
            totalTrial += 1
            vignette.image = 'stimuli/' + trials[b][t] + '.png'
            tstart = trialClock.getTime()
            vignette.draw()
            win.flip()
            RTtime = trialClock.getTime()
            response = event.waitKeys(keyList = acceptedKeys)
            tend = trialClock.getTime()
            if trials[b][t] in Awful: trialtype = 'awful'
            elif trials[b][t] in Ambig: trialtype = 'ambig'
            else: trialtype = 'harmless' 
            RT = tend - RTtime
            # add extra norm info
            with open('normed_ethics.csv') as vignettes: 
                vignettes = csv.reader(vignettes, delimiter = ',')
                line = 0
                for row in vignettes:
                    if line > 0: 
                       if row[1] == trials[b][t]: 
                        norm_mean = row[11]
                        norm_sd = row[12]
                        norm_se = row[13]
                        break 
                    line += 1
            addOutput(filename, [id, age, cn, cb, b+1, t+1, totalTrial, tstart, tend, vignette.image, response[0], RT, trialtype, ethicsfreq[b], norm_mean, norm_sd, norm_se])
            fix.draw()
            win.flip()
            core.wait(fixTime)
            win.flip()
        if b != 9: 
            endBlocktxt.draw()
            win.flip()
            core.wait(4)
    
    #show end screen
    if inBattery: 
        end.text = 'You are finished with this task! \n\nPress SPACE when you are ready to move on to the next task.'
    end.draw()
    win.flip()
    event.waitKeys(keyList = advance)
    win.close()