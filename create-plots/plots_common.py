from pathlib import Path

MARKERS_BY_N = {16 : '+', 24 : 'x', 32: '1', 48 : '*', 64 : '3', 96 : '2', 128 : '.', 'large' : '4'}
MARKERS_BY_L = {2 : '+', 3 : 'x', 4 : '1'}
COLOURS_BY_N = {16 : '#0000a7', 24 : 'orange', 32 : '#008176', 48: '#000000', 64 : '#c1272d', 96 : '#b3b3b3', 128 : 'sienna', 'large': '#424287'}
COLOURS_BY_L = {2 :'#0000a7', 3 : 'orange', 4 : '#008176' }
COLOURS = {'pred' : '#c1272d', 'fit' : '#008176'}
COLOURS_BY_N_T = {16 : '#0000a7', 24 : 'orange', 32 : '#b3b3b3', 'large' : '#008176'}
FIG_SIZE = (10,10 / 1.618)

SAVE_FOLDER = Path("products", "plots")
SAVE_NUMBERS_FOLDER = Path("products", "numbers")