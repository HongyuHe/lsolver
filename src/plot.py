#!/usr/local/bin/python3 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy import interpolate
from scipy import optimize 
from scipy import integrate
from datetime import datetime
from datetime import timedelta
import functools as fct
import matplotlib
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
from typing import *
import math
import warnings
import random
warnings.filterwarnings('ignore')
plt.style.use('default')

def plot_stacked_bar(data, series_labels, category_labels=None, 
                     show_values=False, value_format="{}", 
                     y_label=None, x_label=None, title=None, 
                     colors=None, grid=True, reverse=False):
    """Plots a stacked bar chart with the data and labels provided.

    Keyword arguments:
    data            -- 2-dimensional numpy array or nested list
                       containing data for each series in rows
    series_labels   -- list of series labels (these appear in
                       the legend)
    category_labels -- list of category labels (these appear
                       on the x-axis)
    show_values     -- If True then numeric value labels will 
                       be shown on each bar
    value_format    -- Format string for numeric value labels
                       (default is "{}")
    y_label         -- Label for y-axis (str)
    x_label         -- Label for x-axis (str)
    title           -- Title of the plot (str)
    colors          -- List of color labels
    grid            -- If True display grid
    reverse         -- If True reverse the order that the
                       series are displayed (left-to-right
                       or right-to-left)
    """

    ny = len(data[0])
    ind = list(range(ny))

    axes = []
    cum_size = np.zeros(ny)

    data = np.array(data)

    if reverse:
        data = np.flip(data, axis=1)
        category_labels = reversed(category_labels)

    for i, row_data in enumerate(data):
        color = colors[i] if colors is not None else None
        axes.append(plt.bar(ind, row_data, bottom=cum_size, 
                            label=series_labels[i], color=color))
        cum_size += row_data

    if category_labels:
        plt.xticks(ind, category_labels)

    if y_label:
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        plt.title (title)

    plt.legend()

    if grid:
        plt.grid()

    if show_values:
        for axis in axes:
            for bar in axis:
                w, h = bar.get_width(), bar.get_height()
                plt.text(bar.get_x() + w/2, bar.get_y() + h/2, 
                         value_format.format(h), ha="center", 
                         va="center")

if __name__ == '__main__':
    plt.figure(figsize=(6, 4))

    series_labels = ['naive_par', 'naive']
    category_labels = [str(i) for i in range(1,9)]

    # test
    data = [
        [0.1, 0.3, 0.35, 0.3],
        [0.8, 0.7, 0.6, 0.5]
    ]

    plot_stacked_bar(
        data, 
        series_labels, 
        category_labels=category_labels, 
        show_values=True, 
        value_format="{:.3f}",
        colors=['tab:blue', 'tab:orange'],
        y_label="Execution Time (ms)",
        x_label="#Threads",
        title="Matrix: torso1"
    )
    plt.show()