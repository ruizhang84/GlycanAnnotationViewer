import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon

def plot_circle(center_x, center_y, color="y", r=0.5):
    circle = plt.Circle((center_x, center_y), r, color=color)
    return circle
def plot_rectangle(center_x, center_y, color="b", r=0.5):
    xy = [[center_x - r, center_y - r], [center_x - r, center_y + r], 
          [center_x + r, center_y + r], [center_x + r, center_y - r]]
    polygon = Polygon(xy, True, color=color)
    return polygon
def plot_triangle(center_x, center_y, orient="left", color="r", r=0.5):
    if orient == "left":
        xy = [[center_x, center_y - r], [center_x, center_y + r], [center_x - r * np.sqrt(3), center_y]]
    else:
        xy = [[center_x, center_y - r], [center_x, center_y + r], [center_x + r * np.sqrt(3), center_y]]
    polygon = Polygon(xy, True, color=color)
    return polygon
def plot_diamond(center_x, center_y, color="m", r=0.5):
    center_y -= r
    r *= 2
    xy = [[center_x, center_y], [center_x + r * np.sqrt(2) / 2, center_y +  r * np.sqrt(2) / 2], 
      [center_x, center_y +  r * np.sqrt(2)], [center_x -  r * np.sqrt(2) / 2, center_y +  r * np.sqrt(2) / 2]]
    polygon = Polygon(xy, True, color=color)
    return polygon

def compute_antenna_y(y, gap, r=0.5, d=0.5):
    """
    y: the center of bottom circle.
    d: the distance to connect two circles.
    gap: the distance between two circles.
    """
    hypotenuse = 2.0 * r + d
    return y + np.sqrt(hypotenuse * hypotenuse - gap * gap / 4.0)

def plot_antenna_line(x, y, gap, orient="left", color="k", width=2, r=0.5, d=0.5):
    """
    r: radius.
    d: bond length.
    """
    hypotenuse = 2.0 * r + d
    y_distance = np.sqrt(hypotenuse * hypotenuse - gap * gap / 4.0)
    x1 = x2 = None
    if orient == "left":
        x1 = x - (gap / 2.0) / hypotenuse * r
        x2 = x - (gap / 2.0) / hypotenuse * (r + d)
    else:
        x1 = x + (gap / 2.0) / hypotenuse * r 
        x2 = x + (gap / 2.0) / hypotenuse * (r + d)
    y1 = y_distance / hypotenuse * r + y
    y2 = y_distance / hypotenuse * (r + d) + y
    plot_line(x1, y1, x2, y2, color, width)
    
def plot_line(x1, y1, x2, y2, color="black", width=2):
    plt.plot([x1, x2], [y1, y2], color=color, linewidth=width)

def compute_gap(d, r, angle=45):
    """
    d = distance between layers
    r = radius betwen circles
    """
    degree = angle / 180.0 * np.pi
    return (d + 2 * r) * np.tan(degree) * 2

def compute_distance(gap, r, angle=45):
    degree = angle / 180.0 * np.pi 
    return gap / np.sin(degree) / 2.0 - 2 * r

def select_orient(i):
    if i%2 == 0:
        return "left"
    return "right"
def plot_antenna_rect_line(x, y, gap, orient="left", color="k", width=2, r=0.5, d=0.5):
    """
    r: radius.
    d: bond length.
    """
    hypotenuse = 2.0 * r + d
    y_distance = np.sqrt(hypotenuse * hypotenuse - gap * gap / 4.0)
    x1 = x2 = None
    if orient == "left":
        x1 = x - (gap / 2.0) / hypotenuse * r
        x2 = x - (gap / 2.0) / hypotenuse * (r + d)
    else:
        x1 = x + (gap / 2.0) / hypotenuse * r 
        x2 = x + (gap / 2.0) / hypotenuse * (r + d)
    y1 = y_distance / hypotenuse * r + y
    y2 = y_distance / hypotenuse * (r + d + r - r * np.sqrt(2)) + y
    plot_line(x1, y1, x2, y2, color, width)
def branch_full(table, i):
    return table[6 + i * 2] > 0 and table[6 + i * 2 + 1] > 0

def glycan_complex_plot(start_x, start_y, line_length, ax, table):
    if len(table) != 26:
        return 
    
    draw_table = table[:]
    
    x, y = start_x, start_y
    d = line_length
    
    # fucose core
    if draw_table[2] > 0:
        obj = plot_triangle(x - 2 * d, y, orient="left", color="r", r=d)
        ax.add_patch(obj)
        plot_line(x - d, y, x - 2 * d, y, color="k", width=2)
        
    # GlcNAc
    while draw_table[0] > 0:
        obj = plot_rectangle(x, y, color="b", r=d)
        ax.add_patch(obj)
        plot_line(x, y + d, x, y + 2 * d, color="k", width=2)
        # update
        y += 3 * d
        draw_table[0] -= 1
    
    # Man
    if draw_table[1] > 0:
        obj = plot_circle(x, y, color="g", r=d)
        ax.add_patch(obj)      
    
    # bisect:
    if draw_table[3] > 0:
        plot_line(x, y + d, x, y + 2 * d, color="k", width=2)
        obj = plot_rectangle(x, y + 3 * d, color="b", r=d)
        ax.add_patch(obj)

        
    # antenna
    gap = compute_gap(d=d/3.0, r=d, angle=61)
    for i in range(2):
        x_bar = x - (-1) ** i * gap / 2.0
        y_bar = compute_antenna_y(y, gap=gap, r=d, d=compute_distance(gap=gap, r=d, angle=61))
        
        # Man
        if table[4 + i] > 0:
            obj = plot_circle(x_bar, y_bar, color="g", r=d)
            ax.add_patch(obj)
            plot_antenna_line(x, y, gap=gap, orient=select_orient(i), color="k", width=2, r=d, 
                              d=compute_distance(gap=gap, r=d, angle=61))
        
        # chain
        gap_bar = compute_gap(d=d, r=d)
        for j in range(2):
            if (draw_table[6 + i * 2 + j] + draw_table[10 + i * 2 + j]
             + draw_table[14 + i * 2 + j] + draw_table[18 + i * 2 + j]) == 0:
                continue
            
            x_bar_bar = x_bar
            y_bar_bar = y_bar
            
            if branch_full(table, i):
                x_bar_bar = x_bar - (-1) ** j * gap_bar / 2.0
                y_bar_bar = compute_antenna_y(y_bar, gap=gap_bar, r=d, d=compute_distance(gap_bar, d))
            
            # GlcNAc
            if branch_full(table, i):
                obj = plot_rectangle(x_bar_bar, y_bar_bar, color="b", r=d)
                ax.add_patch(obj)
                plot_antenna_rect_line(x_bar, y_bar, gap=gap_bar, 
                                       orient=select_orient(j), color="k", width=2, r=d,
                                       d=compute_distance(gap_bar, d))
                draw_table[6 + i * 2 + j] -= 1
            elif draw_table[6 + i * 2 + j] > 0:
                plot_line(x_bar_bar, y_bar_bar + d, x_bar_bar, y_bar_bar + 2 * d, color="k", width=2)
                obj = plot_rectangle(x_bar_bar, y_bar_bar + 3 * d, color="b", r=d)
                ax.add_patch(obj)
                draw_table[6 + i * 2 + j] -= 1
                y_bar_bar += 3 * d
                
                
            # chain enlongation
            x_bar_bar_bar = x_bar_bar
            y_bar_bar_bar = y_bar_bar

            while True:
                total = 0
                # fucose
                if draw_table[14 + i * 2 + j] > 0 and draw_table[6 + i * 2 + j] == 0:
                    plot_line(x_bar_bar_bar - (-1) ** i * d, y_bar_bar_bar, 
                              x_bar_bar_bar - (-1) ** i * 2 * d, y_bar_bar_bar, color="k", width=2)
                    obj = plot_triangle(x_bar_bar_bar - (-1) ** i * 2 * d, y_bar_bar_bar, 
                                        orient=select_orient(i), color="r", r=d)
                    ax.add_patch(obj)
                    draw_table[14 + i * 2 + j] -= 1
                    total += draw_table[14 + i * 2 + j]
                # Gal
                if draw_table[10 + i * 2 + j] > 0:
                    plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                              x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_circle(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="y", r=d)
                    ax.add_patch(obj)
                    draw_table[10 + i * 2 + j] -= 1
                    y_bar_bar_bar += 3 * d
                    total += draw_table[10 + i * 2 + j]
                # GlcNAc
                if draw_table[6 + i * 2 + j] > 0:
                    plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                              x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_rectangle(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="b", r=d)
                    ax.add_patch(obj)
                    draw_table[6 + i * 2 + j] -= 1
                    y_bar_bar_bar += 3 * d
                    total += draw_table[6 + i * 2 + j]
                # terminal
                elif draw_table[18 + i * 2 + j] > 0:
                    plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                              x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_diamond(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="m", r=d)
                    ax.add_patch(obj) 
                    draw_table[18 + i * 2 + j] -= 1
                    y_bar_bar_bar += 3 * d
                    total += draw_table[18 + i * 2 + j]
                if total == 0:
                    break 


def glycan_hybrid_plot(start_x, start_y, line_length, ax, table):
    if len(table) != 18:
        return 
    draw_table = table[:]

    x, y = start_x, start_y
    d = line_length

    # fucose core
    if draw_table[2] > 0:
        obj = plot_triangle(x - 2 * d, y, orient="left", color="r", r=d)
        ax.add_patch(obj)
        plot_line(x - d, y, x - 2 * d, y, color="k", width=2)

    # GlcNAc
    while draw_table[0] > 0:
        obj = plot_rectangle(x, y, color="b", r=d)
        ax.add_patch(obj)
        plot_line(x, y + d, x, y + 2 * d, color="k", width=2)
        # update
        y += 3 * d
        draw_table[0] -= 1

    # Man
    if draw_table[1] > 0:
        obj = plot_circle(x, y, color="g", r=d)
        ax.add_patch(obj)      

    # bisect:
    if draw_table[3] > 0:
        plot_line(x, y + d, x, y + 2 * d, color="k", width=2)
        obj = plot_rectangle(x, y + 3 * d, color="b", r=d)
        ax.add_patch(obj)

    # antenna
    gap = compute_gap(d=d/3.0, r=d, angle=61)
    for i in range(2):
        x_bar = x - (-1) ** i * gap / 2.0
        y_bar = compute_antenna_y(y, gap=gap, r=d, d=compute_distance(gap=gap, r=d, angle=61))

        # Man
        if table[4 + i] > 0:
            obj = plot_circle(x_bar, y_bar, color="g", r=d)
            ax.add_patch(obj)
            plot_antenna_line(x, y, gap=gap, orient=select_orient(i), color="k", width=2, r=d, 
                              d=compute_distance(gap=gap, r=d, angle=61))

        gap_bar = compute_gap(d=d, r=d)

        # Man branch
        if i == 0:
            branch_full = (table[6] > 0 and table[7] > 0)
            for j in range(2):
                if draw_table[6 + j] == 0:
                    continue

                x_bar_bar = x_bar
                y_bar_bar = y_bar

                if branch_full:
                    x_bar_bar = x_bar - (-1) ** j * gap_bar / 2.0
                    y_bar_bar = compute_antenna_y(y_bar, gap=gap_bar, r=d, d=compute_distance(gap_bar, d))

                # Man
                if branch_full:
                    obj = plot_circle(x_bar_bar, y_bar_bar, color="g", r=d)
                    ax.add_patch(obj)
                    plot_antenna_line(x_bar, y_bar, gap=gap_bar, 
                                           orient=select_orient(j), color="k", width=2, r=d,
                                           d=compute_distance(gap_bar, d))
                    draw_table[6 + j] -= 1
                else:
                    plot_line(x_bar_bar, y_bar_bar + d, x_bar_bar, y_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_circle(x_bar_bar, y_bar_bar + 3 * d, color="g", r=d)
                    ax.add_patch(obj)
                    draw_table[6 + j] -= 1
                    y_bar_bar += 3 * d

                # chain enlongation
                x_bar_bar_bar = x_bar_bar
                y_bar_bar_bar = y_bar_bar

                # Man
                while draw_table[6 + j] > 0:
                    plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                              x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_circle(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="g", r=d)
                    ax.add_patch(obj)
                    draw_table[6 + j] -= 1
                    y_bar_bar_bar += 3 * d
        # GlcNAc branch
        elif i == 1:
            branch_full = (draw_table[8] > 0 and draw_table[9] > 0)
            for j in range(2):
                if draw_table[8 + j] + draw_table[10 + j] + draw_table[12 + j] + draw_table[14 + j] == 0:
                    continue

                x_bar_bar = x_bar
                y_bar_bar = y_bar

                if branch_full:
                    x_bar_bar = x_bar - (-1) ** j * gap_bar / 2.0
                    y_bar_bar = compute_antenna_y(y_bar, gap=gap_bar, r=d, d=compute_distance(gap_bar, d))

                # GlcNAc
                if branch_full:
                    obj = plot_rectangle(x_bar_bar, y_bar_bar, color="b", r=d)
                    ax.add_patch(obj)
                    plot_antenna_rect_line(x_bar, y_bar, gap=gap_bar, 
                                           orient=select_orient(j), color="k", width=2, r=d,
                                           d=compute_distance(gap_bar, d))
                    draw_table[8 + j] -= 1
                else:
                    plot_line(x_bar_bar, y_bar_bar + d, x_bar_bar, y_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_rectangle(x_bar_bar, y_bar_bar + 3 * d, color="b", r=d)
                    ax.add_patch(obj)
                    draw_table[8 + j] -= 1
                    y_bar_bar += 3 * d


                # chain enlongation
                x_bar_bar_bar = x_bar_bar
                y_bar_bar_bar = y_bar_bar

                while True:
                    total = 0
                    # fucose
                    if draw_table[12 + j] > 0 and draw_table[8 + j] == 0:
                        plot_line(x_bar_bar_bar - (-1) ** i * d, y_bar_bar_bar, 
                                  x_bar_bar_bar - (-1) ** i * 2 * d, y_bar_bar_bar, color="k", width=2)
                        obj = plot_triangle(x_bar_bar_bar - (-1) ** i * 2 * d, y_bar_bar_bar, 
                                            orient=select_orient(i), color="r", r=d)
                        ax.add_patch(obj)
                        draw_table[12 + j] -= 1
                        total += draw_table[12 + j]
                    # Gal
                    if draw_table[10 + j] > 0:
                        plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                                  x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                        obj = plot_circle(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="yellow", r=d)
                        ax.add_patch(obj)
                        draw_table[10 + j] -= 1
                        y_bar_bar_bar += 3 * d
                        total += draw_table[10 + j]
                    # GlcNAc
                    if draw_table[8 + j] > 0:
                        plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                                  x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                        obj = plot_rectangle(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="b", r=d)
                        ax.add_patch(obj)
                        draw_table[8 + j] -= 1
                        y_bar_bar_bar += 3 * d
                        total += draw_table[8 + j]
                    # terminal
                    if draw_table[14 + j] > 0:
                        plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                                  x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                        obj = plot_diamond(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="m", r=d)
                        ax.add_patch(obj)
                        draw_table[14 + j] -= 1
                        y_bar_bar_bar += 3 * d
                        total += draw_table[14 + j]
                    if total == 0:
                        break    

def glycan_highmannose_plot(start_x, start_y, line_length, ax, table):
    if len(table) != 8:
        return
    draw_table = table[:]

    x, y = start_x, start_y
    d = line_length

    # fucose core
    if draw_table[2] > 0:
        obj = plot_triangle(x - 2 * d, y, orient="left", color="r", r=d)
        ax.add_patch(obj)
        plot_line(x - d, y, x - 2 * d, y, color="k", width=2)

    # GlcNAc
    while draw_table[0] > 0:
        obj = plot_rectangle(x, y, color="b", r=d)
        ax.add_patch(obj)
        plot_line(x, y + d, x, y + 2 * d, color="k", width=2)
        # update
        y += 3 * d
        draw_table[0] -= 1

    # Man
    if draw_table[1] > 0:
        obj = plot_circle(x, y, color="g", r=d)
        ax.add_patch(obj)      

    # antenna
    gap = compute_gap(d=d/3.0, r=d, angle=61)
    x_bar = x
    y_bar = y
    for i in range(2):
        x_bar = x - (-1) ** i * gap / 2.0
        y_bar = compute_antenna_y(y, gap=gap, r=d, d=compute_distance(gap=gap, r=d, angle=61))

        # Man
        if table[3 + i] > 0:
            obj = plot_circle(x_bar, y_bar, color="g", r=d)
            ax.add_patch(obj)
            plot_antenna_line(x, y, gap=gap, orient=select_orient(i), color="k", width=2, r=d, 
                              d=compute_distance(gap=gap, r=d, angle=61))

        # Man 1, 2 branch
        gap_bar = compute_gap(d=d, r=d)
        if i == 0 and table[3 + i] > 0:
            x_bar_bar = x_bar
            y_bar_bar = y_bar

            branch_full = (table[5] > 0 and table[6] > 0)
            for j in range(2):
                if branch_full:
                    x_bar_bar = x_bar - (-1) ** j * gap_bar / 2.0
                    y_bar_bar = compute_antenna_y(y_bar, gap=gap_bar, r=d, d=compute_distance(gap_bar, d))

                # Man
                if branch_full:
                    obj = plot_circle(x_bar_bar, y_bar_bar, color="g", r=d)
                    ax.add_patch(obj)
                    plot_antenna_line(x_bar, y_bar, gap=gap_bar, 
                                           orient=select_orient(j), color="k", width=2, r=d,
                                           d=compute_distance(gap_bar, d))
                    draw_table[5 + j] -= 1
                else:
                    plot_line(x_bar_bar, y_bar_bar + d, x_bar_bar, y_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_circle(x_bar_bar, y_bar_bar + 3 * d, color="g", r=d)
                    ax.add_patch(obj)
                    draw_table[5 + j] -= 1
                    y_bar_bar += 3 * d

                # chain enlongation
                x_bar_bar_bar = x_bar_bar
                y_bar_bar_bar = y_bar_bar

                # Man
                while draw_table[5 + j] > 0:
                    plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                              x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                    obj = plot_circle(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="g", r=d)
                    ax.add_patch(obj)
                    draw_table[5 + j] -= 1
                    y_bar_bar_bar += 3 * d   
        elif i == 1 and table[7] > 0:
            x_bar_bar = x_bar
            y_bar_bar = y_bar

            plot_line(x_bar_bar, y_bar_bar + d, x_bar_bar, y_bar_bar + 2 * d, color="k", width=2)
            obj = plot_circle(x_bar_bar, y_bar_bar + 3 * d, color="g", r=d)
            ax.add_patch(obj)
            draw_table[7] -= 1
            y_bar_bar += 3 * d
            # chain enlongation
            x_bar_bar_bar = x_bar_bar
            y_bar_bar_bar = y_bar_bar

            # Man
            while draw_table[7] > 0:
                plot_line(x_bar_bar_bar, y_bar_bar_bar + d, 
                          x_bar_bar_bar, y_bar_bar_bar + 2 * d, color="k", width=2)
                obj = plot_circle(x_bar_bar_bar, y_bar_bar_bar + 3 * d, color="g", r=d)
                ax.add_patch(obj)
                draw_table[7] -= 1
                y_bar_bar_bar += 3 * d   

def glycan_table_plot(start_x, start_y, line_length, ax, table):
    if len(table) == 26:
        glycan_complex_plot(start_x, start_y, line_length, ax, table)
    elif len(table) == 18:
        glycan_hybrid_plot(start_x, start_y, line_length, ax, table)
    elif len(table) == 8:
        glycan_highmannose_plot(start_x, start_y, line_length, ax, table)