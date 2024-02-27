#!/usr/bin/env python

# Import party
print("Import neccessary libraries...")
import pycircos
import matplotlib.pyplot as plt
import collections
print("Import complete.")
print("")

# Set up classes
Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

#Set chromosomes
print("Set chromosomes...")

circle = Gcircle(figsize=(8,8)) 
with open("Rsp4_names.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split(",") 
        name   = line[0]
        length = int(line[-1]) 
        arc    = Garc(arc_id=name, size=length, interspace=1, raxis_range=(850,950), labelposition=120, 
                      label_visible=True, facecolor="#577590") #
        circle.add_garc(arc)

with open("MJ_names_rev.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split(",") 
        name   = line[0]
        length = int(line[-1]) 
        arc2    = Garc(arc_id=name, size=length, interspace=1, raxis_range=(850,950), labelposition=120, 
                       label_visible=True) #facecolor = '#4297A0'
        circle.add_garc(arc2)
        
circle.set_garcs(0,360)

with open("Rsp4_names.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split(",") 
        name   = line[0]
        length = int(line[-1]) 
        circle.tickplot(garc_id=name, tickinterval=5000000, tickwidth='2', raxis_range=(950, 965))
        circle.tickplot(garc_id=name, tickinterval=1000000, raxis_range=(950, 965))

with open("MJ_names_rev.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split(",") 
        name   = line[0]
        length = int(line[-1]) 
        circle.tickplot(garc_id=name, tickinterval=5000000, tickwidth='2', raxis_range=(950, 965))
        circle.tickplot(garc_id=name, tickinterval=1000000, raxis_range=(950, 965))

print("Set chromosomes complete.")
print("")

# Repeats bar plot
print("Set repeat plots...")

values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("Rsp4_repeats.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))

vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin, vmax],
                   raxis_range=[740,840], facecolor="#F94144", spine=True)

values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("MJ_repeats.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))

vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin, vmax],
                   raxis_range=[740,840], facecolor="#F94144", spine=True)

print("Set repeat plots complete.")
print("")

#Genes bar plot
print("Set gene plots...")

values_genes   = [] 
genedata_dict = collections.defaultdict(dict)
with open("Rsp4_genes.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in genedata_dict:
            genedata_dict[name]["positions"] = []
            genedata_dict[name]["widths"]    = [] 
            genedata_dict[name]["values"]    = [] 
        genedata_dict[name]["positions"].append(start) 
        genedata_dict[name]["widths"].append(width)
        genedata_dict[name]["values"].append(float(line[-1]))
        values_genes.append(float(line[-1]))

vmin, vmax = min(values_genes), max(values_genes) 
for key in genedata_dict:  
    circle.barplot(key, data=genedata_dict[key]["values"], positions=genedata_dict[key]["positions"], 
                   width=genedata_dict[key]["widths"], base_value=0.0, rlim=[vmin, vmax],
                   raxis_range=[630,730], facecolor="#90BE6D", spine=True)

values_genes   = [] 
genedata_dict = collections.defaultdict(dict)
with open("MJ_genes.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in genedata_dict:
            genedata_dict[name]["positions"] = []
            genedata_dict[name]["widths"]    = [] 
            genedata_dict[name]["values"]    = [] 
        genedata_dict[name]["positions"].append(start) 
        genedata_dict[name]["widths"].append(width)
        genedata_dict[name]["values"].append(float(line[-1]))
        values_genes.append(float(line[-1]))

vmin, vmax = min(values_genes), max(values_genes) 
for key in genedata_dict:  
    circle.barplot(key, data=genedata_dict[key]["values"], positions=genedata_dict[key]["positions"], 
                   width=genedata_dict[key]["widths"], base_value=0.0, rlim=[vmin, vmax],
                   raxis_range=[630,730], facecolor="#90BE6D", spine=True)

print("Set gene plots complete.")
print("")

# Linkplot
print("Set link plot...")

values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open("Rsp4_vs_MJ_link.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name1  = line[0]     
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = line[3]     
        start2 = int(line[4])-1
        end2   = int(line[5])
        source = (name1, start1, end1, 620)
        destination = (name2, start2, end2, 620)
        circle.chord_plot(source, destination, facecolor=circle.garc_dict[name2].facecolor)

print("Set link plot complete.")
print("")

# Save output
print("Exporting figures...")

circle.save(file_name="circos", format="jpg")
#circle.save(file_name="circos", format="pdf")

print("Circos plot complete.")