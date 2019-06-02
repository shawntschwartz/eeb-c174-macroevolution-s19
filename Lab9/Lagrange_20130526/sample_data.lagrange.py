#!/usr/bin/env python2.7
import os
import lagrange
data = """\
### begin data
{'area_adjacency': [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 0], [1, 1, 0, 1]],
 'area_dispersal': [[[1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0, 1.0]]],
 'area_labels': ['A', 'B', 'C', 'D'],
 'base_rates': '__estimate__',
 'dispersal_durations': [11.0],
 'dm_symmetric_entry': True,
 'excluded_ranges': [],
 'lagrange_version': '20130526',
 'max_range_size': 3,
 'model_name': 'sample_data',
 'newick_trees': [{'included': '__all__',
                   'name': 'Tree0',
                   'newick': '(((t1:6.681697201, (((t6:1.943860938, t7:1.943860938)N3:0.2564040428, (t11:1.619556123, ((t17:0.8605859135, t18:0.8605859135)N7:0.07537970442, t16:0.9359656179)N9:0.6835905054)N10:0.5807088571)N11:2.036408271, t2:4.236673252)N13:2.445023949)N14:1.756342702, (((t10:1.730346866, ((t14:1.019987943, t15:1.019987943)N18:0.009134121693, t13:1.029122065)N20:0.7012248017)N21:0.9332721617, t5:2.663619028)N23:1.271652841, (((t19:0.6301399509, t20:0.6301399509)N26:0.8521645762, t12:1.482304527)N28:1.401751983, t4:2.88405651)N30:1.051215359)N31:4.502768034)N32:1.561960098, (t3:3.240695461, (t8:1.812839413, t9:1.812839413)N36:1.427856048)N37:6.759304539)N38;',
                   'root_age': 10.0000000014}],
 'ranges': [(),
            (0,),
            (0, 1),
            (0, 1, 2),
            (0, 1, 3),
            (0, 2),
            (0, 2, 3),
            (0, 3),
            (1,),
            (1, 2),
            (1, 2, 3),
            (1, 3),
            (2,),
            (3,)],
 'taxa': ['t1',
          't2',
          't3',
          't4',
          't5',
          't6',
          't7',
          't8',
          't9',
          't10',
          't11',
          't12',
          't13',
          't14',
          't15',
          't16',
          't17',
          't18',
          't19',
          't20'],
 'taxon_range_data': {'t1': (0,),
                      't10': (0,),
                      't11': (1,),
                      't12': (0,),
                      't13': (0,),
                      't14': (0, 1),
                      't15': (0,),
                      't16': (1,),
                      't17': (1,),
                      't18': (3,),
                      't19': (3,),
                      't2': (0, 1),
                      't20': (3,),
                      't3': (0,),
                      't4': (0,),
                      't5': (0,),
                      't6': (1, 2),
                      't7': (1,),
                      't8': (1,),
                      't9': (0, 1)}}
### end data
"""

i = 0
while 1:
    if not i:
        outfname = "sample_data.results.txt"
    else:
        outfname = "sample_data.results-"+str(i)+".txt"
    if not os.path.exists(outfname): break
    i += 1
outfile = open(outfname, "w")
lagrange.output.log(lagrange.msg, outfile, tee=True)
model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(data)
lagrange.output.ascii_tree(outfile, tree, model, data, tee=True)
if base_rates != "__estimate__":
    d, e = base_rates
else:
    d, e = lagrange.output.optimize_dispersal_extinction(outfile, tree, model, tee=True)
if nodelabels:
    if nodelabels == "__all__":
        nodelabels = None
    lagrange.output.ancsplits(outfile, tree, model, d, e, nodelabels=nodelabels, tee=True)
