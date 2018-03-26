counter = 0
for p in phenotypes:
    pg = [x[1] for x in ph_gn_links if x[0]==p]
    counter += len(set([x[0] for x in go_gn_links if x[1] in pg]))
print counter
 

 
