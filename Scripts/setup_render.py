import sys

# Define the URL of the 'info.txt' file
infoUrl = "./../FrustraEvo_"+sys.argv[1]+"/OutPutFiles/IC_"+sys.argv[2]+"_"+sys.argv[1]

# Fetch the contents of the file using urllib.request.urlopen()
def create_coloring(infoUrl):
    
    with open(infoUrl, 'r') as f:
        infoFile = f.read()
    num='number'
    color='color'
    r='r'
    g='g'
    b='b'
    residues = [
        {num.strip('\''):line.split('\t')[2].strip('\''),color.strip('\''):{r.strip('\''):0,g.strip('\''):0,b.strip('\''):0} if float(line.split('\t')[13]) < 0.5 
         else ({r.strip('\''):0,g.strip('\''):255,b.strip('\''):0} if line.split('\t')[14] == 'MIN' 
               else ({r.strip('\''):255,g.strip('\''):0,b.strip('\''):0} if line.split('\t')[14] == 'MAX' 
                     else ({r.strip('\''):128,g.strip('\''):128,b.strip('\''):128} if line.split('\t')[14] == 'NEU' else None
                    )
                )
            )
        }
        for line in infoFile.strip().split('\n')[1:]
    ]
    return residues
    
with open (str(infoUrl)+'.txt', 'w') as o:
    o.write(str(create_coloring(infoUrl)))
