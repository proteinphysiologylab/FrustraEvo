#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder
import numpy as np
import os

#inputs

job_id=sys.argv[1]
structure_pdb='./../FrustraEvo_'+job_id+'/Data/'+sys.argv[2]+'.done/VisualizationScrips/'+sys.argv[2]+'.pdb'
pdb_id=sys.argv[2]
IC_SR='./../FrustraEvo_'+job_id+'/OutPutFiles/'+sys.argv[3]
IC_MUT='./../FrustraEvo_'+job_id+'/OutPutFiles/'+sys.argv[4]
IC_CONF='./../FrustraEvo_'+job_id+'/OutPutFiles/'+sys.argv[5]

# Fetch the contents of the file using urllib.request.urlopen()
def create_coloring(infoUrl):
    
    with open(infoUrl, 'r') as f:
        infoFile = f.read()

    residues = [
        {'number':line.split('\t')[2],'color':{'r':0,'g':0,'b':0} if float(line.split('\t')[13]) < 0.5 
         else ({'r':0,'g':255,'b':0} if line.split('\t')[14] == 'MIN' 
               else ({'r':255,'g':0,'b':0} if line.split('\t')[14] == 'MAX' 
                     else ({'r':128,'g':128,'b':128} if line.split('\t')[14] == 'NEU' else None
                    )
                )
            )
        }
        for line in infoFile.strip().split('\n')[1:]
    ]
    return residues

colorString=str(create_coloring(IC_SR))

parser = PDBParser()
structure = parser.get_structure(pdb_id, structure_pdb)

#this retrieves the coordinates of a residue
def get_residue_coordinates(structure, chain_id, residue_id):
   for model in structure:
       for chain in model:
           if chain.id == chain_id:
               for residue in chain:
                  if residue.id[1] == residue_id:
                      return residue['CA'].get_coord()
   return None

def format_row(row):
   print(row)
   a = ', '.join(map(str, row['coord_Res1']))
   b = ', '.join(map(str, row['coord_Res2']))
   color = row['contact_color']
   return f"\t\t\t{{a: Vec3.create({a}), b: Vec3.create({b}), color: Color({color})}}"

#process Configurational data
IC_conf=pd.read_csv(IC_CONF, sep='\t')
IC_conf=IC_conf[(IC_conf['FreqConts'].astype(float)>0.5)&(IC_conf['ICtotal'].astype(float)>0.5)][['NumRes1_Ref','NumRes2_Ref','Chain1_Ref','Chain2_Ref','FstConserved']]
# Apply the function to each row and create the new columns
IC_conf['coord_Res1'] = IC_conf.apply(lambda row: get_residue_coordinates(structure, row['Chain1_Ref'] ,row['NumRes1_Ref']), axis=1)
IC_conf['coord_Res2'] = IC_conf.apply(lambda row: get_residue_coordinates(structure,row['Chain2_Ref'] ,row['NumRes2_Ref']), axis=1)
IC_conf['contact_color'] = np.where(IC_conf['FstConserved'] == 'MAX', '0xff0000',
                            np.where(IC_conf['FstConserved'] == 'NEU', '0x808080',
                                    np.where(IC_conf['FstConserved'] == 'MIN', '0x00ff00', 'Unknown')))

#create the coloring vector
IC_conf['formatted_string'] = IC_conf.apply(format_row, axis=1)

formatted_string_concat_conf = '\n'.join(IC_conf['formatted_string'].tolist())

formatted_strings_conf = []

for index, row in IC_conf.iterrows():
    formatted_strings_conf.append(row['formatted_string'] + ',')

formatted_string_concat_conf = '\n'.join(formatted_strings_conf)

#process Mutational data
IC_mut=pd.read_csv(IC_MUT, sep='\t')
IC_mut=IC_mut[(IC_mut['FreqConts'].astype(float)>0.5)&(IC_mut['ICtotal'].astype(float)>0.5)][['NumRes1_Ref','NumRes2_Ref','Chain1_Ref','Chain2_Ref','FstConserved']]
# Apply the function to each row and create the new columns
IC_mut['coord_Res1'] = IC_mut.apply(lambda row: get_residue_coordinates(structure, row['Chain1_Ref'], row['NumRes1_Ref']), axis=1)
IC_mut['coord_Res2'] = IC_mut.apply(lambda row: get_residue_coordinates(structure, row['Chain2_Ref'], row['NumRes2_Ref']), axis=1)
IC_mut['contact_color'] = np.where(IC_mut['FstConserved'] == 'MAX', '0xff0000',
                            np.where(IC_mut['FstConserved'] == 'NEU', '0x808080',
                                    np.where(IC_mut['FstConserved'] == 'MIN', '0x00ff00', 'Unknown')))

#create the coloring vector
IC_mut['formatted_string'] = IC_mut.apply(format_row, axis=1)

formatted_string_concat_mut = '\n'.join(IC_mut['formatted_string'].tolist())

formatted_strings_mut = []

for index, row in IC_mut.iterrows():
    formatted_strings_mut.append(row['formatted_string'] + ',')

formatted_string_concat_mut = '\n'.join(formatted_strings_mut)

htmlString=('''
    
<!DOCTYPE html>
<html lang="en">
    <head>
    	<link rel="stylesheet" type="text/css" href="https://www.ebi.ac.uk/pdbe/pdb-component-library/css/pdbe-molstar-3.1.0.css">
    	<script type="text/javascript" src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdbe-molstar-plugin-3.1.0.js"></script>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">
        <title>Mol* Lighting Demo</title>
        <style>
            * {
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }
           #visualizers {
  	        	display: flex; /* Add this line to apply Flexbox layout */
    	        height: 100vh;
  		        gap: 40px;
	         }
            #app {
                position: relative;
                width: 400px;
                height: 400px;
            }
            #app2 {
                position: relative;
                width: 400px;
                height: 400px;
            }
            #app3 {
                position: relative;
                width: 400px;
                height: 400px;
            }
        </style>
        <link rel="stylesheet" type="text/css" href="molstar.css" />
        <script type="text/javascript" src="./molstar.js"></script>
    </head>
    <body>
        <div id='controls'></div>
        <div id="visualizers">
            <div id="app"></div>
            <div id="app2"></div>
            <div id="app3"></div>
        </div>
        <script>
 
        //Create plugin instance
    		var viewerInstance = new PDBeMolstarPlugin();

    	//Set options (Checkout available options list in the documentation)
    		var options = {
    	 	 customData: {
    	 	   url: './'''+pdb_id+'''.pdb',
    	  	  format: 'pdb'
    		  },
    		  alphafoldView: true,
    	 	 bgColor: { r : 255, g : 255, b : 255 },
   		 };


    	//Get element from HTML/Template to place the viewer 
  	 	 var viewerContainer = document.getElementById('app');

   		 //Call render method to display the 3D view
   		 viewerInstance.render(viewerContainer, options);


   		 function colorResidues(residueData) {
    		  viewerInstance.visual.select({
    		 	   data: residueData.map(function(residue) {
    		      return { residue_number: residue.number, color: residue.color };
     		   })
     		 });
    		};
    		
    		
    
   		 window.onload = function() {
   	 	  setTimeout(function() {
    		    colorResidues('''+colorString.replace("'","")+''');
   		 }, 500);
   		 }

	    	LightingDemo.init('app2').then(() => {
			LightingDemo.load(
		    	{ url: './'''+pdb_id+'''.pdb', format: 'pdb', isBinary: false, assemblyId: 0 }, 
		    	5, 
		    	1.3, 
		    	['''+formatted_string_concat_mut+'''
		    	]		    	
			);
	    	});
	    	LightingDemo2.init('app3').then(() => {
			LightingDemo2.load(
		    	{ url: './'''+pdb_id+'''.pdb', format: 'pdb', isBinary: false, assemblyId: 0 }, 
		    	5, 
		    	1.3, 
		    	['''+formatted_string_concat_conf+'''
		    	]		    	
			);
	    	});
        </script>
    </body>
</html>
''')    
    
os.system('mkdir ./../FrustraEvo_'+job_id+'/molstar')


with open('./../FrustraEvo_'+job_id+'/molstar/molstar_mut_aux.txt','w') as file:
    file.write(formatted_string_concat_mut)
    
with open('./../FrustraEvo_'+job_id+'/molstar/molstar_mut_aux.txt', 'r') as file:
    lines = [line.strip('\t\n') for line in file]

# Concatenas todas las líneas en una sola cadena
concatenated_string = ''.join(lines)

# Luego, puedes escribir la cadena concatenada en otro archivo o utilizarla como desees
with open('./../FrustraEvo_'+job_id+'/molstar/molstar_mut.txt', 'w') as output_file:
    output_file.write(concatenated_string)

with open('./../FrustraEvo_'+job_id+'/molstar/molstar_conf_aux.txt','w') as file:
    file.write(formatted_string_concat_conf)
    
with open('./../FrustraEvo_'+job_id+'/molstar/molstar_conf_aux.txt', 'r') as file:
    lines = [line.strip('\t\n') for line in file]

# Concatenas todas las líneas en una sola cadena
concatenated_string = ''.join(lines)

# Luego, puedes escribir la cadena concatenada en otro archivo o utilizarla como desees
with open('./../FrustraEvo_'+job_id+'/molstar/molstar_conf.txt', 'w') as output_file:
    output_file.write(concatenated_string)

os.system('rm ./../FrustraEvo_'+job_id+'/molstar/molstar_conf_aux.txt')
os.system('rm ./../FrustraEvo_'+job_id+'/molstar/molstar_mut_aux.txt')
