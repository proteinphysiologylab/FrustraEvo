import os
import sys
import time

# Registra el tiempo de inicio
start_time = time.time()

sys.path.append('')#Path to Functions.py file
import Functions
import argparse

parser = argparse.ArgumentParser(description='Calculation of the frustration logo.')
parser.add_argument("--JobId", help="The name of the job")
parser.add_argument("--RPath", default='Scripts', help="Path to R script files (Default: Scripts)")
parser.add_argument("--fasta", help="Name of the fasta")
parser.add_argument("--ref", default='None', help="Id of the reference protein of your logo (Default: None)")
parser.add_argument("--pdb_db", default='None', help="Path to the PDBs folder (Default: None)")
parser.add_argument("--cmaps", default='None', help="Put 'yes' for contactmaps calculation (Default: None)")


#How to run the pipeline in linux terminal:
#python run_logo.py JobId RPath fasta list_file ref pdb_db 

args = parser.parse_args()
list_file=''

#out_log=open('FrustraEvo_'+args.JobId+'/CheckPoints_FE','w')
#the parameter Scripts is the path to the .r and .py files for plots
seq=0
#out_log.write('Copyng files for Frustration Logo\n')
Functions.copyfiles(args.JobId,args.RPath,args.pdb_db)
list_file=Functions.pdb_list(args.fasta)
#out_log.write('Changes in MSA Frustration Logo\n')
Functions.changes(args.JobId,args.fasta)
#out_log.write('Running Checks in sequence\n')
Functions.checks_seq(list_file,args.JobId,args.pdb_db)
path_file='FrustraEvo_'+args.JobId+'/ErrorSeq.log'
out_log=open('results/CheckPoints_FE','a')
out_log.write(' Reading and Preparing Files...\n')
out_log.close()
if not os.path.exists(path_file):
  out_log=open('results/CheckPoints_FE','a')
  out_log.write(' Calculating Single Residue Frustration Index...\n')
  out_log.close()
  Functions.FrustraPDB(list_file,args.JobId,args.pdb_db)
  out_log=open('results/CheckPoints_FE','a')
  out_log.write(' Running Internal Checks...\n')
  out_log.close()
  Functions.checks(args.JobId)
  out_log=open('results/CheckPoints_FE','a')
  out_log.write(' Preparing MSA Files...\n')
  out_log.close()
  Functions.prepare_file(args.JobId,args.ref)
  Functions.FinalAlign(args.JobId)
#  out_log=open('results/CheckPoints_FE','a')
#  out_log.write('Running Equivalences...\n')
#  out_log.close()
  Functions.Equivalences(args.JobId)
  out_log=open('results/CheckPoints_FE','a')
  out_log.write(' Calculating SeqIC and FrustIC...\n')
  out_log.close()
  Functions.FastaMod(args.JobId)
#  out_log=open('results/CheckPoints_FE','a')
#  out_log.write('Running Checks\n')
#  out_log.close()
  Functions.LogoCheck(args.JobId)
  out_log=open('results/CheckPoints_FE','a')
  out_log.write(' Generating Sequence and Frustration Logos Plots...\n')
  out_log.close()
  Functions.plots_logo(args.JobId,args.ref,args.RPath) 
  if args.cmaps == 'yes':
        out_log=open('results/CheckPoints_FE','a')
        out_log.write(' Calculating Mutational Frustration Index and Contact Maps...\n')
        out_log.close()
        Functions.CMaps_Mutational(args.JobId,args.RPath,args.ref)#Genera los mapas de contacto para el indice mutational
        out_log=open('results/CheckPoints_FE','a')
        out_log.write(' Calculating Configurational Frustration Index and Contact Maps...\n')
        out_log.close()
        Functions.CMaps_Configurational(args.JobId,args.RPath,args.ref)#Genera los mapas de contacto para el indice configurational
  Functions.clean_files(args.JobId,args.RPath,args.ref, cmaps) #add cmaps  
  end_time = time.time()
  elapsed_time = end_time - start_time
#path_direc='FrustraEvo_'+args.JobId
#out=open(path_direc+'/'+JodID+'.log','w')
#out.write(str(elapsed_time)+'\n')
#out.close()
  out_log=open('results/CheckPoints_FE','a')
  out_log.write(' Generating Output Files...\n')
  out_log.close()
  Functions.VScript(args.JobId,elapsed_time)
#  os.system('cd '+args.RPath+';python3 setup_render.py '+args.JobId+' SingleRes')
#  os.system('cd '+args.RPath+';python3 Seq_IC.py '+args.JobId)
  os.system(f'python3 {args.RPath}/setup_render.py {args.JobId} SinfleRes') #avoid cd
  os.system(f'python3 {args.RPath}/Seq_IC.py {args.JobId}') #avoid cd
#contact_maps.py	
  if args.cmaps == 'yes': #add clause to check for cmaps
  #os.system('cd '+args.RPath+';python3 contact_maps.py '+args.JobId+' '+args.ref+' IC_SingleRes_'+args.JobId+' IC_Mut_'+args.JobId+' IC_Conf_'+args.JobId)
        os.system(f'python3 {args.RPath}/contact_maps.py {args.JobId} {args.ref} IC_SingleRes_{args.JobId} IC_Mut_{args.JobId} IC_Conf_{args.JobId}')#avoid changing dir
out_log=open('results/CheckPoints_FE','a')
out_log.write(' Job finished')
out_log.close()

