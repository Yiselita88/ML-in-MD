import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import Bio.PDB as bp
from Bio.PDB import *
import os, re
from colabdesign import mk_afdesign_model, clear_mem
from colabdesign.af.alphafold.common import residue_constants
from colabdesign.shared.protein import _np_get_cb
import pickle
from colabdesign import af
#from IPython.display import HTML
#from google.colab import files
import numpy as np
from scipy.special import softmax
import sys
import tqdm.notebook
import argparse

#TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

def parse_args():                              #in line argument parser with help 
    parser = argparse.ArgumentParser()
    parser.add_argument('-seq',type=str,default='',help='input sequence file')
    parser.add_argument('-recycle',type=int,help='number of recycles')
    parser.add_argument('-binder_len',type=int,help='length of binder')
    parser.add_argument('-pdb',default=None,type=str,help='pdb path')
    parser.add_argument('-target_chain',type=str,help='index of target chain, e.g. A')
    parser.add_argument('-target_hotspot',type=str,help='residue index of selected hot spot to design binder around, e.g. 1-7,12,15')
    parser.add_argument('-flexible',action="store_true",help='predefine if the binder is flexible or not')
    parser.add_argument('-use_multimer',action="store_true",help='define if using multimer or not')
    parser.add_argument('-num_models',type=str,default=5,help='select number of AF models')
    parser.add_argument('-protocol',type=str,default="pssm_semigreedy",help='peptide design protocol')
    parser.add_argument('-optimizer',type=str,default="sgd",help='peptide design optimizer')
    parser.add_argument('-lr',type=float,default=0.1,help='optimization learning rate')
    parser.add_argument('-dropout',action="store_true",help='dropout or not')
    return parser.parse_args()
 
args=parse_args()

def add_cyclic_offset(self, offset_type=2):
  '''add cyclic offset to connect N and C term'''
  def cyclic_offset(L):
    i = np.arange(L)
    ij = np.stack([i,i+L],-1)
    offset = i[:,None] - i[None,:]
    c_offset = np.abs(ij[:,None,:,None] - ij[None,:,None,:]).min((2,3))
    if offset_type == 1:
      c_offset = c_offset
    elif offset_type >= 2:
      a = c_offset < np.abs(offset)
      c_offset[a] = -c_offset[a]
    if offset_type == 3:
      idx = np.abs(c_offset) > 2
      c_offset[idx] = (32 * c_offset[idx] )/  abs(c_offset[idx])
    return c_offset * np.sign(offset)
  idx = self._inputs["residue_index"]
  offset = np.array(idx[:,None] - idx[None,:])

  if self.protocol == "binder":
    c_offset = cyclic_offset(self._binder_len)
    offset[self._target_len:,self._target_len:] = c_offset

  if self.protocol in ["fixbb","partial","hallucination"]:
    Ln = 0
    for L in self._lengths:
      offset[Ln:Ln+L,Ln:Ln+L] = cyclic_offset(L)
      Ln += L
  self._inputs["offset"] = offset

def get_pdb(pdb_code=""):
  if pdb_code is None or pdb_code == "":
    upload_dict = files.upload()
    pdb_string = upload_dict[list(upload_dict.keys())[0]]
    with open("tmp.pdb","wb") as out: out.write(pdb_string)
    return "tmp.pdb"
  elif os.path.isfile(pdb_code):
    return pdb_code
  elif len(pdb_code) == 4:
    os.system(f"wget -qnc https://files.rcsb.org/view/{pdb_code}.pdb")
    return f"{pdb_code}.pdb"
  else:
    os.system(f"wget -qnc https://alphafold.ebi.ac.uk/files/AF-{pdb_code}-F1-model_v3.pdb")
    return f"AF-{pdb_code}-F1-model_v3.pdb"


#@title initialize the model
binder_len = args.binder_len #@param {type:"integer"}
#@markdown Provide a starting point (optional)
binder_seq = open(args.seq).readlines()[-1]
binder_seq = re.sub("[^A-Z]", "", binder_seq.upper())
print(f"run seq {binder_seq} with recycle {args.recycle}")
target_chain = args.target_chain
target_hotspot = args.target_hotspot
target_flexible = args.flexible
pdb = args.pdb
if target_hotspot == "": target_hotspot = None

from colabdesign.af.alphafold.common import protein, residue_constants
aa_order = residue_constants.restype_order
#@markdown Experimental options
num_models = "" #@param ["1", "2", "3", "4", "5", "all"]
num_models = 5 if num_models == "all" else int(num_models or 5) # from chatGPT
#num_models = 5 if num_models == "all" else int(num_models) # temporary commented in case I need to put it back
use_multimer = args.use_multimer 
num_recycles = args.recycle
#template_path = args.template_path
##@markdown - `xyz` - use structure output as template input
#@markdown - `dgram` - use distogram output as template input
#@markdown - `dgram_retrain` - replace distogram head from AlphaFold with one retrained to map output bins to template bins.

if len(binder_seq) > 0:
  binder_len = len(binder_seq)

binder_chain = "" #@param {type:"string"}
if binder_chain == "": binder_chain = None

x = {"pdb_filename":pdb,
     "chain":target_chain,
     "binder_len":binder_len,
     "binder_chain":binder_chain,
     "hotspot":target_hotspot,
     "use_multimer":use_multimer,
     "rm_target_seq":target_flexible}
     

if "x_prev" not in dir() or x != x_prev:
  clear_mem()
  model = mk_afdesign_model(protocol="binder",
                            use_multimer=x["use_multimer"],
                            num_recycles=num_recycles,
                            recycle_mode="sample")
  model.prep_inputs(**x,
                    ignore_missing=False)
  add_cyclic_offset(af_model, offset_type=2)
  x_prev = copy_dict(x)
  print("target length:", model._target_len)
  print("binder length:", model._binder_len)
  binder_len = model._binder_len

#@title **run AfDesign**
from scipy.special import softmax
optimizer = args.protocol
#@markdown #### advanced GD settings
GD_method = args.optimizer 
learning_rate = args.lr 
norm_seq_grad = True 
dropout = True 

model.restart(seq=binder_seq)
model.set_optimizer(optimizer=GD_method,
                    learning_rate=learning_rate,
                    norm_seq_grad=norm_seq_grad)
models = model._model_names[:num_models]

flags = {"num_recycles":num_recycles,
         "models":models,
         "dropout":dropout}

if optimizer == "3stage":
  model.design_3stage(120, 60, 10, **flags)
  pssm = softmax(model._tmp["seq_logits"],-1)

if optimizer == "pssm_semigreedy":
  model.design_pssm_semigreedy(120, 32, **flags)
  pssm = softmax(model._tmp["seq_logits"],1)

if optimizer == "semigreedy":
  model.design_pssm_semigreedy(0, 32, **flags)
  pssm = None

if optimizer == "pssm":
  model.design_logits(120, e_soft=1.0, num_models=1, ramp_recycles=True, **flags)
  model.design_soft(32, num_models=1, **flags)
  flags.update({"dropout":False,"save_best":True})
  model.design_soft(10, num_models=num_models, **flags)
  pssm = softmax(model.aux["seq"]["logits"],-1)

O = {"logits":model.design_logits,
     "soft":model.design_soft,
     "hard":model.design_hard}

if optimizer in O:
  O[optimizer](120, num_models=1, ramp_recycles=True, **flags)
  flags.update({"dropout":False,"save_best":True})
  O[optimizer](10, num_models=num_models, **flags)
  pssm = softmax(model.aux["seq"]["logits"],-1)

model.save_pdb(f"{model.protocol}.pdb")
