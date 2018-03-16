#!/usr/bin/env python

import os, sys, glob, csv
from numpy import *
from scipy.stats import binom_test
from scipy.stats import binom
import argparse
import logging

__author__ = 'James Sun'
__contact__ = 'jsun@foundationmedicine.com'
__version__ = '1.0.0'
__doc__ = '''FMI SGZ method to evaluate CNA origin and zygosity.'''


# Functions for the CNA model
def cn2lr_bl(p,bl,cn):
   lr = log2((p*cn + 2*(1-p)) / bl)
   return lr


def core_SGZ(data_CNA, short_variants):

   ALPHA = 0.01
   
   y      = []
   
   for sv in short_variants:
      if sv.chr_.upper() in ['X', 'Y']:
         continue

      if sv.frequency>1:
         sv.frequency = 1.0
   
      p = M = C = nan
      
      segFound = 0

      clonality = 'NA'
      in_tumor  = 'NA'
      
      for seg in data_CNA:
         
         if seg.chr_== sv.chr_:
            if seg.segStart<=sv.position<=seg.segEnd:
               
               segFound = 1
               
               p  = seg.purity
               M  = seg.numMAtumorPred
               C  = seg.CN
               bL = seg.baseLevel
               
               model_MAF = seg.mafPred
               data_MAF  = seg.segMAF
               error_MAF = abs(model_MAF-data_MAF)
               
               model_LR  = cn2lr_bl(p,bL,C)
               data_LR   = seg.segLR
               error_LR  = abs(model_LR-data_LR)
               
               if M==0:
                  M=C-M
               
               P_G1 = P_G2 = P_S1 = P_S2 = 0
               
               AF_G1 = (p*M + 1*(1-p))/(p*C+2*(1-p))
               AF_S1 = (p*M + 0*(1-p))/(p*C+2*(1-p))
               
               P_G1  = binom_test(round(sv.depth*sv.frequency), sv.depth, AF_G1)
               P_S1  = binom_test(round(sv.depth*sv.frequency), sv.depth, AF_S1)
               
               if M!=C-M:
                  AF_G2 = (p*(C-M) + 1*(1-p))/(p*C+2*(1-p))
                  P_G2  = binom_test(round(sv.depth*sv.frequency), sv.depth, AF_G2)
                  
                  if C-M!=0:
                     AF_S2 = (p*(C-M) + 0*(1-p))/(p*C+2*(1-p))
                     P_S2  = binom_test(round(sv.depth*sv.frequency), sv.depth, AF_S2)
                  else:
                     AF_S2 = P_S2 = nan
                     
               else:
                  AF_G2 = AF_S2 = P_G2 = P_S2 = nan
               
               
               
               max_prob_germline = max(P_G1, P_G2)
               max_prob_somatic  = max(P_S1, P_S2)
               
               if   max_prob_somatic  == 0:  logodds =  inf
               elif max_prob_germline == 0:  logodds = -inf
               else:                         logodds = log10(max_prob_germline) - log10(max_prob_somatic)
               
              
               if 'CYP2D6' in sv.mutation:
                  SG_status = 'nocall_CYP2D6'
                  
               elif 'HLA' in sv.mutation:
                  SG_status = 'nocall_HLA'
                  
               elif sv.frequency<0.05 and 0.2<p<0.9:
                  SG_status = 'subclonal somatic'
                  clonality = 'subclonal'
                  in_tumor  = 'in tumor'
               
               elif isnan(C) or isnan(M):
                  SG_status = 'ambiguous_CNA_model'
                  
               elif p>0.95:
                  SG_status = 'nocall_purity>95%'
                  
               elif sv.frequency>0.95:
                  SG_status = 'germline'
                  clonality = 'clonal'
                  in_tumor  = 'in tumor'
                  
               elif error_MAF>0.06 or error_LR>0.4:
                  SG_status = 'ambiguous_CNA_model'
               
               elif max_prob_germline>ALPHA and max_prob_somatic<ALPHA:
                  
                  if logodds<2:
                     SG_status = 'probable germline'
                  else:
                     SG_status = 'germline'
                  
                  clonality = 'clonal'
                  
                  if nanargmax([P_G1, P_G2])==0:
                     in_tumor = 'in tumor'
                  else:
                     M = C-M
                     if M==0:
                        in_tumor = 'not in tumor'
                  
                  
               elif max_prob_germline<ALPHA and max_prob_somatic>ALPHA:

                  if logodds>-2:
                     SG_status = 'probable somatic'
                  else:
                     SG_status = 'somatic'
                     
                  clonality = 'clonal'
                  in_tumor  = 'in tumor'
                  
                  if nanargmax([P_S1, P_S2])==1:
                     M = C-M
                  
                  
               elif max_prob_germline>ALPHA and max_prob_somatic>ALPHA:
                  SG_status = 'ambiguous_both_G_and_S'
   
               elif max_prob_germline<ALPHA and max_prob_somatic<ALPHA:
                  
                  min_soma_EAF = min(AF_S1, AF_S2)
                  min_germ_EAF = min(AF_G1, AF_G2)
                  
                  if p>=0.3 and sv.frequency<0.25 and sv.frequency<min_soma_EAF/1.5 and min_soma_EAF<=min_germ_EAF:
                     SG_status = 'subclonal somatic'
                     clonality = 'subclonal'
                     in_tumor  = 'in tumor'
                     
                  elif p>=0.3 and sv.frequency<0.25 and sv.frequency<min_germ_EAF/2.0 and min_germ_EAF<min_soma_EAF:
                     SG_status = 'subclonal somatic'
                     clonality = 'subclonal'
                     in_tumor  = 'in tumor'

                  elif logodds<-5 and max_prob_somatic>1e-10:
                     SG_status = 'somatic'
                     clonality = 'clonal'
                     in_tumor  = 'in tumor'
                     
                     if nanargmax([P_S1, P_S2])==1:
                        M = C-M

                  elif logodds>5 and max_prob_germline>1e-4:
                     SG_status = 'germline'
                     clonality = 'clonal'
                     
                     if nanargmax([P_G1, P_G2])==0:
                        in_tumor = 'in tumor'
                     else:
                        M = C-M
                        if M==0:
                           in_tumor = 'not in tumor'
                     
                  else:
                     SG_status = 'ambiguous_neither_G_nor_S'
                  
               else:
                  SG_status = 'unknown'
               
               break
   

      if segFound==0:
         AF_G1 = AF_G2 = AF_S1 = AF_S2 = p = C = M = error_LR = error_MAF = max_prob_germline = max_prob_somatic = logodds = nan
         SG_status = 'nocall_segmentMissing'
      
      
      zygosity = 'NA'
      if isnan(C) or isnan(M) or p<0.19:
         zygosity = 'NA'
         if 'germline' in SG_status:
            in_tumor = 'NA'
            
      elif C==0 and M==0:   zygosity = 'homoDel'
      elif C==M and M==1:   zygosity = 'homozygous'
      elif C==M and M>=2:   zygosity = 'homozygous'
      elif C>=1 and M==0:   zygosity = 'not in tumor'
      elif C!=M and M!=0:   zygosity = 'het'
      
      # calculate allele burden (quantitative clonality)
      allele_burden = nan
      burden_CI     = [nan,nan]
      if ('germline' in SG_status) or ('somatic' in SG_status):
         
         if 'germline' in SG_status:
            EAF = (p*M+1-p)/(p*C+2*(1-p))
         elif 'somatic' in SG_status:
            EAF = p*M/(p*C+2*(1-p))
         
         allele_burden = sv.frequency/EAF
         
         try:
            # sim 1000 times and get 95% CI of the allele burden
            # EAF_sim   = binom.rvs(sv.depth,EAF, size=1000)/sv.depth
            EAF_sim   = binom.rvs(int(sv.depth),EAF, size=1000)/sv.depth
            ab_sim    = sv.frequency/EAF_sim          # take ratio
            burden_CI = percentile(ab_sim, [2.5,97.5])
         except:
            print 'problem with EAF_sim: ', sv.depth, EAF
         

      result = {'mutation' : sv.mutation,
                'position' : 'chr%s:%d' %(sv.chr_, sv.position),
                'depth'    : '%d' %sv.depth,
                'AF_obs'   : '%0.2f' %sv.frequency,
                'AF_E[G]'  : [AF_G1, AF_G2],
                'AF_E[S]'  : [AF_S1, AF_S2],
                'purity'   : p,
                'CN'       : C,
                'M'        : M,
                'errMAF'   : error_MAF,
                'errLR'    : error_LR,
                'P(G)'     : max_prob_germline,
                'P(S)'     : max_prob_somatic,
                'log(G/S)' : logodds,
                'call'     : SG_status,
                'zygosity' : zygosity,
                'clonality': clonality,
                'in_tumor' : in_tumor,
                'burden'   : allele_burden,
                'burdenCI' : burden_CI}
         
      y.append(result)
   
   return y

def read_cna_model_file(fname):
   data_CNA = []
   with open(fname, 'r') as fin:
      reader = csv.DictReader(fin, dialect = 'excel-tab')
      for line in reader:
         seg = Segment(line['CHR'], line['segStart'], line['segEnd'], line['mafPred'], line['CN'], 
            line['segLR'], line['segMAF'], line['numMAtumorPred'], line['purity'], line['baseLevel'])
         data_CNA.append(seg)

   return data_CNA

def read_mut_aggr_full(fname):
   short_variants = []
   with open(fname, 'r') as fin:
      # if 'Notes' in fin:
         # next(fin)
      reader = csv.DictReader(fin, dialect = 'excel-tab')
      for line in reader:
         sv = SV(line['mutation'], line['frequency'], line['depth'], line['pos'])
         short_variants.append(sv)

   return short_variants

class Segment(object):
   """docstring for Segment"""
   def __init__(self, chr_, segStart, segEnd, mafPred, CN, segLR, segMAF, numMAtumorPred, purity, baseLevel):
      super(Segment, self).__init__()
      self.chr_           = chr_.strip("chr")
      self.segStart       = int(segStart)
      self.segEnd         = int(segEnd)
      self.mafPred        = float(mafPred) if mafPred != 'NA' else nan

      self.CN             = int(float(CN)) if CN != 'NA' else nan
      self.segLR          = float(segLR) if segLR != 'NA' else nan
      self.segMAF         = float(segMAF) if segMAF != 'NA' else nan
      self.numMAtumorPred = float(numMAtumorPred) if numMAtumorPred != 'NA' else nan
      self.purity         = float(purity)
      self.baseLevel      = float(baseLevel)
      
class SV(object):
   """lightweight structure for short variants"""
   def __init__(self, mutation, frequency, depth, position):
      super(SV, self).__init__()
      self.mutation  = mutation
      self.frequency = float(frequency)
      self.depth     = float(depth)
      
      self.chr_ = position.split(':')[0].strip('chr')
      self.position  = int(position.split(':')[1])
      

def _arg_parser():
   """
   Handle input arguments.

   Returns:
     Argparse parser object.
   """
   script_version = globals().get('__version__')
   script_description = globals().get('__doc__')
   script_epilog = None
   script_usage = '''%(prog)s [options] aggregated_mutations_file cna_model_file
                 %(prog)s [-h|--help]
                 %(prog)s [--version]'''

   parser = argparse.ArgumentParser(usage=script_usage,
                                  description=script_description,
                                  epilog=script_epilog,
                                  add_help=False)

   g = parser.add_argument_group('Program Help')
   g.add_argument('-h', '--help', action='help',
                help='show this help message and exit')
   g.add_argument('--version', action='version', version=script_version)

   g = parser.add_argument_group('Required Arguments')
   g.add_argument('aggregated_mutations_file', action='store', type=str,
                help='Full aggregated mutations text file')
   g.add_argument('cna_model_file', action='store', type=str,
                help='CNA model calls')

   g = parser.add_argument_group('Optional Arguments')
   g.add_argument('-o', dest='output_header', action='store', type=str, default='',
                help='Output file folder and header')

   return parser

def main(args):
   
   random.seed(10)
   fname_muts    = args.aggregated_mutations_file  # vars/sample11.mut_aggr.full.txt
   fname_CNmodel = args.cna_model_file
   out_header = args.output_header

   if len(out_header)==0:
      out_header = fname_muts[0:-4]
   sgz_full = out_header + '.fmi.sgz.full.txt'
   sgz_out = out_header + '.fmi.sgz.txt'

   # Sometimes we can not calculate a CNA model, so an empty file is supplied.  We need
   # to ensure that this script will produce its anticipated output file in this case.
   if not os.path.getsize(fname_CNmodel):
       logger.warn("Empty cna_calls.txt file.  Not running SGZ code.")
       for filename in (sgz_full, sgz_out):
           open(filename, 'w').close()
       sys.exit(0)
   else:
      #############################################################
      # 1) Read CN model file
      data_CNA = read_cna_model_file(fname_CNmodel)
      # 2) Read mutations file
      short_variants = read_mut_aggr_full(fname_muts)

      # 3) Compute SGZ
      y_obj = core_SGZ(data_CNA, short_variants)

      # 4) Write output
      fout      = open(sgz_full, 'w')
      fout2     = open(sgz_out, 'w')   

      header = ['mutation', 'pos', 'depth', 'frequency', 'afG1', 'afS1', 'afG2', 'afS2',
                'p', 'C', 'M', 'logOR_G', 'clonality', 'clonality_CI_low', 'clonality_CI_high',
                'germline/somatic', 'zygosity']
      fout.write('\t'.join(header) + '\n')
      
   
      header = ['mutation', 'pos', 'depth', 'frequency', 'C', 'germline/somatic', 'zygosity']
      fout2.write('\t'.join(header) + '\n')
   
      for sv in y_obj:
         
         out_list = [sv['mutation'],
                     sv['position'],
                     sv['depth'],
                     sv['AF_obs'],
                     '%.2f' % sv['AF_E[G]'][0],
                     '%.2f' % sv['AF_E[S]'][0],
                     '%.2f' % sv['AF_E[G]'][1],
                     '%.2f' % sv['AF_E[S]'][1],
                     '%.3f' % sv['purity'],
                     '%g'   % sv['CN'],
                     '%g'   % sv['M'],
                     '%.1f' % sv['log(G/S)'],
                     '%.2f' % sv['burden'], 
                     '%.2f' % sv['burdenCI'][0], 
                     '%.2f' % sv['burdenCI'][1], 
                     sv['call'],
                     sv['zygosity']]
      
         fout.write('\t'.join(out_list)+'\n')
      

         out_list = [sv['mutation'],
                     sv['position'],
                     sv['depth'],
                     sv['AF_obs'],
                     '%g'   % sv['CN'],
                     sv['call'],
                     sv['zygosity']]
         
         fout2.write('\t'.join(out_list)+'\n')


      fout.close()
      fout2.close()
      
if __name__ == "__main__":
   args = None
   logger = None

   logging.basicConfig(format='%(levelname)s: [%(asctime)s] %(message)s',
                     datefmt='%Y-%m-%d %H:%M:%S %Z')
   logger = logging.getLogger()
   logger.setLevel(logging.INFO)

   # Here, we collect the arguments into a dictionary keyed on argname.
   args = _arg_parser().parse_args(sys.argv[1:])

   if 'debug' in args:
      logger.setLevel(logging.DEBUG)

   sys.exit(main(args))

