"""
"""

import os
import sys
import shutil
import tempfile
import subprocess

class ClustalWrun:
    """
    """
    def __init__(self,opts=None):
        self.opts = opts
        self.iname = 'infile_copy'
        shutil.copy(self.opts.input,self.iname) 

    def run(self):
        tlf = open(self.opts.outlog,'w')
        cl = ['clustalw2 -INFILE=%s -OUTFILE=%s -OUTORDER=%s -TYPE=%s -OUTPUT=%s' % (self.iname,self.opts.output,self.opts.out_order,self.opts.dnarna,self.opts.outform)]
        if self.opts.seq_range_end <> None and self.opts.seq_range_start <> None:
            cl.append('-RANGE=%s,%s' % (self.opts.seq_range_start,self.opts.seq_range_end))
        if self.opts.outform=='CLUSTAL' and self.opts.outseqnos <> None:
            cl.append('-SEQNOS=ON')
        process = subprocess.Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf)
        rval = process.wait()
    

proc = ClustalWrun(opts)
proc.run()
    
            

