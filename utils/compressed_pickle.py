"""
save and load compressed pickled objects in python 
"""

import bz2 
import sys 
import cPickle

def save(filename, myObj):
    """
    save the object to a file using compressed pickle in bz2 

    @args filename: result file name 
    @type filename: str 
    @args myObj: object to save 
    @type myObj: python pickleable object
    """

    try:
        fh = bz2.BZ2File(filename, "wb")
    except IOError, details:
        sys.stderr.write('File %s cannot be written\n' % filename)
        sys.stderr.write(details) 
        return 

    cPickle.dump(myObj, fh, protocol=2) 
    fh.close()


def load(filename):
    """
    load a compressed pickled file in bz2 format

    @args filename: name of the file to load 
    @type filename: str
    """
    
    try:
        fh = bz2.BZ2File(filename, 'rb') 
    except IOError, details:
        sys.stderr.write('File %s cannot be read\n' % filename)
        sys.stderr.write(details) 
        return 

    myObj = cPickle.load(fh) 
    fh.close() 

    return myObj
        

