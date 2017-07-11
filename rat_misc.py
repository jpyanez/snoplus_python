import os, rat

def openRat(infile_dir='', infile_name=None, *ratreader):
    try:
        ratreader.close()
        print 'Closing file before reopenning'
    except:
        print 'No ratreader file. Opening it for the first time'
    if infile_name == None:
        full_name = infile_dir
    else:
        full_name = os.path.join(infile_dir, infile_name)
    if os.path.exists(full_name):
        ratreader = rat.dsreader(full_name)
    else:
        print 'File does not exist'
        print full_name
    return ratreader
