def RNAfoldseqcentroid(fullRNAsequence):
    # call RNAfold program in python using centroid , ***change full path of program accordingly
    import subprocess
    import os
    input = fullRNAsequence
    input = input.encode('utf-8')
    from subprocess import Popen,PIPE, STDOUT
    RNAfoldpath = os.path.join(os.path.dirname(__file__), "Vienna\\RNAfold.exe")
    result = Popen([RNAfoldpath, "-p"], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    index = val.find('d')
    val = val[index-10-len(fullRNAsequence):index-10]
    return val

print(RNAfoldseqcentroid('AUCUUUAUCUACUUAUCUAUCUAUCU'))