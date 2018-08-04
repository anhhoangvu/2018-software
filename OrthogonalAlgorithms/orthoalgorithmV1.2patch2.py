# orthoribalalgorithm.py
# Determine Shine Dalgarno (SD) and Anti-Shine Dalgarno sequences for strains of bacteria to design orthogonal ribosomes
# Rice University iGEM 2018
# Vu Hoang Anh and Sai Sriram
# Correspondence: ss162@rice.edu; ahv1@rice.edu


# ---------------------------------------Step 1----------------------------------------------------------
# -------------------------------- Initialize library --------------------------------------------------


def libbuild(rRNA): 
    """
    Build the initial library of 4096 pairs including
    all possible combinations from randomizing 6 bases in the SD/ASD region.
    Additionally, formats the initial dictionary so that the first SD has the key "SD1", and so on.
    :param: the 16s rRNA sequence of the bacteria strain
    :return: The initial dictionary containing all possible SD/ASD sequences
    """
    wtASD = getASD(rRNA)[0]
    # Build a library of 4096 pairs
    Lib = []
    basecode = ['A','U','C','G']
    for b1 in range(1, 5):
        for b2 in range(1, 5):
            for b3 in range(1, 5):
                for b4 in range(1, 5):
                    for b5 in range(1, 5):
                        for b6 in range(1, 5):
                            seqcode = str(b1) + str(b2) + str(b3) + str(b4) + str(b5) + str(b6)
                            SDseq = basecode[int(seqcode[0])-1] + basecode[int(seqcode[1])-1] + basecode[int(seqcode[2])-1] + basecode[int(seqcode[3])-1] + basecode[int(seqcode[4])-1] + basecode[int(seqcode[5])-1]
                            ASDseq = wtASD[:4] + revcomp(SDseq) + wtASD[10:]
                            Lib.append([SDseq, ASDseq])
    Libdict = {}
    for x in range(0, 4096):
        pair = Lib[x]
        Libdict['SD{0}'.format(x + 1)] = pair[0]
        Libdict['ASD{0}'.format(x + 1)] = pair[1]
    return Libdict
    return IniLibrary


# --------------------------------- Supporting functions -------------------------------------------------------

def convertTtoU(sequence):
    """
    Convert thymine to uracil
    :param sequence: The DNA sequence that is being converted to mRNA, in string form
    :return: Another string with T's replaced with U's.
    """
    return sequence.replace('T', 'U')


def revcomp(sequence):
    """
    Find reverse complementary sequence
    :param sequence: The RNA sequence in string form
    :return: The reverse complement sequence in string form
    """
    complement = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}
    revcompseq = ""
    sequence_list = list(sequence)
    sequence_list.reverse()
    for letter in sequence_list:
        revcompseq += complement[letter.upper()]
    return revcompseq


def getASD(rRNA):
    """
    :param rRNA: the 16s rRNA sequence
    :return: ASD: the ASD region
    :return: index: the index of the first nucleotide in the ASD region
    """
    import re
    for pseudoindex in re.finditer('ACCUCC',rRNA):
        index = pseudoindex.start()
    index = index - 3
    ASD = rRNA[index:index+12]
    return ASD,index


def revcompDNA(sequence):
    """
    Find reverse complementary sequence
    :param sequence: The RNA sequence in string form
    :return: The reverse complement sequence in string form
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    revcompseq = ""
    sequence_list = list(sequence)
    sequence_list.reverse()
    for letter in sequence_list:
        revcompseq += complement[letter.upper()]
    return revcompseq


def RNAduplexval(firstseq, secondseq):
    """
    Call RNAduplex in python shell; change RNAduplex path used below corresponding to position of RNAduplex.exe
    Pipes results from RNAduplex back into the algorithm
    :param firstseq: The ASD sequence as a string
    :param secondseq: The SD sequence as a string
    :return: A value corresponding to the binding energy of the ASD and SD in kcal/mol
    """
    import subprocess
    import os
    input = firstseq + '\n' + secondseq + '\n'
    input = input.encode('utf-8')
    from subprocess import Popen, PIPE, STDOUT
    RNAduplexpath = os.path.join(os.path.dirname(__file__), "Vienna\\RNAduplex.exe")
    result = Popen([RNAduplexpath], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    n = val.find(':')
    val = val[n + 11:-3]
    return val


def updatelibindex(library):
    """
    Update the index of the sequence in the library to be continuous
    After narrowing down the dictionary through the various steps of the algorithm, certain terms are deleted. Thus, the
    dictionary is no longer continuous, and could be ordered "ASD1" and the next term "ASD9". This function changes the
    "ASD9" into an "ASD2"
    :param library: The dictionary whose keys need to be updated.
    :return: The dictionary with updated keys.
    """
    liblen = len(library) / 2
    n = 1
    while (n < liblen + 1):
        testidx = n
        while 'ASD' + str(testidx) not in library:
            testidx = testidx + 1
        library['ASD{0}'.format(n)] = library.pop('ASD{0}'.format(testidx))
        library['SD{0}'.format(n)] = library.pop('SD{0}'.format(testidx))
        n = n + 1
    return library


def RNAfoldseq(fullRNAsequence):
    """
    Call RNAfold program in python using the default option utilizing MFE structure.
    :param fullRNAsequence: The RNA sequence in string form.
    :return: The free energy of the folded structure.
    """
    import subprocess
    import os
    input = fullRNAsequence
    input = input.encode('utf-8')
    from subprocess import Popen, PIPE, STDOUT
    RNAfoldpath = os.path.join(os.path.dirname(__file__), "Vienna\\RNAfold.exe")
    result = Popen([RNAfoldpath], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    val = val[len(fullRNAsequence) + 1:len(fullRNAsequence) * 2 + 1]
    return val


def RNAfoldcentroidseq(fullRNAsequence):
    """
    Call RNAfold program in python using the centroid option.
    :param fullRNAsequence: The RNA sequence in string form.
    :return: The free energy of the folded structure.
    """
    import subprocess
    import os
    input = fullRNAsequence
    input = input.encode('utf-8')
    from subprocess import Popen, PIPE, STDOUT
    RNAfoldpath = os.path.join(os.path.dirname(__file__), "Vienna\\RNAfold.exe")
    result = Popen([RNAfoldpath, "-p"], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    index = val.find('d')
    val = val[index - 10 - len(fullRNAsequence):index - 10]
    return val


def formatCDSfile(CDSfile, nameoforganism, path):
    """
    Formats the CDS file to be usable with the algorithm. This function creates a new .txt file with the new CDS file
    :param CDSfile: The file path of the CDS file from the NCBI database
    :param nameoforganism: The name of the organism for the purpose of creating the file path name
    :param path: The path where you want the formatted file to end up. If I wanted my file to end up in a folder called
    "Python Code", I would use C:/users/siddu/Desktop/Python Code/
    :return: filepath, the path of the new file with the formatted CDS.
    """

    with open(CDSfile, 'r+') as cdsfile:
        filepath = path + nameoforganism + 'CDS_formatted.txt'
        with open(filepath, 'w') as formatted:
            # ---------------------- Format the CDS file to be usable with the algorithm-----------------------
            previousline = ''
            for line in cdsfile:
                currentandprevious = previousline.rstrip() + line.rstrip()
                idx1 = currentandprevious.find('location=')
                if idx1 == -1:  # gets rid of all lines that don't actually matter
                    pass
                else:
                    if currentandprevious[idx1 + 9] == 'c':  # complement cases
                        idx2 = line.find(')]')
                        if idx2 > idx1:
                            line = line[0: idx2 + 1] + '\n'
                            formatted.write(line)
                        else:  # some cases have ')]' (what we are using to localize the CDS in each line) repeated, so
                            # this treats them by ignoring the first instance.
                            line.replace(')]', 'ss', 1).find(')]')
                            line = line[0: idx2 + 1] + '\n'
                            formatted.write(line)
                    elif currentandprevious[idx1 + 9] == 'j':  # Join cases are ignored
                        pass
                    else:
                        idx2 = line.find('..')
                        if idx2 > idx1:
                            line = line[0: idx2] + '\n'
                            formatted.write(line)
                        else:
                            line.replace('..', 'ss', 1).find('..')
                            line = line[0: idx2] + '\n'
                            formatted.write(line)
    return filepath


def findCDSfromformatted(filepath):
    """
    Gives the lists of CDSs in a particular genome given the formatted filepath from formatCDSfile.
    :param filepath: The filepath of the formatted CDS file from the function formatCDSfile
    :return: cdslist1 - The list of CDSs corresponding to the cases where the CDS is on the genome strand
    :return: cdslist2 - The list of CDSs corresponding to when the CDS is on the complement strand.
    """
    with open(filepath, 'r') as cdsfile:  # get all CDSs
        # Note that the CDS file sometimes skips some numbers -- for example goes from 54 to 56.
        cdslist1 = []
        cdslist2 = []
        for line in cdsfile:
            idx1 = line.find('[location=')
            if idx1 != -1:
                idx2 = line.find('..')
                if idx2 < idx1:
                    idx2 = line.replace('..', 'ss', 1).find('..')  # in the case that there is another .. in the line,
                    # we replace it with ss then find the actual '..' to localize the CDS.
                if line[idx1 + 10: idx1 + 13] == 'com':
                    # complement cases -- take the second number then reverse complement
                    if line[idx1 + 21] == 'j':  # join cases -- unaccounted for
                        pass
                    else:
                        idx2_1 = line.find(')]\n')
                        cdslist2.append(int(line[idx2 + 2: idx2_1 - 1]) - 1)  # -1 to account for Python indexing
                        pass
                else:
                    cdslist1.append(int(line[idx1 + 10: idx2]) - 1)
            else:
                pass
    return cdslist1, cdslist2


def getallTIRs(CDSfile, genome, nameoforganism, path):
    """
    Builds a dictionary of strings that contains all the TIRS, accessible via
    the TIR number.
    eg TIRdict['TIR1'] returns the first TIR in the genome. TIRs are determined from the CDS file of the species as well
    as a .txt of the genome file, containing only base pairs.
    :param CDSfile: The file path to the the species CDS file from NCBI
    :param genome: The file path to the species genome file, which is purely the genome itself.
    :param nameoforganism: The name of the organism
    :param path: The file path where you want the formatted CDS to end up (e.g. 'C:/users/siddu/Desktop/Python Code/').
    :return: A dictionary containing all translation initiation regions, starting 20 bps in front of the start codon and
    including the start codon itself. Formatted to access the first TIR of the genome by key "TIR1"
    """
    filepath = formatCDSfile(CDSfile, nameoforganism, path)
    cdslist1, cdslist2 = findCDSfromformatted(filepath)
    print('Number of CDS = ' + str(len(cdslist1)))
    print('Number of complement CDS = ' + str(len(cdslist2)))
    # MUST ACCOUNT FOR THE JOIN CASES
    # ------------------------ Regular cases when TIRs are found on the genome strand ---------------------------------
    with open(genome, 'r') as genomefile:
        index = 0
        i = 0
        previousline = ''
        TIRdict = {}
        for line in genomefile:  # cdslist1 deals with the forward sequences
            line = line.rstrip()
            currentandprevious = previousline.rstrip() + line.rstrip()
            # we concatenate the current and previous lines to account for
            # CDS that overlap over two lines
            if i < len(cdslist1):
                idx3 = cdslist1[i] - index
                if idx3 > len(line) - 4:  # 66
                    # this is the case when the end 3 bps overlap into the next line or we haven't yet reached
                    # the line of the CDS yet
                    previousline = line
                    index += len(line)  # each line is 70 bps long
                    pass

                else:
                    # here we have the case when we already passed the start of this TIR, so we have to go back to the
                    # previous line
                    a = convertTtoU(currentandprevious[idx3 + (len(line) - 21):idx3 + len(line) + 3])  # 49, 73

                    TIRdict['TIR' + str(i + 1)] = a  # 59 to 73
                    i = i + 1
                    previousline = line.rstrip()
                    index += len(line)  # each line is 70 bps long
                    pass

    # ------------------ Complement cases where TIRs are found on complementary strand -------------------
    with open(genome, 'r') as genomefile:
        index = 0
        j = i
        i = 0
        previousline = ''
        for line in genomefile:  # This loop handles cdslist2, which contains the reverse complement sequences.
            line = line.rstrip()
            currentandprevious = previousline.rstrip() + line.rstrip()
            # we concatenate the current and previous lines to account for
            # CDS that overlap over two lines
            if i < len(cdslist2):
                idx3 = cdslist2[i] - index
                if idx3 > len(line) - 21:
                    # this is the case when the end 3 bps overlap into the next line or we haven't yet reached
                    # the line of the CDS yet
                    previousline = line.rstrip()
                    index += len(line)  # each line is 70 bps long
                    pass
                else:
                    # here we have the case when we already passed the start of this TIR, so we have to go back to the
                    # previous line
                    a = convertTtoU(revcompDNA(currentandprevious[idx3 + (len(line) - 2):idx3 + (len(line) + 22)]))
                    TIRdict['TIR' + str(i + 1 + j)] = a
                    i = i + 1
                    previousline = line.rstrip()
                    index += len(line)  # each line is 70 bps long
                    pass
    return TIRdict


def secstructurelib(Library,filename,start16sseq,stop16sseq,direction):
    """
    Build a library of secondary structure of 16s rRNA corresponding to a given library of ASD and SD sequence.
    :param Library: The library of ASD SD sequence
    :param: filename: the name of the file containing the full genome sequence
    :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
    :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
    :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
    :return: A library containing secondary structure of rRNA with each ASD
    """
    full16srRNAseq = get16srRNAseq(filename,start16sseq,stop16sseq,direction)
    NewLib = {}
    index = getASD(full16srRNAseq)[1]
    for n in range(0, int(len(Library) / 2)):
        full16srRNAseq = full16srRNAseq[0:index] + Library['ASD{0}'.format(n + 1)] + full16srRNAseq[index+12:]
        secstructure = RNAfoldcentroidseq(full16srRNAseq)
        NewLib['foldedseq{0}'.format(n + 1)] = str(secstructure)
    return NewLib


def get16srRNAseq(filename,start16sseq,stop16sseq,direction):
    """
    Get the full 16s rRNA sequence from the genome
    :param: filename: the name of the file containing the full genome sequence
    :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
    :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
    :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
    :return: the full 16s rRNA sequence from the genome
    """
    import os
    filepath = os.path.join(os.path.dirname(__file__), filename)
    with open(filepath, 'r') as myfile:  # determine rRNA sequence
        data = myfile.read().replace('\n', '')
    full16srRNAseq = convertTtoU(data[start16sseq-1:stop16sseq])
    if direction == 'reverse':
        listfull16srRNAseq = list(reversed(full16srRNAseq))
        pseudosequence = ''
        for i in listfull16srRNAseq:
            pseudosequence += str(i)
        full16srRNAseq = pseudosequence
    return full16srRNAseq


# ---------------------- Define a class to model the library and apply each step to that object -----------------


class Library:
    """A model for the library of ASD/ SD sequences that will be narrowed as more steps are applied"""

    def __init__(self, library):
        """
        :param genome: The string corresponding to the file path with the strain's genome file
        :param library: The library itself
        """
        # ******************* Still have to make self.genome a parameter by moving get16srRNAseq into the class
        self.library = library

    # -----------------------------------------Steps 2-3----------------------------------------------------------------

    def narrow_binding(self,filename,start16sseq,stop16sseq,direction):
        """
        Input the wildtype ASD and SD, output is the narrowed library with pairs that fall within 0.5 kcal range of the
        wildtype binding energy
        :param: filename: the name of the file containing the full genome sequence
        :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
        :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
        :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
        :return: The narrowed dictionary with updated key indices only consisting of ASD/ SD pairs that have binding
        energies of -0.5 <= x <= 0.5
        """
        rRNA = get16srRNAseq(filename,start16sseq,stop16sseq,direction)

        WtASD = getASD(rRNA)[0]
        WtSD = revcomp(WtASD[4:10])
        Wildval = RNAduplexval(WtASD, WtSD)  # Calculate the binding energies of Wild-type ASD/SD pair
        # self.library is the initial library after step 1
        for n in range(0, 4096):
            ASD = self.library['ASD' + str(n + 1)]
            SD = self.library['SD' + str(n + 1)]
            val = RNAduplexval(ASD, SD)
            if float(Wildval) - 0.5 >= float(val) or float(val) >= float(Wildval) + 0.5:  # Compare each ASD/SD binding
                # energy to wild-type binding energy
                del self.library['ASD' + str(n + 1)]
                del self.library['SD' + str(n + 1)]
        self.library = updatelibindex(self.library)
        print('Library after step 3: ', self.library)
        return self.library

    # -------------------------------Steps 4-5--------------------------

    def narrow_crossbinding(self,filename,start16sseq,stop16sseq,direction):
        """
        Narrow down the library more by eliminating pairs with ASD that has < 1 kcal/mol binding energy with
        wild type SD
        :param: filename: the name of the file containing the full genome sequence
        :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
        :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
        :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
        :return: Further narrowed dictionary.
        """
        rRNA = get16srRNAseq(filename,start16sseq,stop16sseq,direction)
        WtASD = getASD(rRNA)[0]
        WtSD = revcomp(WtASD[4:10])
        for n in range(0, int(len(self.library) / 2)):
            ASD = 'ASD' + str(n + 1)
            SD = 'SD' + str(n + 1)
            if float(RNAduplexval(self.library[ASD], WtSD)) < -1:
                del self.library[ASD]
                del self.library[SD]
        self.library = updatelibindex(self.library)
        print('Library after step 5:', self.library)
        return self.library

    # --------------------------Steps 6-7-----------------------------------

    def ASD_2rystructure_narrow(self,filename,start16sseq,stop16sseq,direction):
        """ Narrow down the library by discarding sequences that forms secondary structure in ASD region
         Import the most narrowed Libdict up to step 5 as well as the secondary structure dictionary.
         :param: filename: the name of the file containing the full genome sequence
         :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
         :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
         :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
         :return: The narrowed dictionary of ASD/ SD that do not have secondary structure in the ASD region.
         """

        full16srRNAseq = get16srRNAseq(filename,start16sseq,stop16sseq,direction)
        secstructureseqlib = secstructurelib(self.library,filename,start16sseq,stop16sseq,direction)
        # Locate the index of the first bp of the ASD for all 16S rRNA sequences
        index = getASD(rRNA)[1]
        """Iterate through the whole secondary structure dictionary, locating a constant positioned
        ASD and determining whether there is secondary structure in that ASD"""
        for i in range(0, int(len(self.library) / 2)):
            secname = "foldedseq" + str(i + 1)
            ASDname = "ASD" + str(i + 1)
            SDname = "SD" + str(i + 1)
            secstructureseq = secstructureseqlib[secname]
            ASDrandomregion = secstructureseq[index + 4:index + 10]
            if "(" in ASDrandomregion or ")" in ASDrandomregion or "{" in ASDrandomregion or "}" in ASDrandomregion or "[" in ASDrandomregion or "]" in ASDrandomregion:
                del self.library[ASDname]
                del self.library[SDname]
        self.library = updatelibindex(self.library)
        print('Library after step 7: ', self.library)
        return self.library

    # ---------------------------------Steps 8-10-----------------------------------------------------

    def allASDTIRpairs(self, CDSfile, genome, nameoforganism, path):
        """
        Iterate through all possible ASD TIR pairs and find the ones with the highest average binding energies with host
        TIRs
        We wish to find the ASDs that do not bind well with the host translation initiation regions to assure orthogonality
        of the ribosomes, so here we choose the ASDs that have the highest binding energies (i.e. don't bind well with host
        TIRs)
        :return: Prints the list of the top ten ASD candidates.
        """
        TIRdict = getallTIRs(CDSfile, genome, nameoforganism, path)
        dictofvals = {}
        print("Number of TIRs: " + str(len(TIRdict)))
        listofaverages = []
        for i in range(0, round(len(self.library) / 2)):  # iterate through all ASDs
            listofvals = []
            ASDname = str('ASD' + str(i + 1))
            for j in range(0, len(TIRdict)):  # for each ASD, iterate through all TIRs in the genome
                TIRname = str('TIR' + str(j + 1))
                if TIRdict[TIRname] == '':
                    pass
                else:
                    val = float(RNAduplexval(self.library[ASDname], TIRdict[TIRname]))
                    listofvals.append(val)
            average = sum(listofvals) / len(listofvals)
            dictofvals[average] = ASDname  # calculate the average binding energy between the ASD
            # and all TIRs; here we store the key as the average so that the we can call the names of the highest ASDs after
            # the list is sorted
            listofaverages.append(average)

        listofaverages.sort(reverse=True)
        print('Here are the 10 top candidates with highest ASD-host binding values:')
        for i in range(0, 10):
            print(self.library[str(dictofvals[listofaverages[i]])])

# # ------------------------ Create the E. coli instance and run all the steps on it----------------------
# ## 1542 bp long 16s rRNA, position 4,166,659 -> 4,168,200, postion of ASD from 1535-1540
# import pickle
# import os
# CDSfile = os.path.join(os.path.dirname(__file__), "E coli K-12 MG1655 CDSs.fasta.txt")
# genome = os.path.join(os.path.dirname(__file__), "E coli K-12 MG1655 genome.fasta.txt")
# nameoforganism = 'E. Coli'
# path = os.path.dirname(__file__)
# start = 4166659
# end = 4168200
# direction = 'forward'
# rRNA = get16srRNAseq(genome,start,end,direction)
# initiallib = libbuild(rRNA)
# ecolistrain = Library(initiallib) #create the E. Coli instance.
# ecolistrain.narrow_binding(genome,start,end)
# ecolistrain.narrow_crossbinding(genome,start,end)
# ecolistrain.ASD_2rystructure_narrow(genome,start,end)
# pickle.dump(ecolistrain, open("E.coli after step 7 ASD-SD lists.p",'wb'))
# ecolistrain.allASDTIRpairs(CDSfile, genome, nameoforganism, path)
#
# # ------------------------ Create the P. Putida instance and run all the steps on it----------------------
# ## 1522 bp long 16s rRNA, position 2,548,690 -> 2,550,211
# import os
# CDSfile = os.path.join(os.path.dirname(__file__), "Pseudomonas putida KT2440 (g-proteobacteria) CDS.txt")
# genome = os.path.join(os.path.dirname(__file__), "Pseudomonas putida KT2440 (g-proteobacteria) genome.txt")
# nameoforganism = 'P. Putida'
# path = os.path.dirname(__file__)
# start = 2548690
# end = 2550211
# direction = 'forward'
# rRNA = get16srRNAseq(genome,start,end,direction)
# initiallib = libbuild(rRNA)
# pputidastrain = Library(initiallib) #create the Pneudomonas Putida instance.
# pputidastrain.narrow_binding(genome,start,end)
# pputidastrain.narrow_crossbinding(genome,start,end)
# pputidastrain.ASD_2rystructure_narrow(genome,start,end)
# pickle.dump(pputidastrain, open("P. Putida after step 7 ASD-SD lists.p",'wb'))
# pputidastrain.allASDTIRpairs(CDSfile, genome, nameoforganism, path)
#
# # ------------------------ Create the L. Lactis instance and run all the steps on it----------------------
# ## 1548 bp long 16s rRNA, position 537,561 -> 539,108
# import os
# CDSfile = os.path.join(os.path.dirname(__file__), "Lactococcus lactis subsp. lactis Il1403 (firmicutes) CDS.txt")
# genome = os.path.join(os.path.dirname(__file__), "Lactococcus lactis subsp. lactis Il1403 (firmicutes) genome.txt")
# nameoforganism = 'L. Lactis'
# path = os.path.dirname(__file__)
# start = 537561
# end = 539108
# direction = 'forward'
# rRNA = get16srRNAseq(genome,start,end,direction)
# initiallib = libbuild(rRNA)
# llactisstrain = Library(initiallib) #create the L. Lactis instance.
# llactisstrain.narrow_binding(genome,start,end)
# llactisstrain.narrow_crossbinding(genome,start,end)
# llactisstrain.ASD_2rystructure_narrow(genome,start,end)
# pickle.dump(llactisstrain, open("L. Lactis after step 7 ASD-SD lists.p",'wb'))
# llactisstrain.allASDTIRpairs(CDSfile, genome, nameoforganism, path)
#
# # ------------------------ Create the B. Subtilis instance and run all the steps on it----------------------
# ## 1555 bp long 16s rRNA, position 9,810 -> 11,364
# import os
# CDSfile = os.path.join(os.path.dirname(__file__), "Bacillus subtilis str. 168 CDS.txt")
# genome = os.path.join(os.path.dirname(__file__), "Bacillus subtilis str. 168 genome.txt")
# nameoforganism = 'B. Subtilis'
# path = os.path.dirname(__file__)
# start = 9810
# end = 11364
# direction = 'forward'
# rRNA = get16srRNAseq(genome,start,end,direction)
# initiallib = libbuild(rRNA)
# bsubtilisstrain = Library(initiallib) #create the B. Subtilis instance.
# bsubtilisstrain.narrow_binding(genome,start,end)
# bsubtilisstrain.narrow_crossbinding(genome,start,end)
# bsubtilisstrain.ASD_2rystructure_narrow(genome,start,end)
# pickle.dump(bsubtilisstrain, open("B. Subtilis after step 7 ASD-SD lists.p",'wb'))
# bsubtilisstrain.allASDTIRpairs(CDSfile, genome, nameoforganism, path)
#
#

# # ------------------------ Create the A. Baylyi instance and run all the steps on it----------------------
# ## 1472 bp long 16s rRNA, position 193,231 <- 194,702
# import os
# import pickle
# CDSfile = os.path.join(os.path.dirname(__file__), "Acinetobacter baylyi CIP 107474 CDSs.txt")
# genome = os.path.join(os.path.dirname(__file__), "Acinetobacter baylyi CIP 107474 genome.txt")
# nameoforganism = 'A. Baylyi'
# path = os.path.dirname(__file__)
# start = 193231
# end = 194702
# direction = 'reverse'
# rRNA = get16srRNAseq(genome,start,end,direction)
# print(rRNA)
# initiallib = libbuild(rRNA)
# abaylyistrain = Library(initiallib) #create the A. Baylyi instance.
# abaylyistrain.narrow_binding(genome,start,end,direction)
# abaylyistrain.narrow_crossbinding(genome,start,end,direction)
# abaylyistrain.ASD_2rystructure_narrow(genome,start,end,direction)
# pickle.dump(abaylyistrain, open("A. Baylyi after step 7 ASD-SD lists.p",'wb'))
# abaylyistrain.allASDTIRpairs(CDSfile, genome, nameoforganism, path)
# 
