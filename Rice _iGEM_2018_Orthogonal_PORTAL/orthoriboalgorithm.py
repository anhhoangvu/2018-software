# orthoribalalgorithm.py
# Determine Shine Dalgarno (SD) and Anti-Shine Dalgarno sequences for strains of bacteria to design orthogonal ribosomes
# Rice University iGEM 2018
# Vu Hoang Anh and Sai Sriram
# Correspondence: ss162@rice.edu; ahv1@rice.edu

import threading
import time
import re
import subprocess
import os
from subprocess import Popen, PIPE, STDOUT
import pickle
import threading
import time
from tkinter import *
from PIL import ImageTk, Image
import queue
import random
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
                            ASDseq = wtASD[:3] + revcomp(SDseq) + wtASD[9:]
                            Lib.append([SDseq, ASDseq])
    Libdict = {}
    for x in range(0, 4096):
        pair = Lib[x]
        Libdict['SD{0}'.format(x + 1)] = pair[0]
        Libdict['ASD{0}'.format(x + 1)] = pair[1]
    return Libdict


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
    for pseudoindex in re.finditer('ACCUCC', rRNA):
        index = pseudoindex.start()
    index = index - 3
    ASD = rRNA[index:index+12]
    return ASD, index


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

    input = firstseq + '\n' + secondseq + '\n'
    input = input.encode('utf-8')

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


    input = fullRNAsequence
    input = input.encode('utf-8')

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

    input = fullRNAsequence
    input = input.encode('utf-8')

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
                if line[idx1 + 10: idx1 + 13] == 'com':
                    # complement cases -- take the second number then reverse complement
                    idx2 = line.find('..')
                    if idx2 < idx1:
                        idx2 = line.replace('..', 'ss', 1).find('..')  # in the case that there is another ..
                        # in the line, we replace it with ss then find the actual '..' to localize the CDS.

                    if line[idx1 + 21] == 'j':  # join cases -- unaccounted for
                        print('This CDS location was not used ' + (line[idx1 + 21:]))
                        pass
                    else:
                        idx2_1 = line.find(')\n')
                        try:
                            cdslist2.append(int(line[idx2 + 2: idx2_1]) - 1)  # -1 to account for Python indexing
                        except ValueError:
                            print('This CDS location was not used ' + (line[idx2 + 2:idx2_1 - 1]))
                            pass
                        pass
                else:
                    idx2 = line.find('\n')
                    try:
                        cdslist1.append(int(line[idx1 + 10: idx2]) - 1)
                    except ValueError:
                        print('This CDS location was not used ' + (line[idx2 + 2:]))
                        pass
                    pass
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
                    # here we have the case when we already passed the start of this TIR, so the entire TIR is
                    # encompassed in the string "currentandprevious"
                    a = convertTtoU(currentandprevious[idx3 + (len(line) - 21):idx3 + len(line) + 3])  # 49, 73

                    TIRdict['TIR' + str(i + 1)] = a  # 59 to 73
                    i = i + 1
                    previousline = line.rstrip()
                    index += len(line)
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
    # In some cases, the TIRdict contains empty TIRs, so this removes them.
    error_keys = []
    for key in TIRdict:
        if TIRdict[key] == '':
            error_keys.append(key)
            print("Error: " + str(key) + " was not used due to a TIR that wasn't found.")
    for errorkey in error_keys:
        del TIRdict[errorkey]

    return TIRdict


def secstructurelib(Library, filename, start16sseq, stop16sseq, direction):
    """
    Build a library of secondary structure of 16s rRNA corresponding to a given library of ASD and SD sequence.
    :param Library: The library of ASD SD sequence
    :param: filename: the name of the file containing the full genome sequence
    :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
    :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
    :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
    :return: A library containing secondary structure of rRNA with each ASD
    """
    full16srRNAseq = get16srRNAseq(filename, start16sseq, stop16sseq, direction)
    NewLib = {}
    index = getASD(full16srRNAseq)[1]
    for n in range(0, int(len(Library) / 2)):
        partial16srRNAseq = full16srRNAseq[index-100:index] + Library['ASD{0}'.format(n + 1)] + full16srRNAseq[index+12:]
        # Assembles 16s rRNA sequence by adding the ASD sequence in the middle. Only goes 100 bps in front of ASD
        secstructure = RNAfoldcentroidseq(partial16srRNAseq)
        NewLib['foldedseq{0}'.format(n + 1)] = str(secstructure)
    return NewLib


def get16srRNAseq(filename, start16sseq, stop16sseq, direction):
    """
    Get the full 16s rRNA sequence from the genome
    :param: filename: the name of the file containing the full genome sequence
    :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
    :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
    :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
    :return: the full 16s rRNA sequence from the genome
    """

    filepath = os.path.join(os.path.dirname(__file__), filename)
    with open(filepath, 'r') as myfile:  # determine rRNA sequence
        data = myfile.read().replace('\n', '')
    if direction == "forward":
        full16srRNAseq = convertTtoU(data[start16sseq-1:stop16sseq])
    if direction == 'reverse':
        full16srRNAseq = convertTtoU(data[stop16sseq - 1:start16sseq])
        listfull16srRNAseq = list(revcomp(full16srRNAseq))
        pseudosequence = ''
        for i in listfull16srRNAseq:
            pseudosequence += str(i)
        full16srRNAseq = pseudosequence
    return full16srRNAseq


def split_up_threadwork(total_threads, thread_number, num_iterations):
    """
    Splits up a specified number of iterations into approximately equal workloads for multiple threads to work on.
    :param total_threads:
    :param thread_number:
    :param num_iterations:
    :return: the iteration numbers that the specific thread number will be working on (e.g. thread 1 will work on
    iterations 0-1024 for 4096 iterations split into 4 total threads.)
    """
    split_thread = round(num_iterations / total_threads)
    if thread_number == 1:
        index1 = 0
        index2 = (thread_number * split_thread) - 1
    elif thread_number != total_threads:  # not the last thread
        not_impt, prev_idx2 = split_up_threadwork(total_threads, thread_number - 1, num_iterations)
        # recursive definition to find the previous thread's end iteration, so the current thread's 1st iteration is
        # just the iteration after that.
        index1 = prev_idx2 + 1
        index2 = index1 + split_thread
    else:
        not_impt, prev_idx2 = split_up_threadwork(total_threads, thread_number - 1, num_iterations)
        index1 = prev_idx2 + 1
        index2 = num_iterations
    if thread_number == total_threads:
        index2 -= 1
    return index1, index2


# ---------------------- Define a class to model the library and apply each step to that object -----------------


class Library:
    """A model for the library of ASD/ SD sequences that will be narrowed as more steps are applied"""

    def __init__(self, library):
        """
        Initialize
        """
        self.library = library

    def deleter_fun(self):
        """
        Deletes all the terms from the library that are marked for deletion with the value set to "DELETE"
        :param library: The library containing terms that need to be deleted
        :return: The library with all the terms deleted.
        """
        for i in range(0, int(len(self.library)/2)):
            if self.library["ASD" + str(i+1)] == "DELETE":
                del self.library["ASD" + str(i+1)]
                del self.library["SD" + str(i+1)]

    def updatelibindex(self):
        """
        Update the index of the sequence in the library to be continuous
        After narrowing down the dictionary through the various steps of the algorithm, certain terms are deleted. Thus, the
        dictionary is no longer continuous, and could be ordered "ASD1" and the next term "ASD9". This function changes the
        "ASD9" into an "ASD2"
        :param library: The dictionary whose keys need to be updated.
        :return: The dictionary with updated keys.
        """
        liblen = len(self.library) / 2
        n = 1
        while n < liblen + 1:
            testidx = n
            while 'ASD' + str(testidx) not in self.library:
                testidx = testidx + 1
            self.library['ASD{0}'.format(n)] = self.library.pop('ASD{0}'.format(testidx))
            self.library['SD{0}'.format(n)] = self.library.pop('SD{0}'.format(testidx))
            n = n + 1
        return self.library

    # -----------------------------------------Steps 2-3----------------------------------------------------------------

    def narrow_binding(self, filename, start16sseq, stop16sseq, direction, total_threads, thread_number):
        """
        Input the wildtype ASD and SD, output is the narrowed library with pairs that fall within 0.5 kcal range of the
        wildtype binding energy
        :param: filename: the name of the file containing the full genome sequence
        :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
        :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
        :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
        :param: thread_number: the thread number that will be started for parallel computation
        :param: total_threads: the total number of threads that will be run in parallel
        :return: The narrowed dictionary with updated key indices only consisting of ASD/ SD pairs that have binding
        energies of -0.5 <= x <= 0.5
        """
        rRNA = get16srRNAseq(filename, start16sseq, stop16sseq, direction)
        outputbox.insert(END, "Starting Algorithm on thread " + str(thread_number) + "...\n")
        WtASD = getASD(rRNA)[0]
        WtSD = revcomp(WtASD[4:10])
        Wildval = RNAduplexval(WtASD, WtSD)  # Calculate the binding energies of Wild-type ASD/SD pair
        # self.library is the initial library after step 1
        index1, index2 = split_up_threadwork(total_threads, thread_number, 4096)
        for n in range(index1, index2):
            ASD = self.library['ASD' + str(n + 1)]
            SD = self.library['SD' + str(n + 1)]
            val = RNAduplexval(ASD, SD)
            if float(Wildval) - 0.5 >= float(val) or float(val) >= float(Wildval) + 0.5:  # Compare each ASD/SD binding
                # energy to wild-type binding energy
                self.library['ASD' + str(n + 1)] = "DELETE"
                self.library['SD' + str(n + 1)] = "DELETE"
        # print('Library after step 3: ', self.library)
        outputbox.insert(END, "Thread " + str(thread_number) + " in Step 3 is done\n")
        return self.library

    # -------------------------------Steps 4-5--------------------------

    def narrow_crossbinding(self, filename, start16sseq, stop16sseq, direction, total_threads, thread_number):
        """
        Narrow down the library more by eliminating pairs with ASD that has < -1 kcal/mol binding energy with
        wild type SD
        :param: filename: the name of the file containing the full genome sequence
        :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
        :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
        :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
        :return: Further narrowed dictionary.
        """
        rRNA = get16srRNAseq(filename, start16sseq, stop16sseq, direction)
        WtASD = getASD(rRNA)[0]
        WtSD = revcomp(WtASD[4:10])
        index1, index2 = split_up_threadwork(total_threads, thread_number, int(len(self.library) / 2))
        for n in range(index1, index2):
            ASD = 'ASD' + str(n + 1)
            SD = 'SD' + str(n + 1)
            if float(RNAduplexval(self.library[ASD], WtSD)) < -1:
                self.library[ASD] = "DELETE"  # Marking it for deletion rather than deleting it outright so that the
                # indexing for the next iteration is not messed up
                self.library[SD] = "DELETE"
        # print('Library after step 5:', self.library)
        outputbox.insert(END, "Thread " + str(thread_number) + " in Step 5 is done\n")
        return self.library

    # --------------------------Steps 6-7-----------------------------------

    def ASD_2rystructure_narrow(self, filename, start16sseq, stop16sseq, direction): #, total_threads, thread_number):
        """ Narrow down the library by discarding sequences that forms secondary structure in ASD region
         Import the most narrowed Libdict up to step 5 as well as the secondary structure dictionary.
         :param: filename: the name of the file containing the full genome sequence
         :param: start16sseq: the position of the first base of the 16s rRNA in the genome.
         :param: stop16sseq: the position of the last base of the 16s rRNA in the genome.
         :param: direction: indicate whether the 16s rRNA is translated in the forward of reverse direction.
         :return: The narrowed dictionary of ASD/ SD that do not have secondary structure in the ASD region.
         """

        full16srRNAseq = get16srRNAseq(filename, start16sseq, stop16sseq, direction)
        secstructureseqlib = secstructurelib(self.library, filename, start16sseq, stop16sseq, direction)
        # Locate the index of the first bp of the ASD for all 16S rRNA sequences
        print("SECSTRUCLIB ACQUIRED")
        index = getASD(full16srRNAseq)[1]

        """Iterate through the whole secondary structure dictionary, locating a constant positioned
        ASD and determining whether there is secondary structure in that ASD"""
        for i in range(0, int(len(self.library) / 2)):
            secname = "foldedseq" + str(i + 1)
            ASDname = "ASD" + str(i + 1)
            SDname = "SD" + str(i + 1)
            secstructureseq = secstructureseqlib[secname]
            ASDrandomregion = secstructureseq[index + 4:index + 10]
            if "(" in ASDrandomregion or ")" in ASDrandomregion or "{" in ASDrandomregion or "}" in ASDrandomregion or "[" in ASDrandomregion or "]" in ASDrandomregion:
                self.library[ASDname] = "DELETE"
                self.library[SDname] = "DELETE"
        outputbox.insert(END,  "Step 7 is done\n")
        return self.library

    # ---------------------------------Steps 8-10-----------------------------------------------------

    def allASDTIRpairs(self, TIRdict, total_threads, thread_number):
        """
        Iterate through all possible ASD TIR pairs and find the ones with the highest average binding energies with host
        TIRs
        We wish to find the ASDs that do not bind well with the host translation initiation regions to assure orthogonality
        of the ribosomes, so here we choose the ASDs that have the highest binding energies (i.e. don't bind well with host
        TIRs)
        :return: Prints the list of the top ten ASD candidates.
        """
        if thread_number == 1:
            print("The following is the number of ASDs that we are testing:")
            print(len(self.library)/2)

        index1, index2 = split_up_threadwork(total_threads, thread_number, (len(self.library) / 2))
        for i in range(int(index1), int(index2+1)):  # iterate through all ASDs
            outputbox.insert(END, "Thread on ASD number: " + str(i+1) + "\n")
            listofvals = []
            ASDname = str('ASD' + str(i + 1))
            for j in range(0, len(TIRdict)):  # for each ASD, iterate through all TIRs in the genome
                TIRname = 'TIR' + str(j+1) # Randomly sample the TIR dictionary for 100 iterations.
                try:
                    val = float(RNAduplexval(self.library[ASDname], TIRdict[TIRname]))
                except KeyError:
                    continue
                listofvals.append(val)

            average = sum(listofvals) / len(listofvals)
            self.library[average] = self.library.pop(ASDname)  # replaces the old key with the average
            # calculate the average binding energy between the ASD and all TIRs; here we store the key as the average
            # so that the we can call the names of the highest ASDs after the list is sorted

# --------------------------Code for GUI begins here-----------------------------------------------

def GUItoCode():
    """
    Gets inputs from the GUI and runs the code on them.
    :return: Puts the output of top 10 candidate ASDs on the GUI.
    """

    # --------- Obtain all variable inputs from user ----------
    CDSfileE = CDSfileEntry.get()
    genomeFileE = genomefileEntry.get()
    start = int(startrRNA.get())
    end = int(endrRNA.get())
    if start < end:
        direction = "forward"
    else:
        direction = "reverse"
    nameoforganism = nameoforgEntry.get()
    total_threads = int(total_threadsEntry.get())

    #---------- Call the algorithm using the user inputs ----------

    CDSfile = os.path.join(os.path.dirname(__file__), CDSfileE)
    genome = os.path.join(os.path.dirname(__file__), genomeFileE)
    path = os.path.dirname(__file__)
    rRNA = get16srRNAseq(genome, start, end, direction)
    initiallib = libbuild(rRNA)
    strain = Library(initiallib)  # create the initial organism library instance.
# -------------------------Step 2-3 ----------------------------
    thread_list = []
    for i in range(1, total_threads + 1):
        # Begin running each thread
        process = threading.Thread(target=strain.narrow_binding, args=(genome, start, end, direction, total_threads, i))
        process.start()
        #process.join()
        thread_list.append(process)

    [t.join() for t in thread_list]
    strain.deleter_fun()  # delete all the ASDs marked for deletion in step 2-3
    strain.updatelibindex()  #Update the indices of the library to be continuous
    thread_list = []
    # -----------------------Step 4-5-------------------------------
    for i in range(1, total_threads + 1):
        process = threading.Thread(target=strain.narrow_crossbinding,
                                   args=(genome, start, end, direction, total_threads, i))
        process.start()
        #process.join()
        thread_list.append(process)

    [t.join() for t in thread_list]

    strain.deleter_fun()
    strain.updatelibindex()

# -------------------------Step 6-7-----------------------------
    #thread_list = []
    # We do not run threads on step 7 because it's super inefficient to recalculate the secondary structure library for
    # each thread.
    # for i in range(1, total_threads + 1):
    #     process = threading.Thread(target=strain.ASD_2rystructure_narrow,
    #                                args=(genome, start, end, 'forward', total_threads, i))
    #     process.start()
    #     process.join()
    #     #thread_list.append(process)
    strain.ASD_2rystructure_narrow(genome, start, end, direction)

    #[t.join() for t in thread_list]
    strain.deleter_fun()
    strain.updatelibindex()

    pickle.dump(strain, open((nameoforganism + ".p"), 'wb'))

    # -------------------------Step 8-10-----------------------------
    thread_list = []
    TIRdict = getallTIRs(CDSfile, genome, nameoforganism, path)
    print(TIRdict)
    # we return dict of vals and listofaverages from the function then combine all the averages gained by the worker
    # threads to print the final answer
    for i in range(1, total_threads + 1):
        process = threading.Thread(target=strain.allASDTIRpairs, args=(TIRdict, total_threads, i))
        process.start()
        # process.join()
        # Get the string version of the dictionary, stored in the rawoutput box in the GUI, then convert to a
        # dictionary using ast
        thread_list.append(process)

    for t in thread_list:
        # dict = q.get()
        # completedict = completedict.update(dict)  # keeps adding the new dictionaries after every iteration
        t.join()
        outputbox.insert(END, "Done with thread " + str(thread_list.index(t)) + "\n")

    print(strain.library)
    listofaverages = list(strain.library.keys())
    # returns all the average RNAduplex values stored in the dictionary in a list
    print(listofaverages)
    listofaverages = listofaverages[int((len(listofaverages) / 2) + 1):]
    print(listofaverages)
    listofaverages.sort(reverse=True)
    outputbox.insert(END, "Here are the top ASD candidates:")
    for i in range(0, len(listofaverages)):
        outputbox.insert(END, str(i + 1) + ") " + str(strain.library[listofaverages[i]]) + "\n")


def time_keeper():
    for i in range(1, 10000000):
        progressbox.insert(END, i)
        time.sleep(2)


def process_manager():
    a = threading.Thread(target=GUItoCode)
    b = threading.Thread(target=time_keeper)
    a.start()
    b.start()

# ------------- Create GUI interface -----------------

from PIL import ImageDraw

root = Tk()
root.title("Orthogonal Ribosome Algorithm, V1.3, Rice University iGEM Team 2018")
image = Image.open(os.path.join(os.path.dirname(__file__), "g1878.png"))
width, height = image.size



draw = ImageDraw.Draw(image)

text_x = 55
text_y = 45
photoimage = ImageTk.PhotoImage(image)
Label(root, image=photoimage).place(x=420,y=150)
Label(root, text="Enter file name for CDS file (include .txt):").grid(row=0, column=0, sticky=W, padx=10, pady=10)
CDSfileEntry = Entry(root)
CDSfileEntry.grid(row=0, column=1, sticky=W, padx=10, pady=10)

Label(root, text="Enter file name for genome file (include .txt):").grid(row=1, column=0, sticky=W, padx=10, pady=10)
genomefileEntry = Entry(root)
genomefileEntry.grid(row=1, column=1, sticky=W, padx=10, pady=10)

Label(root, text="Enter the name of the organism:").grid(row=2, column=0, sticky=W, padx=10, pady=10)
nameoforgEntry = Entry(root)
nameoforgEntry.grid(row=2, column=1, sticky=W, padx=10, pady=10)

Label(root, text="Enter the start and end indices of the 16s rRNA:").grid(row=3, column=0, sticky=W, padx=10, pady=10)
Label(root, text="Start:").grid(row=4, column=0, sticky=W, padx=10, pady=10)
startrRNA = Entry(root)
startrRNA.grid(row=4, column=1, sticky=W, padx=10, pady=10)
Label(root, text="End:").grid(row=4, column=2, sticky=W, padx=10, pady=10)
endrRNA = Entry(root)
endrRNA.grid(row=4, column=3, sticky=W, padx=10, pady=10)

Label(root, text="Number of threads for parallel processing:").grid(row=5, column=2, sticky=W, padx=10, pady=10)
total_threadsEntry = Entry(root)
total_threadsEntry.grid(row=5, column=3, sticky=W, padx=10, pady=10)

Button(root, width=20, text="Start Algorithm", command=process_manager).grid(row=6, column=0, sticky=W, padx=10, pady=10)

Label(root, text='Here are the 10 top candidates with highest ASD-host binding values:').grid(row=7, column=0, sticky=W, padx=10, pady=10)
outputbox = Text(root, width=45, height=5, background="#819bc4")
outputbox.grid(row=8, column=0, sticky=W, padx=10, pady=10)

progressbox = Text(root, width=100, height=5, background="#819bb4")
progressbox.grid(row=10, column=0, sticky=W, padx=10, pady=10)
progressbox.grid_remove()







root.mainloop()
